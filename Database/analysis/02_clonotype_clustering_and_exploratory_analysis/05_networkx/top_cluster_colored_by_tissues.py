# Network of top clusters shared across patients (nodes colored by tissue)
# - Node colors represent tissue
# - Pie nodes for sequences shared across multiple tissues
# - All nodes are circles

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

# ============================ CONFIG =========================================
input_tsv   = '/doctorai/chiarba/analysis/MS_HC_naive_h0h1.tsv'
cluster_col = 'cluster_h1'
seq_col     = 'v_cdr3_j'
patient_col = 'subject_number'
study_col   = 'study'      # optional; if missing, treated as 'Unknown'
tissue_col  = 'tissue'     # column with tissue information
top_n       = 106
alphabet    = list("ACDEFGHIKLMNPQRSTVWY")  # AA alphabet for Hamming-1

# ===================== Fixed tissue universe and colors ======================
TISSUE_ORDER = [
    "blood",
    "CSF",
    "Brain lesion",
    "Pia mater",
    "Choroid plexus",
    "CLN",
    "Unknown"
]

TISSUE_COLORS = {
    "blood":           "#d62728",  # red
    "CSF":             "#2089e4",  # blue
    "Brain lesion":    "#9467bd",  # purple
    "Pia mater":       "#2ca02c",  # green
    "Choroid plexus":  "#ff7f0e",  # orange
    "CLN":             "#17becf",  # cyan
    "Unknown":         "#bdbdbd",  # gray fallback
}

DEFAULT_TISSUE_COLOR = TISSUE_COLORS["Unknown"]

# ===================== 1) Load & patient-only selection ======================
# ATTENTION: selection of patients with removal of HC (subjects starting with "H")

df = pd.read_csv(input_tsv, sep='\t')

# Identify clusters that contain at least one subject starting with 'H'
clusters_with_H = df[df[patient_col].astype(str).str.startswith('H')][cluster_col].unique()

# Remove entire clusters if they contain ANY 'H' subject
df_noH = df[~df[cluster_col].isin(clusters_with_H)]

# Keep clusters shared by >1 subject_number (now all non-H by definition)
shared_clusters = (
    df_noH.groupby(cluster_col)[patient_col]
          .nunique()
          .loc[lambda x: x > 1]
).index
df_shared = df_noH[df_noH[cluster_col].isin(shared_clusters)].copy()

# Light normalization (strings, strip spaces)
for col in [cluster_col, seq_col, patient_col, study_col, tissue_col]:
    if col in df_shared.columns:
        df_shared[col] = df_shared[col].astype(str).str.strip()

# Ensure required columns
required = {cluster_col, seq_col, patient_col, tissue_col}
missing = required - set(df_shared.columns)
if missing:
    raise ValueError(
        f"Missing required columns: {missing}. "
        f"Need at least {sorted(required)}"
    )

# Fill optional columns
if study_col not in df_shared.columns:
    df_shared[study_col] = 'Unknown'
if tissue_col not in df_shared.columns:
    df_shared[tissue_col] = 'Unknown'
df_shared[tissue_col] = (
    df_shared[tissue_col]
    .replace({'': 'Unknown'})
    .fillna('Unknown')
)

# ===================== 2) Rank clusters and keep top-N =======================
cluster_summary = (
    df_shared.groupby(cluster_col, as_index=False)
             .agg(n_patients=(patient_col, 'nunique'),
                  n_sequences=(seq_col, 'nunique'))
)

ranked = cluster_summary.sort_values(
    ['n_patients', 'n_sequences'],
    ascending=[False, False]
)

top_clusters = ranked.head(top_n)[cluster_col].tolist()
print("Top clusters (first 10 ids):", top_clusters[:10])

# Final working dataframe
dfw = df_shared[df_shared[cluster_col].isin(top_clusters)].copy()

# ===================== 3) Build graph (within-cluster Hamming-1) =============
G = nx.Graph()

for cid, sub in dfw.groupby(cluster_col):
    # Per-cluster sequence counts
    seq_counts = Counter(sub[seq_col])

    # Build Hamming-1 edges within this cluster (only compare sequences of equal length)
    edges = set()
    by_len = {}
    for s in seq_counts:
        by_len.setdefault(len(s), []).append(s)

    for length, seqs in by_len.items():
        seqset = set(seqs)
        for seq in seqs:
            for pos in range(length):
                base = seq[:pos] + "{}" + seq[pos+1:]
                for aa in alphabet:
                    if aa == seq[pos]:
                        continue
                    nbr = base.format(aa)
                    if nbr in seqset and seq < nbr:  # avoid duplicates
                        edges.add((seq, nbr))

    # Add nodes with metadata + per-node tissue distribution (for pies)
    for seq, cnt in seq_counts.items():
        rows_seq = sub[sub[seq_col] == seq].copy()

        # Clean / standardize tissue values
        rows_seq[tissue_col] = (
            rows_seq[tissue_col]
            .replace({'': 'Unknown'})
            .fillna('Unknown')
        )
        rows_seq[tissue_col] = rows_seq[tissue_col].where(
            rows_seq[tissue_col].isin(TISSUE_COLORS.keys()),
            'Unknown'
        )

        # tissue frequency for this sequence (node)
        tissue_counts = rows_seq[tissue_col].value_counts().to_dict()

        # Store a representative tissue for single-tissue nodes
        rep_tissue = next(iter(tissue_counts)) if len(tissue_counts) == 1 else None

        G.add_node(
            seq,
            count=int(cnt),
            cluster=cid,
            tissue_counts=tissue_counts,  # {tissue -> count}
            rep_tissue=rep_tissue,        # None if multi-tissue
        )

    G.add_edges_from(edges)

print("Graph built:", G.number_of_nodes(), "nodes,", G.number_of_edges(), "edges")

# Remove singletons for cleaner visualization
singles = [n for n, d in G.degree() if d == 0]
G.remove_nodes_from(singles)
print(
    "After pruning singletons:",
    G.number_of_nodes(), "nodes,",
    G.number_of_edges(), "edges"
)
if G.number_of_nodes() == 0:
    raise RuntimeError("Graph is empty after pruning; relax filters or check input.")

# ===================== 4) Node sizes =========================================
nodelist = list(G.nodes())
counts = np.array([G.nodes[n]['count'] for n in nodelist], dtype=float)
logc   = np.log1p(counts)
scaled = logc / logc.max() if logc.max() > 0 else np.zeros_like(logc)

# radius scaling for drawing patches
r_min, r_max = 0.4, 0.8

# ===================== 5) Colors by tissue ===================================
def tissue_color(name: str) -> str:
    """Return hex color for a tissue name with fallback to gray."""
    return TISSUE_COLORS.get(name, DEFAULT_TISSUE_COLOR)

# Prepare representative color for single-tissue nodes
node_colors = []
for n in nodelist:
    rep = G.nodes[n].get('rep_tissue')
    node_colors.append(tissue_color(rep if rep else "Unknown"))

# ===================== 6) Layout (per-cluster grid) ==========================
INTERNAL_K     = 3.0
INTERNAL_ITER  = 200
INTERNAL_SCALE = 5.0
GRID_SPACING   = 12.0

nodes_by_cluster = defaultdict(list)
for n, data in G.nodes(data=True):
    nodes_by_cluster[data.get('cluster')].append(n)

cluster_ids = list(nodes_by_cluster.keys())
cols = int(np.ceil(np.sqrt(len(cluster_ids))))
pos  = {}

for idx, cid in enumerate(cluster_ids):
    nodes_c = nodes_by_cluster[cid]
    if not nodes_c:
        continue
    subG = G.subgraph(nodes_c)
    subpos = nx.spring_layout(
        subG,
        k=INTERNAL_K,
        iterations=INTERNAL_ITER,
        scale=INTERNAL_SCALE,
        seed=42
    )
    r, c = divmod(idx, cols)
    cx, cy = c * GRID_SPACING, -r * GRID_SPACING
    for node, (x, y) in subpos.items():
        pos[node] = (cx + x, cy + y)

# fill any missing positions
missing = [n for n in G.nodes() if n not in pos]
if missing:
    pos.update(
        nx.spring_layout(
            G.subgraph(missing),
            seed=123,
            k=2.0,
            iterations=200,
            scale=INTERNAL_SCALE
        )
    )

# ===================== 7) Draw network =======================================
fig, ax = plt.subplots(figsize=(12, 12))

# edges
nx.draw_networkx_edges(
    G,
    pos,
    edgelist=list(G.edges()),
    edge_color='black',
    alpha=0.28,
    width=0.8,
    ax=ax
)

def draw_pie_node(ax, center, values, colors, radius, zorder=3):
    """Draw a circular pie at center with given values/colors."""
    values = np.asarray(values, dtype=float)
    total = values.sum()
    if total <= 0:
        return
    angle = 90.0
    for v, c in zip(values, colors):
        theta1 = angle
        theta2 = angle + (v / total) * 360.0
        wedge = mpatches.Wedge(
            center,
            radius,
            theta1,
            theta2,
            facecolor=c,
            edgecolor='k',
            linewidth=0.25,
            zorder=zorder
        )
        ax.add_patch(wedge)
        angle = theta2
    # crisp border
    outline = mpatches.Circle(
        center,
        radius,
        facecolor='none',
        edgecolor='k',
        linewidth=0.6,
        zorder=zorder + 0.5
    )
    ax.add_patch(outline)

idx_map = {n: i for i, n in enumerate(nodelist)}

for n in nodelist:
    i = idx_map[n]
    r = float(r_min + (r_max - r_min) * scaled[i])
    center = pos[n]

    tc = G.nodes[n].get('tissue_counts', {})
    is_pie = (tc is not None) and (len(tc) > 1)

    if is_pie:
        # sort wedges by predefined TISSUE_ORDER for consistent order
        tissues_present = [t for t in TISSUE_ORDER if t in tc] + [
            t for t in sorted(tc.keys()) if t not in TISSUE_ORDER
        ]
        vals = [tc[t] for t in tissues_present]
        cols = [tissue_color(t) for t in tissues_present]
        draw_pie_node(ax, center, vals, cols, radius=r, zorder=3)
    else:
        # single-tissue node
        node_color = node_colors[i]
        circ = mpatches.Circle(
            center,
            r,
            facecolor=node_color,
            edgecolor='k',
            linewidth=0.6,
            alpha=0.95,
            zorder=3
        )
        ax.add_patch(circ)

# axes cosmetics
xs = np.array([xy[0] for xy in pos.values()], dtype=float)
ys = np.array([xy[1] for xy in pos.values()], dtype=float)
pad = INTERNAL_SCALE * 0.8
ax.set_xlim(xs.min() - pad, xs.max() + pad)
ax.set_ylim(ys.min() - pad, ys.max() + pad)
ax.set_aspect('equal', adjustable='box')
plt.axis('off')
plt.tight_layout()
plt.show()

# Horizontal Separate legend 

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# Fallbacks in case variables are not in scope
if 'TISSUE_ORDER' not in globals():
    TISSUE_ORDER = []
if 'TISSUE_COLORS' not in globals():
    TISSUE_COLORS = {}
fallback_color = '#b0b0b0'  # default color for unknown tissues

# Collect tissues present in the graph
present_tissues = set()
for n in G.nodes():
    tcounts = G.nodes[n].get('tissue_counts', {})
    present_tissues.update(tcounts.keys())

# Order: preferred order first, then alphabetical for the rest
labels = [t for t in TISSUE_ORDER if t in present_tissues] + \
         [t for t in sorted(present_tissues) if t not in TISSUE_ORDER]

# Build compact marker handles (circles) to avoid text overlap
handles = [
    Line2D([0], [0],
           marker='o', linestyle='None',
           markerfacecolor=TISSUE_COLORS.get(t, fallback_color),
           markeredgecolor='k',
           markersize=9)   # tweak size if needed
    for t in labels
]

# Figure size scales with number of items to prevent crowding
fig_w = max(6, 1.0 * max(1, len(labels)))  # increase multiplier for more spacing
fig_h = 1.2
fig_leg, ax_leg = plt.subplots(figsize=(fig_w, fig_h))
ax_leg.axis('off')

leg = ax_leg.legend(
    handles=handles,
    labels=labels,
    title="Tissue",
    ncol=len(labels) if labels else 1,  # single row
    loc='center',
    mode='expand',
    frameon=True,
    borderaxespad=0.3,
    handlelength=1.0,     # handle spacing
    handletextpad=0.5,    # space between marker and text
    columnspacing=1.0,    # space between legend entries
    fontsize=10,
    title_fontsize=10,
)

if leg:
    leg.get_frame().set_alpha(0.9)

plt.tight_layout()
plt.show()

# Optional: save to file
# fig_leg.savefig("legend_horizontal_compact.png", dpi=300, bbox_inches="tight")
