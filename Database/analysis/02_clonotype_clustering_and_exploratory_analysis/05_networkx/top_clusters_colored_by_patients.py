# Network of top clusters shared across patients (tissue ignored)
# - Node colors represent subject_number
# - Pie nodes for sequences shared across multiple subjects
# - All nodes are circles

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import colorsys

# ============================ CONFIG =========================================
input_tsv   = '/doctorai/chiarba/analysis/MS_HC_naive_h0h1.tsv'
cluster_col = 'cluster_h1'
seq_col     = 'v_cdr3_j'
patient_col = 'subject_number'
study_col   = 'study'            # optional; if missing, treated as 'Unknown'
top_n       = 106
alphabet    = list("ACDEFGHIKLMNPQRSTVWY")  # AA alphabet for Hamming-1

# ===================== 1) Load & patient-only selection ======================
# NOTE: This block is exactly your requested selection logic
df = pd.read_csv(input_tsv, sep='\t')

# Identify clusters that contain at least one subject starting with 'H'
clusters_with_H = (
    df[df[patient_col].astype(str).str.startswith('H')][cluster_col].unique()
)

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
for col in [cluster_col, seq_col, patient_col, study_col]:
    if col in df_shared.columns:
        df_shared[col] = df_shared[col].astype(str).str.strip()

# Ensure required columns
required = {cluster_col, seq_col, patient_col}
missing = required - set(df_shared.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}. "
                     f"Need at least {sorted(required)}")

# ===================== Rank clusters and keep top-N ==========================
# Rank by (n_patients desc, n_sequences desc)
cluster_summary = (
    df_shared.groupby(cluster_col, as_index=False)
             .agg(n_patients=(patient_col, 'nunique'),
                  n_sequences=(seq_col, 'nunique'))
)
ranked = cluster_summary.sort_values(
    ['n_patients', 'n_sequences'], ascending=[False, False]
)
top_clusters = ranked.head(top_n)[cluster_col].tolist()
print("Top clusters (first 10 ids):", top_clusters[:10])

# Final working dataframe
dfw = df_shared[df_shared[cluster_col].isin(top_clusters)].copy()

# Ensure 'study' exists
if study_col not in dfw.columns:
    dfw[study_col] = 'Unknown'

# ===================== Build graph (within-cluster Hamming-1) ================
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

    # Add nodes with metadata + per-node subject distribution (for pies)
    for seq, cnt in seq_counts.items():
        rows_seq = sub[sub[seq_col] == seq]

        meta_subject = rows_seq[patient_col].iloc[0]
        meta_study   = rows_seq[study_col].iloc[0] if study_col in rows_seq else 'Unknown'

        subj_counts = rows_seq[patient_col].value_counts().to_dict()
        subj_studies = (
            rows_seq[[patient_col, study_col]]
            .drop_duplicates(subset=[patient_col])
            .set_index(patient_col)[study_col]
            .to_dict()
        )

        G.add_node(
            seq,
            count=int(cnt),
            cluster=cid,
            subject_counts=subj_counts,     # {subject -> count}
            subject_studies=subj_studies,   # {subject -> study}
            **{patient_col: meta_subject, study_col: meta_study}
        )

    G.add_edges_from(edges)

print("Graph built:", G.number_of_nodes(), "nodes,", G.number_of_edges(), "edges")

# Remove singletons for cleaner visualization
singles = [n for n, d in G.degree() if d == 0]
G.remove_nodes_from(singles)
print("After pruning singletons:", G.number_of_nodes(), "nodes,", G.number_of_edges(), "edges")
if G.number_of_nodes() == 0:
    raise RuntimeError("Graph is empty after pruning; relax filters or check input.")

# ===================== Node sizes ============================================
nodelist = list(G.nodes())
counts = np.array([G.nodes[n]['count'] for n in nodelist], dtype=float)
logc   = np.log1p(counts)
scaled = logc / logc.max() if logc.max() > 0 else np.zeros_like(logc)

# radius scaling for drawing patches
r_min, r_max = 0.4, 0.8

# ===================== Per-subject color palette =============================
NAMED_COLORS = {
    "Black":   "#000000", "White":  "#FFFFFF", "Red":   "#FF0000", "Lime":  "#00FF2F",
    "Blue":    "#0202FF", "Yellow": "#FFFF00", "Cyan":  "#00FFFF", "Magenta": "#FF00A2",
    "Silver":  "#C0C0C0", "Gray":   "#808080", "Maroon":"#800000", "Olive":   "#808000",
    "Green":   "#129A00", "Purple": "#800069", "Teal":  "#008080", "Navy":    "#3F2156",
}
STUDY_COLORNAME = {
    "Greenfield": "Blue", "Stern": "Magenta", "Ramesh": "Gray", "Palanichemy": "Green",
    "Agrafiotis": "Gray", "Laurent": "Maroon", "lomakin": "Teal", "Perez_Salvidar": "Yellow",
    "Ruschill": "Navy", "Unknown": "Gray",
}
STUDY_BASE_COLORS = {
    study: NAMED_COLORS.get(name, NAMED_COLORS["Gray"])
    for study, name in STUDY_COLORNAME.items()
}

def shades_from_hex(hexstr, n, min_l=0.20, max_l=0.90):
    """Keep hue & saturation, vary lightness in HLS space to get n shades."""
    r, g, b = mcolors.to_rgb(hexstr)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    if n <= 1:
        ls = [0.60]
    else:
        ls = [min_l + (max_l - min_l) * i / (n - 1) for i in range(n)]
    return [colorsys.hls_to_rgb(h, li, s) for li in ls]

# subjects grouped by study to generate distinct shades
subjects_per_study = {}
all_subjects = set()
for n in nodelist:
    sc = G.nodes[n].get('subject_counts')
    ss = G.nodes[n].get('subject_studies')
    if sc and ss:
        for subj in sc:
            study = ss.get(subj, G.nodes[n].get(study_col, 'Unknown'))
            study = 'Unknown' if pd.isna(study) else study
            subjects_per_study.setdefault(study, set()).add(str(subj))
            all_subjects.add(str(subj))
    else:
        subj  = str(G.nodes[n].get(patient_col))
        study = G.nodes[n].get(study_col, 'Unknown')
        study = 'Unknown' if pd.isna(study) else study
        subjects_per_study.setdefault(study, set()).add(subj)
        all_subjects.add(subj)

group_color_map = {}
for study, subs in subjects_per_study.items():
    subs = sorted(subs, key=str)
    base_hex = STUDY_BASE_COLORS.get(study, NAMED_COLORS["Gray"])
    shades = shades_from_hex(base_hex, len(subs))
    for subj, col in zip(subs, shades):
        group_color_map[subj] = col

# representative color selection for single-subject nodes
rep_subjects = []
for n in nodelist:
    sc = G.nodes[n].get('subject_counts')
    if sc and len(sc) == 1:
        rep_subjects.append(next(iter(sc.keys())))
    else:
        rep_subjects.append(G.nodes[n].get(patient_col))
node_colors = [group_color_map.get(str(s), 'lightgray') for s in rep_subjects]

# ===================== Layout (per-cluster grid) =============================
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
        subG, k=INTERNAL_K, iterations=INTERNAL_ITER, scale=INTERNAL_SCALE, seed=42
    )
    r, c = divmod(idx, cols)
    cx, cy = c * GRID_SPACING, -r * GRID_SPACING
    for node, (x, y) in subpos.items():
        pos[node] = (cx + x, cy + y)

# fill any missing positions
missing = [n for n in G.nodes() if n not in pos]
if missing:
    pos.update(nx.spring_layout(G.subgraph(missing), seed=123, k=2.0,
                                iterations=200, scale=INTERNAL_SCALE))

# ===================== Draw ===================================================
fig, ax = plt.subplots(figsize=(12, 12))

# edges
nx.draw_networkx_edges(
    G, pos, edgelist=list(G.edges()), edge_color='black', alpha=0.28, width=0.8, ax=ax
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
        wedge = mpatches.Wedge(center, radius, theta1, theta2,
                               facecolor=c, edgecolor='k', linewidth=0.25, zorder=zorder)
        ax.add_patch(wedge)
        angle = theta2
    # crisp border
    outline = mpatches.Circle(center, radius, facecolor='none', edgecolor='k', linewidth=0.6, zorder=zorder+0.5)
    ax.add_patch(outline)

idx = {n:i for i,n in enumerate(nodelist)}

for n in nodelist:
    i = idx[n]
    r = float(r_min + (r_max - r_min) * scaled[i])
    center = pos[n]

    sc = G.nodes[n].get('subject_counts', {})
    is_pie = (sc is not None) and (len(sc) > 1)

    if is_pie:
        subs = sorted(sc.keys(), key=str)
        vals = [sc[s] for s in subs]
        cols = [group_color_map.get(str(s), 'lightgray') for s in subs]
        draw_pie_node(ax, center, vals, cols, radius=r, zorder=3)
    else:
        node_color = node_colors[i]
        circ = mpatches.Circle(center, r, facecolor=node_color, edgecolor='k', linewidth=0.6, alpha=0.95, zorder=3)
        ax.add_patch(circ)

# axis cosmetics
xs = np.array([xy[0] for xy in pos.values()], dtype=float)
ys = np.array([xy[1] for xy in pos.values()], dtype=float)
pad = INTERNAL_SCALE * 0.8
ax.set_xlim(xs.min() - pad, xs.max() + pad)
ax.set_ylim(ys.min() - pad, ys.max() + pad)
ax.set_aspect('equal', adjustable='box')
plt.axis('off')
plt.tight_layout()
plt.show()


# --- Legend for node coloring (subjects) ---
from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(10, 2))

handles, labels = [], []
for subj in sorted(all_subjects, key=str):
    handles.append(
        Line2D([], [], marker='o', linestyle='',
               markerfacecolor=group_color_map.get(str(subj), 'lightgray'),
               markeredgecolor='k', markersize=6, alpha=0.95)
    )
    labels.append(str(subj))

legend = ax.legend(
    handles, labels,
    title="Node colors represent individual subjects (pie = shared across subjects)",
    loc='center', fontsize=7, title_fontsize=8,
    borderpad=0.3, labelspacing=0.3, handletextpad=0.5,
    ncol=8, frameon=False
)

ax.axis("off")
plt.tight_layout()
plt.show()
