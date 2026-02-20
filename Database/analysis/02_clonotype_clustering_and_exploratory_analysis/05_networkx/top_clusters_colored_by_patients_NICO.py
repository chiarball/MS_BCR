# Network of top clusters shared across patients (tissue ignored)
# - Node colors represent subject_number
# - Pie nodes for sequences shared across multiple subjects
# - All nodes are circles

# Choose scoring method: "A" or "B"
score_method = "B"   # options: "A" (size-weighted) or "B" (distribution-weighted)
score_method = "A" 
score_method = "C"
score_method = "G"
score_method = "E"


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
tissue_col  = 'tissue'
top_n       = 106
alphabet    = list("ACDEFGHIKLMNPQRSTVWY")  # AA alphabet for Hamming-1
base_folder = '/doctorai/niccoloc/MS_db/MS_BCR'



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
    "blood":           "#d62728",  # red
    "Unknown":         "#bdbdbd"   # gray fallback
}

DEFAULT_TISSUE_COLOR = TISSUE_COLORS["Unknown"]

def plot_network(dfw, cluster_col, seq_col, patient_col, study_col, tissue_col, alphabet, base_folder, score_method, color_by='patients'):
    # ===================== Build graph (within-cluster Hamming-1) ================
    from matplotlib.lines import Line2D

    G = nx.Graph()
    cluster_avg_pgen = {}  # New dict to store avg_pgen per cluster

    for cid, sub in dfw.groupby(cluster_col):
        # Per-cluster sequence counts
        seq_counts = Counter(sub[seq_col])
        # Compute avg_pgen excluding NaNs, without dropping rows
        if 'Pgen' in sub.columns:
            avg_pgen = sub['Pgen'].mean(skipna=True)
        else:
            avg_pgen = None
        cluster_avg_pgen[cid] = avg_pgen  # Store avg_pgen for this cluster

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

        # Add nodes with metadata + per-node subject and tissue distribution (for pies)
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

            meta_subject = rows_seq[patient_col].iloc[0]
            meta_study   = rows_seq[study_col].iloc[0] if study_col in rows_seq else 'Unknown'

            subj_counts = rows_seq[patient_col].value_counts().to_dict()
            subj_studies = (
                rows_seq[[patient_col, study_col]]
                .drop_duplicates(subset=[patient_col])
                .set_index(patient_col)[study_col]
                .to_dict()
            )

            # tissue frequency for this sequence (node)
            tissue_counts = rows_seq[tissue_col].value_counts().to_dict()

            # Store a representative tissue for single-tissue nodes
            rep_tissue = next(iter(tissue_counts)) if len(tissue_counts) == 1 else None

            G.add_node(
                seq,
                count=int(cnt),
                avg_pgen=avg_pgen,
                cluster=cid,
                subject_counts=subj_counts,     # {subject -> count}
                subject_studies=subj_studies,   # {subject -> study}
                tissue_counts=tissue_counts,    # {tissue -> count}
                rep_tissue=rep_tissue,          # None if multi-tissue
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

    if color_by == 'patients':
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

    elif color_by == 'tissues':
        def tissue_color(name: str) -> str:
            """Return hex color for a tissue name with fallback to gray."""
            return TISSUE_COLORS.get(name, DEFAULT_TISSUE_COLOR)

        # Prepare representative color for single-tissue nodes
        node_colors = []
        for n in nodelist:
            rep = G.nodes[n].get('rep_tissue')
            node_colors.append(tissue_color(rep if rep else "Unknown"))

    else:
        raise ValueError("color_by must be 'patients' or 'tissues'")

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
    #print the first 10 cluster ids
    print("Cluster IDs:", cluster_ids[:10])
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

    # Compute cluster centroids for labeling
    cluster_centers = {}
    for cid, nodes_c in nodes_by_cluster.items():
        if not nodes_c:
            continue
        cluster_positions = [pos[n] for n in nodes_c if n in pos]
        if cluster_positions:
            avg_x = sum(x for x, y in cluster_positions) / len(cluster_positions)
            avg_y = sum(y for x, y in cluster_positions) / len(cluster_positions)
            cluster_centers[cid] = (avg_x, avg_y)

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

        if color_by == 'patients':
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

        elif color_by == 'tissues':
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
                circ = mpatches.Circle(center, r, facecolor=node_color, edgecolor='k', linewidth=0.6, alpha=0.95, zorder=3)
                ax.add_patch(circ)

    # After drawing nodes, add text labels for average Pgen at cluster centroids
    for cid, center in cluster_centers.items():
        avg_pgen = cluster_avg_pgen.get(cid)
        if avg_pgen is not None:
            ax.text(center[0], center[1] + 6.0, f"{avg_pgen:.2f}", ha='center', va='center', fontsize=7, color='grey')

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
    plt.savefig(f"{base_folder}/Network_clones_method_{score_method}_{color_by.upper()}.png", dpi=300, bbox_inches="tight")

    if color_by == 'patients':
        # --- Legend for node coloring (subjects) ---
        fig, ax = plt.subplots(figsize=(10, 2)) 
        # Horizontal Separate legend 

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
            title="Score method: " + score_method,
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

        if legend:
            legend.get_frame().set_alpha(0.9)

        plt.tight_layout()
        plt.show()
        fig.savefig(f"{base_folder}/Network_clones_method_{score_method}_{color_by.upper()}_legend.png", dpi=300, bbox_inches="tight")

    elif color_by == 'tissues':
        # Horizontal Separate legend 

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
                   markerfacecolor=tissue_color(t),
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
            title="Score method: " + score_method,
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
        fig_leg.savefig(f"{base_folder}/Network_clones_method_{score_method}_{color_by.upper()}_legend.png", dpi=300, bbox_inches="tight")


        def plot_network_legend_cols(dfw, cluster_col, seq_col, patient_col, study_col, tissue_col, alphabet, base_folder, score_method, color_by='patients', column_width = 0.10, save_svg=False):
            # ===================== Build graph (within-cluster Hamming-1) ================
            from matplotlib.lines import Line2D

            G = nx.Graph()
            cluster_avg_pgen = {}  # New dict to store avg_pgen per cluster

            for cid, sub in dfw.groupby(cluster_col):
                # Per-cluster sequence counts
                seq_counts = Counter(sub[seq_col])
                # Compute avg_pgen excluding NaNs, without dropping rows
                if 'Pgen' in sub.columns:
                    avg_pgen = sub['Pgen'].mean(skipna=True)
                else:
                    avg_pgen = None
                cluster_avg_pgen[cid] = avg_pgen  # Store avg_pgen for this cluster

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

                # Add nodes with metadata + per-node subject and tissue distribution (for pies)
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

                    meta_subject = rows_seq[patient_col].iloc[0]
                    meta_study   = rows_seq[study_col].iloc[0] if study_col in rows_seq else 'Unknown'

                    subj_counts = rows_seq[patient_col].value_counts().to_dict()
                    subj_studies = (
                        rows_seq[[patient_col, study_col]]
                        .drop_duplicates(subset=[patient_col])
                        .set_index(patient_col)[study_col]
                        .to_dict()
                    )

                    # tissue frequency for this sequence (node)
                    tissue_counts = rows_seq[tissue_col].value_counts().to_dict()

                    # Store a representative tissue for single-tissue nodes
                    rep_tissue = next(iter(tissue_counts)) if len(tissue_counts) == 1 else None

                    G.add_node(
                        seq,
                        count=int(cnt),
                        avg_pgen=avg_pgen,
                        cluster=cid,
                        subject_counts=subj_counts,     # {subject -> count}
                        subject_studies=subj_studies,   # {subject -> study}
                        tissue_counts=tissue_counts,    # {tissue -> count}
                        rep_tissue=rep_tissue,          # None if multi-tissue
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

            if color_by == 'patients':
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

            elif color_by == 'tissues':
                def tissue_color(name: str) -> str:
                    """Return hex color for a tissue name with fallback to gray."""
                    return TISSUE_COLORS.get(name, DEFAULT_TISSUE_COLOR)

                # Prepare representative color for single-tissue nodes
                node_colors = []
                for n in nodelist:
                    rep = G.nodes[n].get('rep_tissue')
                    node_colors.append(tissue_color(rep if rep else "Unknown"))

            else:
                raise ValueError("color_by must be 'patients' or 'tissues'")

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
            #print the first 10 cluster ids
            print("Cluster IDs:", cluster_ids[:10])
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

            # Compute cluster centroids for labeling
            cluster_centers = {}
            for cid, nodes_c in nodes_by_cluster.items():
                if not nodes_c:
                    continue
                cluster_positions = [pos[n] for n in nodes_c if n in pos]
                if cluster_positions:
                    avg_x = sum(x for x, y in cluster_positions) / len(cluster_positions)
                    avg_y = sum(y for x, y in cluster_positions) / len(cluster_positions)
                    cluster_centers[cid] = (avg_x, avg_y)

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

                if color_by == 'patients':
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

                elif color_by == 'tissues':
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
                        circ = mpatches.Circle(center, r, facecolor=node_color, edgecolor='k', linewidth=0.6, alpha=0.95, zorder=3)
                        ax.add_patch(circ)

            # After drawing nodes, add text labels for average Pgen at cluster centroids
            for cid, center in cluster_centers.items():
                avg_pgen = cluster_avg_pgen.get(cid)
                if avg_pgen is not None:
                    ax.text(center[0], center[1] + 6.0, f"{avg_pgen:.2f}", ha='center', va='center', fontsize=7, color='grey')

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
            file_png = f"{base_folder}/Network_clones_method_{score_method}_{color_by.upper()}.png"
            fig.savefig(file_png, dpi=300, bbox_inches="tight")
            if save_svg:
                file_svg = f"{base_folder}/Network_clones_method_{score_method}_{color_by.upper()}.svg"
                fig.savefig(file_svg, format='svg', bbox_inches="tight")

            if color_by == 'patients':
                studies = sorted(subjects_per_study.keys(), key=str)
                n = len(studies)

                fig_w = 1.6 * n
                fig_h = 4.0
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                #fig title "Patient Legend"
                fig.suptitle("Patient Legend", fontsize=10)
                ax.axis('off')

                # x positions: centers of equal-width columns
                n = len(studies)
                total_width = n * column_width
                left = (1.0 - total_width) / 2.0   # center all columns

                x_positions = [
                    left + column_width * (i + 0.5)
                    for i in range(n)
                ]

                for x, study in zip(x_positions, studies):
                    subs = sorted(subjects_per_study[study], key=str)

                    handles = [
                        Line2D(
                            [], [], marker='o', linestyle='',
                            markerfacecolor=group_color_map[sub],
                            markeredgecolor='k',
                            markersize=7
                        )
                        for sub in subs
                    ]

                    leg = ax.legend(
                        handles,
                        subs,
                        loc='upper center',
                        bbox_to_anchor=(x, 0.98),
                        frameon=False,          # as requested
                        handletextpad=0.5,
                        labelspacing=0.4,
                        fontsize=10,
                    )

                    ax.add_artist(leg)

                plt.tight_layout(pad=0.3)
                plt.show()

                patient_legend_png = f"{base_folder}/Network_clones_method_{score_method}_PATIENTS_legend.png"
                fig.savefig(patient_legend_png, dpi=300, bbox_inches="tight")
                if save_svg:
                    patient_legend_svg = f"{base_folder}/Network_clones_method_{score_method}_PATIENTS_legend.svg"
                    fig.savefig(patient_legend_svg, format='svg', bbox_inches="tight")
            elif color_by == 'tissues':
                # Horizontal Separate legend 

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
                           markerfacecolor=tissue_color(t),
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
                    title="Tissues",
                    ncol=len(labels) if labels else 1,  # single row
                    loc='center',
                    mode='expand',
                    frameon=False,
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
                tissue_legend_png = f"{base_folder}/Network_clones_method_{score_method}_{color_by.upper()}_legend.png"
                fig_leg.savefig(tissue_legend_png, dpi=300, bbox_inches="tight")
                if save_svg:
                    tissue_legend_svg = f"{base_folder}/Network_clones_method_{score_method}_{color_by.upper()}_legend.svg"
                    fig_leg.savefig(tissue_legend_svg, format='svg', bbox_inches="tight")


def compute_pgen_with_sonnia(df,
                             seq_col='junction_aa',
                             vgene_col='v_gene',
                             jgene_col='j_gene',
                             pgen_model='humanIGH',
                             max_length=40):
    """
    Validate input columns, run SoNNia Pgen evaluation and return the filtered dataframe
    with an added 'Pgen' column. Raises ValueError if required columns are missing or
    end up not in the expected order after renaming.
    """
    import pandas as pd
    from sonnia.sonnia import SoNNia
    from sonnia.processing import Processing

    # check presence of required columns
    required = [seq_col, vgene_col, jgene_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # select and rename to the expected order/names
    input_df = df[required].copy().drop_duplicates()
    input_df.columns = ['amino_acid', 'v_gene', 'j_gene']

    # enforce column order / names
    if list(input_df.columns) != ['amino_acid', 'v_gene', 'j_gene']:
        raise ValueError("Input dataframe must have columns in order: "
                         "['amino_acid', 'v_gene', 'j_gene'] after renaming.")

    # initialize SoNNia and processor
    qm = SoNNia(pgen_model=pgen_model)
    processor = Processing(pgen_model=pgen_model, max_length=max_length)

    # filter dataframe using sonnia's processor
    filtered = processor.filter_dataframe(input_df)
    if filtered is None or len(filtered) == 0:
        # return empty dataframe with expected columns + Pgen
        print("No valid sequences after filtering.")
        return pd.DataFrame(columns=df.columns.tolist() + ['Pgen'])

    # prepare array for evaluation (matching prior usage)
    arr = filtered.values.astype(str)

    # evaluate sequences: get pgen_data and attach to filtered dataframe
    Q_data, pgen_data, ppost_data = qm.evaluate_seqs(arr)

    filtered = filtered.reset_index(drop=True).copy()
    filtered['Pgen'] = pgen_data

    return filtered


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
cluster_summary = (
    df_shared.groupby(cluster_col, as_index=False)
             .agg(n_patients=(patient_col, 'nunique'),
                  n_sequences=(seq_col, 'nunique'))
)
# Sharedness
sharedness = (
    df_shared.groupby(cluster_col)[patient_col]
             .nunique()
             .div(df_shared[patient_col].nunique())
             .rename("sharedness")
)

# Dissemination
dissemination = (
    df_shared.groupby(cluster_col)[tissue_col]
             .nunique()
             .div(df_shared[tissue_col].nunique())
             .rename("dissemination")
)

# Cluster size
cluster_size = (
    df_shared.groupby(cluster_col)[seq_col]
             .nunique()
             .rename("cluster_size")
)

# Log-scaled size weight
size_weight = (
    np.log1p(cluster_size) /
    np.log1p(cluster_size.max())
).rename("size_weight")

# Log-scaled size weight
size_weight = (
    (cluster_size) /
    (cluster_size.max())
).rename("size_weight")

# sequences per (cluster, patient)
seq_per_patient = (
    df_shared.groupby([cluster_col, patient_col])[seq_col]
             .nunique()
)

# total sequences per cluster
total_seq = seq_per_patient.groupby(level=0).sum()

# fraction per patient, summed
weighted_sharedness = (
    (seq_per_patient / total_seq)
    .groupby(level=0)
    .sum()
    .rename("weighted_sharedness")
)



def pielou_evenness(counts):
    counts = np.asarray(counts)
    counts = counts[counts > 0]
    if len(counts) <= 1:
        return 0.0
    p = counts / counts.sum()
    H = -np.sum(p * np.log(p))
    return H / np.log(len(counts))

patient_evenness = (
    df_shared
    .groupby(cluster_col)
    .apply(lambda x: pielou_evenness(
        x[patient_col].value_counts().values
    ))
    .rename("patient_evenness")
)

tissue_evenness = (
    df_shared
    .groupby(cluster_col)
    .apply(lambda x: pielou_evenness(
        x[tissue_col].value_counts().values
    ))
    .rename("tissue_evenness")
)



def gini_evenness(counts):
    counts = np.asarray(counts)
    if counts.sum() == 0:
        return 0.0
    sorted_counts = np.sort(counts)
    n = len(counts)
    cum = np.cumsum(sorted_counts)
    gini = (n + 1 - 2 * np.sum(cum) / cum[-1]) / n
    return 1 - gini

import numpy as np

sharedness_gini = (
    df_shared
    .groupby(cluster_col)
    .apply(lambda x: gini_evenness(
        x[patient_col].value_counts().values
    ))
    .rename("sharedness_gini")
)

dissemination_gini = (
    df_shared
    .groupby(cluster_col)
    .apply(lambda x: gini_evenness(
        x[tissue_col].value_counts().values
    ))
    .rename("dissemination_gini")
)

alfa = 0.40
beta = 0.40
gamma = 0.2

# Final score
score_A = (
    sharedness *
    dissemination *
    size_weight
).rename("score_size_weighted")

score_B = (
    weighted_sharedness *
    dissemination
).rename("score_distribution_weighted")

score_C = (
    sharedness *
    dissemination
).rename("score_unweighted")

score_D = (
    (sharedness * patient_evenness) *
    (dissemination * tissue_evenness)*
    size_weight
).rename("score_evenness_weighted")

score_E = (
    (alfa * sharedness * patient_evenness) *
    (beta * dissemination * tissue_evenness)*
    (gamma * cluster_size)
).rename("score_evenness_weighted_abg")

score_G = (
    (alfa * sharedness * sharedness_gini) *
    (beta * dissemination * dissemination_gini) *
    (gamma * cluster_size)
).rename("score_gini_evenness_weighted")


cluster_metrics = (
    cluster_summary
    .set_index(cluster_col)
    .join([
        sharedness,
        dissemination,
        cluster_size,
        size_weight,
        weighted_sharedness,
        patient_evenness,
        tissue_evenness,
        sharedness_gini,
        dissemination_gini,
        score_A,
        score_B,
        score_C,
        score_D,
        score_G,
        score_E
    ])
    .reset_index()
)


if score_method == "A":
    score_col = "score_size_weighted"
elif score_method == "B":
    score_col = "score_distribution_weighted"
elif score_method == "C":
    score_col = "score_unweighted"
elif score_method == "D":
    score_col = "score_evenness_weighted"
elif score_method == "G":
    score_col = "score_gini_evenness_weighted"
elif score_method == "E":
    score_col = "score_evenness_weighted_abg"
else:
    raise ValueError("score_method must be 'A', 'B', 'C', or 'D'")

ranked = (
    cluster_metrics
    .sort_values(
        by=[score_col, "n_patients", "n_sequences"],
        ascending=[False, False, False]
    )
)

#print ranked and only the relevant columns
print(ranked[[cluster_col,
            #    score_col,
                 "n_patients", 
            #    "sharedness", "dissemination",
              "tissue_evenness", "patient_evenness", 
                "sharedness_gini", "dissemination_gini",
              "cluster_size",  "size_weight" ]
              ].head(20))

print(ranked[[cluster_col,
            #    score_col,
                    "n_sequences",
              "cluster_size",  "size_weight" ]
              ].head(20))

top_clusters = ranked.head(top_n)[cluster_col].tolist()


 
pgens=compute_pgen_with_sonnia(df_shared,
                                 seq_col='junction_aa',
                                 vgene_col='v_gene',
                                 jgene_col='j_gene',
                                 pgen_model='humanIGH',
                                 max_length=40)
#rename junction_aa back to v_cdr3_j
pgens = pgens.rename(columns={'amino_acid': 'junction_aa'})
#check pgens columns for zero or negative values
print("Pgen values summary:")
print(pgens['Pgen'].describe())
#convert pgens to log10 scale, if 0 return na
pgens['Pgen'] = pgens['Pgen'].apply(lambda x: np.log10(x) if x > 0 else np.nan)

df_shared1 = df_shared.merge(pgens, on=['junction_aa','v_gene','j_gene'], how='left')
#read pgens
# pgens=pd.read_csv('/doctorai/chiarba/analysis/PER_CHIARA_df_sharing_malati_withpgen.tsv', sep='\t')
# antigens = pd.read_csv('/doctorai/chiarba/AbAg_database/fisher_results_annot_matches_enriched.tsv', sep='\t')
# #create a new column in antigens dataframe called junction_aa which is equal to CDR3_AA_sequence, + with C and W at the beginning and end
# antigens['junction_aa']  = 'C' + antigens['cdr3_ms'].astype(str).str.strip() + 'W'
# #deduplicate the dataframe
# antigens = antigens.drop_duplicates( )

# #show the duplicated rows based on junction_aa and antigen columns
# dups=antigens[antigens.duplicated(subset=['junction_aa', 'antigen'], keep=False)]
# print("Duplicated rows in antigens dataframe based on junction_aa and antigen columns:")
# print(dups)
# dups.to_csv('/doctorai/niccoloc/AGDB_duplicated_rows_antigens.tsv', sep='\t', index=False)



# ag_clusters= antigens['cluster_h1'].unique().tolist()

# df_shared1 = df_shared.merge(antigens[['junction_aa', 'antigen']].drop_duplicates(),  on=['junction_aa'], how='left')

# Final working dataframe
dfw = df_shared1[df_shared1[cluster_col].isin(top_clusters)].copy()
 

# Ensure 'study' exists
if study_col not in dfw.columns:
    dfw[study_col] = 'Unknown'

#sort dfw based on  the top clusters order
dfw[cluster_col] = pd.Categorical(dfw[cluster_col], categories=top_clusters, ordered=True)
dfw = dfw.sort_values(by=[cluster_col])

dfw['Pgen'].mean(skipna=True)
 
# #show cluster 626680 in dfw
#2506097 is mog
dfw[dfw[cluster_col]=='542137']['Pgen'] .unique()
dfw[dfw[cluster_col]=='542137']['junction_aa'].unique().tolist()
dfw[dfw[cluster_col]=='626680']['v_gene'].unique().tolist()
dfw[dfw[cluster_col]=='626680']['j_gene'].unique().tolist()

compute_pgen_with_sonnia(dfw[dfw[cluster_col]=='626680'],
                                 seq_col='junction_aa',
                                 vgene_col='v_gene',
                                 jgene_col='j_gene',
                                 pgen_model='humanIGH',
                                 max_length=40)


plot_network_legend_cols(dfw, cluster_col, seq_col, patient_col, study_col, tissue_col, alphabet, base_folder,
                          score_method,color_by='patients', column_width=0.07)
plot_network_legend_cols(dfw, cluster_col, seq_col, patient_col, study_col, tissue_col, alphabet, base_folder, score_method,color_by='tissues') 

#read antigen matches and check if any of the top clusters are present in the antigen matches
antigens = pd.read_csv('/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/agab_vs_ms_lev_matches1_withAgabMetadata_IGHV.tsv', sep='\t')
antigens['junction_aa']  = 'C' + antigens['cdr3_ms'].astype(str).str.strip() + 'W'
dfw_antigen = dfw.merge(antigens[['junction_aa', 'single_annotation','lev']].drop_duplicates(),  on=['junction_aa'], how='left')
dfw_antigen_all = df_shared1.merge(antigens[['junction_aa', 'single_annotation','lev']].drop_duplicates(),  on=['junction_aa'], how='left')

dfw_antigen[dfw_antigen['single_annotation'].notna()]['single_annotation'].value_counts()
dfw_antigen_all[dfw_antigen_all['single_annotation'].notna()]['single_annotation'].value_counts()
dfw_antigen_all[dfw_antigen_all['single_annotation'].notna()]['single_annotation'].tolist()
dfw_antigen_all[dfw_antigen_all['single_annotation']=='Myelin-oligodendrocyte glycoprotein']

dfw_antigen_all[dfw_antigen_all[cluster_col]=='2506097']['junction_aa'] 

dfw.columns





#------------ LD check


import pandas as pd
from rapidfuzz.distance import Levenshtein
from tqdm import tqdm

def prepare_reference_cdr3(df_ref: pd.DataFrame,
                           cdr3_col: str = "IGH CDR3 sequence",
                           vgene_col: str = "IGH variable gene") -> pd.DataFrame:
    """Return a cleaned reference dataframe with C/W-flanked CDR3 and length."""
    df_ref = df_ref.copy()
    df_ref["cdr3_full"] = "C" + df_ref[cdr3_col].str.strip() + "W"
    df_ref["cdr3_len"] = df_ref["cdr3_full"].str.len()
    if vgene_col not in df_ref.columns:
        df_ref[vgene_col] = None
    return df_ref

def find_cdr3_matches_fast(
    df_ref: pd.DataFrame,
    df_query: pd.DataFrame,
    seq_col_query: str = "sequence_aa",
    max_distance: int = 2,
    constrain_vgene: bool = False,
    vgene_ref_col: str = "IGH variable gene",
    vgene_query_col: str = "VGENE",
) -> pd.DataFrame:
    """
    Faster matching using RapidFuzz Levenshtein and some cheap prefilters:
    - group by V gene if constrain_vgene=True
    - filter by sequence length

    Returns one row per match.
    """
    df_ref = df_ref.copy()
    if "cdr3_full" not in df_ref.columns:
        df_ref["cdr3_full"] = "C" + df_ref["IGH CDR3 sequence"].str.strip() + "W"
    if "cdr3_len" not in df_ref.columns:
        df_ref["cdr3_len"] = df_ref["cdr3_full"].str.len()

    df_query = df_query.copy()
    df_query = df_query[df_query[seq_col_query].notna()]
    df_query["query_len"] = df_query[seq_col_query].str.len()

    results = []

    # If we constrain by V gene, group both tables by V gene
    if constrain_vgene:
        vgene_groups = list(df_ref.groupby(vgene_ref_col))
        for vgene, ref_block in tqdm(vgene_groups, desc="V-gene groups"):
            if pd.isna(vgene):
                continue
            query_block = df_query[df_query[vgene_query_col] == vgene]
            if query_block.empty:
                continue

            for ref_idx, ref_row in tqdm(ref_block.iterrows(), total=len(ref_block), leave=False):
                ref_seq = ref_row["cdr3_full"]
                ref_len = ref_row["cdr3_len"]

                # Length prefilter at dataframe level
                candidates = query_block[
                    (query_block["query_len"] >= ref_len - max_distance) &
                    (query_block["query_len"] <= ref_len + max_distance)
                ]

                for q_idx, q_row in candidates.iterrows():
                    q_seq = q_row[seq_col_query]
                    dist = Levenshtein.distance(ref_seq, q_seq)
                    if dist <= max_distance:
                        results.append({
                            "ref_index": ref_idx,
                            "query_index": q_idx,
                            "ref_cdr3_full": ref_seq,
                            "query_seq": q_seq,
                            "distance": dist,
                            "ref_vgene": vgene,
                            "query_vgene": q_row[vgene_query_col],
                        })

    else:
        # No V-gene constraint: still do length prefilter
        for ref_idx, ref_row in tqdm(df_ref.iterrows(), total=len(df_ref), desc="Reference sequences"):
            ref_seq = ref_row["cdr3_full"]
            ref_len = ref_row["cdr3_len"]

            candidates = df_query[
                (df_query["query_len"] >= ref_len - max_distance) &
                (df_query["query_len"] <= ref_len + max_distance)
            ]

            for q_idx, q_row in candidates.iterrows():
                q_seq = q_row[seq_col_query]
                dist = Levenshtein.distance(ref_seq, q_seq)
                if dist <= max_distance:
                    results.append({
                        "ref_index": ref_idx,
                        "query_index": q_idx,
                        "ref_cdr3_full": ref_seq,
                        "query_seq": q_seq,
                        "distance": dist,
                        "ref_vgene": ref_row.get(vgene_ref_col, None),
                        "query_vgene": q_row.get(vgene_query_col, None),
                    })

    return pd.DataFrame(results)

# Example usage:
# df_ref = prepare_reference_cdr3(df_ref)
# matches_fast = find_cdr3_matches_fast(
#     df_ref,
#     df_query,
#     seq_col_query="sequence_aa",
#     max_distance=2,
#     constrain_vgene=True,
#     vgene_ref_col="IGH variable gene",
#     vgene_query_col="VGENE",
# )



def annotate_query_with_min_distance(
    df_ref: pd.DataFrame,
    df_query: pd.DataFrame,
    seq_col_query: str = "sequence_aa",
    max_distance: int = 2,
    constrain_vgene: bool = True,
    vgene_ref_col: str = "IGH variable gene",
    vgene_query_col: str = "VGENE",
) -> pd.DataFrame:
    matches = find_cdr3_matches_fast(
        df_ref,
        df_query,
        seq_col_query=seq_col_query,
        max_distance=max_distance,
        constrain_vgene=constrain_vgene,
        vgene_ref_col=vgene_ref_col,
        vgene_query_col=vgene_query_col,
    )

    df_query = df_query.copy()
    df_query["min_distance_to_ref"] = max_distance + 1  # default: "no match"

    if not matches.empty:
        min_dist_per_query = matches.groupby("query_index")["distance"].min()
        for q_idx, d in min_dist_per_query.items():
            df_query.loc[q_idx, "min_distance_to_ref"] = d

    # Convenience flags for 0/1/2
    df_query["match_dist_0"] = df_query["min_distance_to_ref"] == 0
    df_query["match_dist_1"] = df_query["min_distance_to_ref"] <= 1
    df_query["match_dist_2"] = df_query["min_distance_to_ref"] <= 2

    return df_query


dfw.columns
cols_to_check = ['v_gene', 'd_gene',
       'j_gene', 'v_cdr3_j', 'study', 'tissue_type', 'junction_aa_length',
        'junction_length',   'junction',  
       'cdr3', 'cdr3_aa']
df_shared[cols_to_check]

#how many nans in cdr3_aa
df_shared['junction_aa'].isna().sum()

path_ms="/doctorai/niccoloc/MS_db/MS_BCR/ms_cdr3_lossius.csv"

path_ms= "/doctorai/chiarba/AbAg_database/aggregating_test/agab_cdr3_annotated1_with_origin_adjustantigen.csv"

ms_seqs = pd.read_csv( path_ms)
ms_seqs_prepared = prepare_reference_cdr3(ms_seqs,
                                         cdr3_col="cdr3_filtered",
                                         vgene_col="IGHV")
ms_seqs_prepared

matches_fast = find_cdr3_matches_fast(
    ms_seqs_prepared,
    df,
    seq_col_query="junction_aa",
    max_distance=3,
    constrain_vgene=False,
    vgene_ref_col="IGHV",
    vgene_query_col="v_gene",
)
matches_fast

#show matches_fast with distance 0,1,2
matches_fast[matches_fast['distance'].isin([2])]
