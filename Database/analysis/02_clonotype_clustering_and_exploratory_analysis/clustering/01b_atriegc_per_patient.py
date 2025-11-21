import pandas as pd
import atriegc

# Robust normalization
def normalize(s):
    return str(s).strip().upper().replace("_", "").replace("–", "-")

# Input / output paths
input_path  = "/doctorai/chiarba/analysis/MS_HC_naive.tsv"
output_path = "/doctorai/chiarba/analysis/MS_HC_naive_per_patient_h0h1.tsv"

# Load and normalize
df = pd.read_table(input_path, sep="\t", engine="python")

# Apply normalization to each column directly
df["v_gene"]      = df["v_gene"].apply(normalize)   # Normalize V genes
df["j_gene"]      = df["j_gene"].apply(normalize)   # Normalize J genes
df["junction_aa"] = df["junction_aa"].apply(normalize)  # Normalize amino-acid junctions

# clustering function unchanged…
def cluster_df(df, hamming_dist, label):
    """
    Build clusters (single linkage via TrieAA) within each subject_number + v/j gene pair.
    """
    cluster_map    = {}
    cluster_id_map = {}
    global_counter = 0

    # Group by subject, V and J
    for (subject, v_gene, j_gene), subdf in df.groupby(
        ["subject_number", "v_gene", "j_gene"]
    ):
        tr = atriegc.TrieAA()
        # Insert all junction sequences for this subject+genes
        for aa in subdf["junction_aa"].dropna():
            tr.insert(aa)
        # Get local clusters at this Hamming cutoff
        local_clusters = tr.clusters(hamming_dist)

        # Assign a unique global ID to each local cluster
        for local_cid in sorted(set(local_clusters.values())):
            global_counter += 1
            cluster_id_map[(subject, v_gene, j_gene, local_cid)] = global_counter

        # Map each sequence to its global cluster ID
        for seq, local_cid in local_clusters.items():
            cluster_map[(subject, v_gene, j_gene, seq)] = cluster_id_map[
                (subject, v_gene, j_gene, local_cid)
            ]

    # Safely fetch the cluster ID (or NA)
    def get_cluster_safe(row):
        if pd.isna(row["junction_aa"]):
            return pd.NA
        key = (row["subject_number"], row["v_gene"], row["j_gene"], row["junction_aa"])
        return cluster_map.get(key, pd.NA)

    return df.apply(get_cluster_safe, axis=1).rename(label)

# Compute clusters at H=0 and H=1, per subject
df["cluster_h0"] = cluster_df(df, hamming_dist=0, label="cluster_h0")
df["cluster_h1"] = cluster_df(df, hamming_dist=1, label="cluster_h1")

# Write out and report
df.to_csv(output_path, sep="\t", index=False)

print("Clustering completed.")
print("Unique clusters (H=0):", df["cluster_h0"].nunique(dropna=True))
print("Unique clusters (H=1):", df["cluster_h1"].nunique(dropna=True))
print("Rows without cluster H0:", df["cluster_h0"].isna().sum())
print("Rows without cluster H1:", df["cluster_h1"].isna().sum())
