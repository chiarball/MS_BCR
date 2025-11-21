import pandas as pd
import atriegc

# Funzione di normalizzazione robusta
def normalize(s):
    return str(s).strip().upper().replace("_", "").replace("–", "-")

# Nuovo percorso input
input_path = "/doctorai/chiarba/analysis/MS_HC_naive.tsv"

# Nuovo percorso output
output_path = "/doctorai/chiarba/analysis/MS_HC_naive_h0h1.tsv"

# Carica e normalizza
df = pd.read_table(input_path, sep="\t", engine="python")
df["v_gene"] = df["v_gene"].apply(normalize)
df["j_gene"] = df["j_gene"].apply(normalize)
df["junction_aa"] = df["junction_aa"].apply(normalize)

# Funzione generica di clustering per un dato hamming
def cluster_df(df, hamming_dist, label):
    cluster_map = {}
    cluster_id_map = {}
    global_counter = 0

    for (v_gene, j_gene), subdf in df.groupby(["v_gene", "j_gene"]):
        tr = atriegc.TrieAA()
        for aa in subdf["junction_aa"].dropna():
            tr.insert(aa)
        local_clusters = tr.clusters(hamming_dist)
        for local_cid in sorted(set(local_clusters.values())):
            global_counter += 1
            cluster_id_map[(v_gene, j_gene, local_cid)] = global_counter
        for seq, local_cid in local_clusters.items():
            cluster_map[(v_gene, j_gene, seq)] = cluster_id_map[(v_gene, j_gene, local_cid)]

    def get_cluster_safe(r):
        if pd.isna(r["junction_aa"]):
            return pd.NA
        key = (r["v_gene"], r["j_gene"], r["junction_aa"])
        return cluster_map.get(key, pd.NA)

    return df.apply(get_cluster_safe, axis=1).rename(label)

# Applica clustering con hamming = 0
df["cluster_h0"] = cluster_df(df, hamming_dist=0, label="cluster_h0")

# Applica clustering con hamming = 1
df["cluster_h1"] = cluster_df(df, hamming_dist=1, label="cluster_h1")

# Esporta
df.to_csv(output_path, sep="\t", index=False)

# Verifica finale
print("Clustering completato.")
print("Cluster unici (H=0):", df["cluster_h0"].nunique(dropna=True))
print("Cluster unici (H=1):", df["cluster_h1"].nunique(dropna=True))
print("Righe senza cluster H0:", df["cluster_h0"].isna().sum())
print("Righe senza cluster H1:", df["cluster_h1"].isna().sum())
