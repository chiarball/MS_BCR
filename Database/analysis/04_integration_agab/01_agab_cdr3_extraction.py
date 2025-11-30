import ast
import pandas as pd

# English comment: input and output paths
input_path = "/doctorai/chiarba/AbAg_database/agab_merged.csv"
output_path = "/doctorai/chiarba/AbAg_database/clean/annotation/agab_with_cdr3.tsv"

# English comment: auto-detect separator to avoid parsing errors
with open(input_path, "r") as f:
    first_line = f.readline()

for sep in ["\t", ",", ";", "|"]:
    if sep in first_line:
        detected_sep = sep
        break
else:
    detected_sep = ","  # fallback

print(f"Detected separator: {repr(detected_sep)}")

# English comment: function to extract heavy-chain CDR3 from metadata string
def extract_heavy_cdr3(meta):
    if pd.isna(meta):
        return None
    try:
        data = ast.literal_eval(meta)         # Python dict parser
        heavy = data.get("heavy_riot_numbering", {})
        return heavy.get("cdr3_aa", None)
    except Exception:
        return None

def main():
    chunksize = 100_000
    first_chunk = True

    # English comment: stream the CSV in chunks
    for chunk in pd.read_csv(input_path, 
                              sep=detected_sep,
                              chunksize=chunksize):

        # English comment: metadata column check
        if "metadata" not in chunk.columns:
            print("ERROR: column 'metadata' not found. Columns are:")
            print(list(chunk.columns))
            return

        # English comment: extract heavy CDR3
        chunk["cdr3_aa"] = chunk["metadata"].apply(extract_heavy_cdr3)

        # English comment: append to output
        chunk.to_csv(
            output_path,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk
        )
        first_chunk = False

if __name__ == "__main__":
    main()
