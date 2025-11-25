import pandas as pd
from infer_lineages import InferLineages
from pathlib import Path

# Your full list of subject IDs
patients = [
  "Ru8"
]

def main():
    # Path to your single TSV
    input_tsv = Path("/doctorai/chiarba/analysis/hilary/MS_HC_sub_last.tsv")
    output_dir = "/doctorai/chiarba/analysis/hilary/output/"

    # Load it once
    df_all = pd.read_csv(input_tsv, sep="\t", dtype=str)
    df_sub = {}
    df_out = {}

    # Loop over each patient
    for subject_id in patients:
        df_sub[subject_id] = df_all[df_all["subject_number"] == subject_id]
        if df_sub[subject_id].empty:
            print(f"Warning: no records found for {subject_id}")
            continue

        # instantiate and call the right method on InferLineages
        model = InferLineages()
        # if your class method is named `infer`, call:
        df_out[subject_id] = model.infer(df_sub[subject_id])
        # or if it’s called `run` or `fit`, use that name instead
        output_file = output_dir + "{}_hilary.tsv".format(subject_id)
        df_out[subject_id].to_csv(output_file, sep="\t", index=False)

    # Uncomment if you want a merged file
    # df_tot = pd.concat([df_out[subject_id] for subject_id in patients], ignore_index=True)
    # df_tot.to_csv(output_dir + "Output_hilary_tot.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
