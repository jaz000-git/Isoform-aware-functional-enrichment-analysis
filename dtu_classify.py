#!/usr/bin/env python3

import pandas as pd
import argparse

def classify(row): #define foreground / background
    if row["padj"] < 0.05:
        if row["logFC"] > 0.5:
            return "cond"
        elif row["logFC"] < -0.5:
            return "ref"
    return "non_sig"

def main():
    parser = argparse.ArgumentParser(description="Classify DTU transcripts")
    parser.add_argument("--dtu", required=True)
    parser.add_argument("--out", default="dtu_classified.tsv")
    args = parser.parse_args()

    dtu = pd.read_csv(args.dtu)

    dtu = dtu.rename(columns={"isoform_id": "transcript"})

    # normalize IDs
    dtu["transcript"] = dtu["transcript"].str.replace(r"\.\d+$", "", regex=True)

    dtu["condition"] = dtu.apply(classify, axis=1)

    dtu[["transcript", "condition"]].to_csv(args.out, sep="\t", index=False)

    print("\nCounts:")
    print(dtu["condition"].value_counts())

if __name__ == "__main__":
    main()