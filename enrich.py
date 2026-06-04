#!/usr/bin/env python3

import pandas as pd
import argparse
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2_contingency, fisher_exact
import numpy as np


def run_enrichment(merged, dtu, target):

    # foreground/background transcript sets
    fg = set(dtu[dtu["condition"] == target]["transcript"])
    bg = set(dtu[dtu["condition"] != target]["transcript"])

    n_fg = len(fg)
    n_bg = len(bg)

    results = []

    for go, sub in merged.groupby("GO"):

        go_set = set(sub["transcript"])

        # contingency table counts
        a = len(go_set & fg)   # FG hits
        b = len(fg - go_set)  # FG not in GO
        c = len(go_set & bg)  # BG hits
        d = len(bg - go_set)  # BG not in GO

        if a + b == 0 or c + d == 0:
            continue

        table = np.array([
            [a, b],
            [c, d]
        ])

        # chi-square
        chi2, p_chi, dof, expected = chi2_contingency(
            table,
            correction=False
        )

        # Fisher if expected counts less than 5
        if (expected < 5).any():
            _, p = fisher_exact(table, alternative="greater")
            test_used = "fisher"
        else:
            p = p_chi
            test_used = "chi2"

        # odds ratio with Haldane correction
        OR = (
            ((a + 0.5) * (d + 0.5)) /
            ((b + 0.5) * (c + 0.5))
        )

        # calculate transcript ratio
        transcript_ratio = a / n_fg if n_fg > 0 else np.nan

        transcript_ratio_str = f"{a}/{n_fg}"

        # calculate background ratio

        bg_ratio = c / n_bg if n_bg > 0 else np.nan

        bg_ratio_str = f"{c}/{n_bg}"

        # Core enrichment gene list

        core_genes = (
            sub[sub["transcript"].isin(fg)]["gene"]
            .dropna()
            .unique()
        )

        core_genes_str = ";".join(core_genes)

        results.append([
            go,
            a,
            c,
            transcript_ratio,
            transcript_ratio_str,
            bg_ratio,
            bg_ratio_str,
            OR,
            p,
            test_used,
            core_genes_str
        ])

    df = pd.DataFrame(results, columns=[
        "GO",
        "fg_hits",
        "bg_hits",
        "TranscriptRatio",
        "TranscriptRatio_raw",
        "BgRatio",
        "BgRatio_raw",
        "OR",
        "p",
        "test",
        "genes"
    ])

    # multiple testing correction

    if len(df):
        df["FDR"] = multipletests(
            df["p"],
            method="fdr_bh"
        )[1]

        df = df.sort_values("FDR")

    return df


def main():

    parser = argparse.ArgumentParser(
        description="Merge + GO enrichment"
    )

    parser.add_argument("--go", required=True)
    parser.add_argument("--dtu", required=True)
    parser.add_argument("--out_prefix", default="go")

    args = parser.parse_args()

    go = pd.read_csv(args.go, sep="\t")
    dtu = pd.read_csv(args.dtu, sep="\t")

    merged = go.merge(
        dtu,
        on="transcript",
        how="inner"
    )

    print("\nMerged rows:", len(merged))

    for cond in ["cond", "ref"]:

        print(f"\nRunning {cond} enrichment...")

        res = run_enrichment(
            merged,
            dtu,
            cond
        )

        out = f"{args.out_prefix}_{cond}.tsv"

        res.to_csv(
            out,
            sep="\t",
            index=False
        )

        print(f"Saved: {out}")


if __name__ == "__main__":
    main()