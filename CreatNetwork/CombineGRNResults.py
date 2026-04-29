import os
import glob
import pandas as pd

BASE_DIR = "/vast/palmer/scratch/noonan/yj345/snRNA/grn"

TISSUES = ["FL", "HL", "PA1", "PA2"]

FREQ_CUTOFF = 0.95

def combine_one_tissue(tissue):
    base = os.path.join(BASE_DIR, tissue, "adj_sub")
    pattern = os.path.join(base, "adj_sub_*.tsv")

    files = sorted(glob.glob(pattern))

    print(f"Combining tissue: {tissue}")
    
    dfs = []

    for f in files:
        df = pd.read_csv(f, sep="\t")

        expected_cols = {"TF", "target", "importance"}

        if not expected_cols.issubset(df.columns):
            raise ValueError(
                f"{f} does not have required columns {expected_cols}. "
                f"Found columns: {list(df.columns)}"
            )

        df["source_file"] = os.path.basename(f)
        dfs.append(df)

    all_edges = pd.concat(dfs, ignore_index=True)
    print(f"{tissue}: total edge rows = {all_edges.shape[0]}")

    summary = (
        all_edges
        .groupby(["TF", "target"], as_index=False)
        .agg(
            freq=("importance", "count"),
            mean_importance=("importance", "mean"),
            max_importance=("importance", "max")
        )
    )

    n_runs = len(files)
    threshold = int(FREQ_CUTOFF * n_runs)

    print(
        f"{tissue}: {n_runs} runs; keeping edges with "
        f"freq >= {threshold} (>= {FREQ_CUTOFF * 100:.0f}% of runs)"
    )

    consensus = summary[summary["freq"] >= threshold].copy()
    consensus.sort_values("mean_importance", ascending=False, inplace=True)

    out_file = os.path.join(
        base,
        f"adj_{tissue}_consensus_{n_runs}runs_{int(FREQ_CUTOFF * 100)}pct.tsv"
    )

    consensus.to_csv(out_file, sep="\t", index=False)

    print(f"{tissue}: wrote consensus network with {consensus.shape[0]} edges to:")
    print(f"  {out_file}")
    print()


def main():
    for tissue in TISSUES:
        combine_one_tissue(tissue)

    print("All tissues finished.")


if __name__ == "__main__":
    main()
