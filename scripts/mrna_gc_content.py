#!/usr/bin/env python3
"""Calculate GC content for yeast mRNA FASTA and plot distribution.

This script is intended to be driven by the one‑shot prompt described in the
homework assignment.  It will:

* download the UCSC sacCer3 mRNA FASTA if it is not already present,
* parse the gzipped FASTA, computing length and GC fraction for every record,
* write a tab‑separated table sorted by gc_content descending,
* draw a histogram with density curve, mean/median lines and a caption,
  and save it at the requested size and dpi.

Usage::

    python3 scripts/mrna_gc_content.py [--in data/mrna.fa.gz]
                                      [--out-table results/mrna_metrics.tsv]
                                      [--out-plot results/gc_content_distribution.png]

The default input file is downloaded automatically.
"""

import argparse
import gzip
import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def download_if_missing(path: Path, url: str) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url} -> {path}")
    subprocess.run(["curl", "-L", "-o", str(path), url], check=True)


def parse_mrna_fasta(fasta_path: Path):
    """Yield (accession, length, gc_content) for each record."""
    with gzip.open(fasta_path, "rt", encoding="utf-8") as handle:
        accession = None
        seq_lines = []
        for line in handle:
            if line.startswith(">"):
                if accession is not None:
                    seq = "".join(seq_lines)
                    length = len(seq)
                    if length > 0:
                        gc = (seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")) / length
                    else:
                        gc = 0.0
                    yield accession, length, gc
                header = line[1:].strip().split()[0]
                accession = header
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        # last record
        if accession is not None:
            seq = "".join(seq_lines)
            length = len(seq)
            gc = (seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")) / length if length > 0 else 0.0
            yield accession, length, gc


def make_plot(gc_array: np.ndarray, outpath: Path) -> None:
    n = len(gc_array)
    mean = float(np.mean(gc_array))
    median = float(np.median(gc_array))
    sd = float(np.std(gc_array, ddof=1)) if n > 1 else 0.0

    fig, ax = plt.subplots(figsize=(16, 9), dpi=200)
    # histogram
    counts, bins, patches = ax.hist(gc_array, bins=50, density=True, color="#bbbbff", alpha=0.7)
    # density curve
    from scipy.stats import gaussian_kde

    kde = gaussian_kde(gc_array)
    xs = np.linspace(0, 1, 200)
    ax.plot(xs, kde(xs), color="navy", lw=1.5)

    # mean/median lines
    ax.axvline(mean, color="red", linestyle="--", label=f"mean {mean:.4f}")
    ax.axvline(median, color="green", linestyle="--", label=f"median {median:.4f}")

    ax.set_xlabel("GC content")
    ax.set_ylabel("Density")
    ax.set_title("GC content distribution of yeast mRNA sequences")

    caption = f"n={n}, mean={mean:.4f}, median={median:.4f}, sd={sd:.4f}"
    fig.text(0.5, -0.05, caption, ha="center")

    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Compute GC content for Saccharomyces cerevisiae mRNA and plot"    )
    parser.add_argument("--in", dest="in_fasta", type=Path, default=Path("data/mrna.fa.gz"),
                        help="input gzipped mRNA fasta")
    parser.add_argument("--out-table", type=Path, default=Path("results/mrna_metrics.tsv"),
                        help="output TSV table")
    parser.add_argument("--out-plot", type=Path, default=Path("results/gc_content_distribution.png"),
                        help="output plot PNG")
    args = parser.parse_args()

    # download if necessary
    download_if_missing(args.in_fasta,
                        "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz")

    # parse records
    rows = []
    for acc, length, gc in parse_mrna_fasta(args.in_fasta):
        rows.append((acc, length, gc))

    # build dataframe
    df = pd.DataFrame(rows, columns=["accession", "length", "gc_content"])
    df = df.sort_values("gc_content", ascending=False)

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_table, sep="\t", index=False, float_format="%.4f")

    gc_array = df["gc_content"].to_numpy(dtype=float)
    args.out_plot.parent.mkdir(parents=True, exist_ok=True)
    make_plot(gc_array, args.out_plot)


if __name__ == "__main__":
    main()
