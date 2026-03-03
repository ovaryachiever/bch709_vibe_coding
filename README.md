# Yeast Feature Distribution Analysis

This project analyzes how genomic features are distributed across *Saccharomyces cerevisiae* chromosomes and tests whether feature density correlates with chromosome size.

## Inputs

- `data/saccharomyces_cerevisiae.gff.gz`
- `data/chrom.sizes`

## Features analyzed

- `gene`
- `exon`
- `tRNA`
- `snoRNA`

Note: in this GFF, `exon` entries are absent (`noncoding_exon` is used instead), so `exon` counts are zero.

## Environment setup

```bash
conda env create -f environment.python.yml
conda activate bch709-yeast-features
```

If environment already exists:

```bash
conda env update -f environment.python.yml --prune
conda activate bch709-yeast-features
```

## Run analysis

```bash
python scripts/yeast_feature_distribution.py
```

## Outputs

- `results/feature_density_by_chromosome.tsv`
  - Per-chromosome size, feature counts, and feature density (per Mb)
- `results/feature_size_correlations.tsv`
  - Pearson correlations between chromosome size and:
    - raw feature count
    - feature density
- `results/feature_distribution_density.png`
  - Composite figure with:
    - stacked counts per chromosome
    - density-vs-size scatterplots with trend lines

## Main findings from current run

- `gene` count vs chromosome size: very strong positive correlation (`r = 0.9988`)
- `tRNA` count vs chromosome size: moderate positive (`r = 0.6126`)
- `snoRNA` count vs chromosome size: moderate positive (`r = 0.6189`)
- `exon`: no variation (all zero in this annotation)
- `tRNA` density tends to decrease with chromosome size (`r = -0.4609`)

## mRNA GC-content analysis (added)

I computed GC content for all UCSC sacCer3 mRNA records and produced `results/mrna_metrics.tsv` and `results/gc_content_distribution.png`.

Interpretation: The GC-content distribution shows (1) a subset of transcripts with relatively high GC, which can reflect selection for translation efficiency and stable mRNA structures in highly expressed genes, and (2) a lower-GC subpopulation that may arise from regional mutational biases and GC-biased gene conversion. These patterns are consistent with a mix of selective and neutral processes shaping transcript GC content.

Files:

- `prompt.txt`: the one-shot prompt used to generate the analysis script.
- `scripts/mrna_gc_content.py`: script that produced the metrics and the plot.
- `results/mrna_metrics.tsv` and `results/gc_content_distribution.png` are included.

## Manual

For a full step-by-step workflow, see `MANUAL.md`.
