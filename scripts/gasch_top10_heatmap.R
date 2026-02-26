#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

input_path <- "data/gasch2000.txt"
output_path <- "results/gene_top10_heatmap.png"
url <- "https://www.shackett.org/files/gasch2000.txt"

if (!file.exists(input_path)) {
  dir.create(dirname(input_path), recursive = TRUE, showWarnings = FALSE)
  download.file(url, destfile = input_path, mode = "wb")
}

raw <- read.delim(
  input_path,
  sep = "\t",
  header = TRUE,
  quote = "\"",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

if (ncol(raw) < 4) {
  stop("Expected at least 4 columns (UID + metadata + conditions).")
}

gene_col <- names(raw)[1]
non_expr_names <- c("UID", "NAME", "GWEIGHT", "EWEIGHT")
candidate_cols <- setdiff(names(raw), c(gene_col, intersect(names(raw), non_expr_names[-1])))

is_num_col <- function(x) {
  suppressWarnings(xn <- as.numeric(x))
  mean(!is.na(xn)) > 0.5
}

expr_cols <- candidate_cols[vapply(raw[candidate_cols], is_num_col, logical(1))]
if (length(expr_cols) == 0) {
  stop("No numeric expression columns detected.")
}
expr_cols <- expr_cols[seq_len(min(30, length(expr_cols)))]

expr <- as.data.frame(lapply(raw[expr_cols], function(x) suppressWarnings(as.numeric(x))))
rownames(expr) <- raw[[gene_col]]

na_count <- sum(is.na(as.matrix(expr)))
has_negatives <- any(as.matrix(expr) < 0, na.rm = TRUE)

gene_ids <- raw[[gene_col]]
unique_gene_ids <- length(unique(gene_ids))
cat(sprintf("Rows: %d\n", nrow(raw)))
cat(sprintf("Gene ID column: %s\n", gene_col))
cat(sprintf("Unique gene IDs: %d\n", unique_gene_ids))
cat(sprintf("Expression columns used: %d\n", length(expr_cols)))
cat(sprintf("Missing values (NA): %d\n", na_count))
cat(sprintf("Any negative values: %s\n", ifelse(has_negatives, "yes", "no")))
if (has_negatives) {
  cat("Interpretation: data already appears to be on a log/ratiometric scale; no log transform applied.\n")
} else {
  cat("Interpretation: no negatives detected; consider log2(x + 1) if raw-scale counts/intensities.\n")
}

mat <- as.matrix(expr)
row_mean <- rowMeans(mat, na.rm = TRUE)
row_sd <- apply(mat, 1, sd, na.rm = TRUE)

near_zero <- abs(row_mean) < 1e-8 | is.na(row_mean)
cv <- row_sd / abs(row_mean)
cv[near_zero] <- NA_real_

ord <- order(cv, decreasing = TRUE, na.last = NA)
if (length(ord) < 10) {
  stop("Fewer than 10 genes with finite CV values.")
}

top_idx <- ord[1:10]
top_genes <- rownames(mat)[top_idx]

top_mat <- mat[top_idx, , drop = FALSE]
long <- as.data.frame(as.table(top_mat), stringsAsFactors = FALSE)
names(long) <- c("gene", "conditions", "expression")
long$gene <- factor(long$gene, levels = top_genes)

p <- ggplot(long, aes(x = gene, y = conditions, fill = expression)) +
  geom_tile(color = "black", linewidth = 0.15) +
  scale_fill_distiller(palette = "PuBuGn", direction = 1, na.value = "grey90") +
  labs(title = "gene top 10", x = "gene", y = "conditions", fill = "expression") +
  theme_minimal(base_size = 11, base_family = "Times New Roman") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = p, width = 5, height = 5, units = "in", dpi = 300, bg = "white")

cat(sprintf("Saved heatmap: %s\n", output_path))
cat("Top 10 genes by CV (sd / abs(mean)):\n")
print(data.frame(gene = top_genes, cv = cv[top_idx], row.names = NULL))
git config --global user.name "ovaryachiever"
git config --global user.email lparker@unr.edu

