# R golden references for the fgsea and goseq ports (M8-12; GIPPA_SPEC.md §11.5).
# Runs only inside the conda env gippa-golden-r (scripts/golden/golden_r.yml), through
# scripts/slurm/golden_r.sbatch or .github/workflows/golden-r.yml; never from src/.
#   usage: make_golden_r.R <inputs_dir> <outputs_dir> <versions_txt> <lock_file> <GOLDEN_SEED>
# Writes outputs/{fgsea,calc_gsea_stat,goseq}.tsv.gz (numbers with format(x, digits = 17)) and
# VERSIONS.txt. Installs nothing: every package comes from conda-forge or bioconda.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("usage: make_golden_r.R <inputs> <outputs> <versions> <lock> <seed>")
inputs <- args[1]
outputs <- args[2]
versions_file <- args[3]
lock_file <- args[4]
GOLDEN_SEED <- as.integer(args[5])
stopifnot(!is.na(GOLDEN_SEED))

suppressPackageStartupMessages({
  library(fgsea)
  library(goseq)
  library(BiasedUrn)
})

warned <- character()
capture <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

read_tsv <- function(name, classes) {
  read.delim(gzfile(file.path(inputs, name)), colClasses = classes, quote = "",
             comment.char = "")
}

num <- function(x) {
  vapply(x, function(v) if (is.na(v)) "NA" else format(v, digits = 17), character(1),
         USE.NAMES = FALSE)
}

write_tsv <- function(df, name) {
  path <- file.path(outputs, name)
  con <- gzfile(path, "wb", compression = 9)
  lines <- c(paste(names(df), collapse = "\t"), do.call(paste, c(df, sep = "\t")))
  writeLines(lines, con, sep = "\n", useBytes = TRUE)
  close(con)
  if (file.size(path) > 5 * 1024^2) stop(name, " exceeds 5 MB (R11)")
}

dir.create(outputs, recursive = TRUE, showWarnings = FALSE)

# --- fgsea ------------------------------------------------------------------------------------
st <- read_tsv("gsea_stats.tsv.gz", c("character", "numeric"))
stats <- setNames(st$score, st$tu_id)
gs <- read_tsv("gsea_sets.tsv.gz", c("character", "character"))
pathways <- split(gs$member, factor(gs$set_id, levels = unique(gs$set_id)))

set.seed(GOLDEN_SEED)
fg <- capture(fgseaMultilevel(pathways = pathways, stats = stats, sampleSize = 101,
                              minSize = 10, maxSize = 500, eps = 0, nPermSimple = 1000,
                              scoreType = "std", nproc = 1))
fg <- fg[match(intersect(names(pathways), fg$pathway), fg$pathway), ]
write_tsv(data.frame(
  pathway = fg$pathway, pval = num(fg$pval), log2err = num(fg$log2err), ES = num(fg$ES),
  NES = num(fg$NES), size = as.character(fg$size),
  leadingEdge = vapply(fg$leadingEdge, paste, character(1), collapse = ";"),
  stringsAsFactors = FALSE
), "fgsea.tsv.gz")

# calcGseaStat for every set with >= 1 member in the scores, on the same ranking fgsea uses.
ranked <- sort(stats, decreasing = TRUE)
matched <- lapply(pathways, function(p) unique(na.omit(match(p, names(ranked)))))
matched <- matched[lengths(matched) > 0]
es <- vapply(matched, function(idx) calcGseaStat(ranked, selectedStats = idx, gseaParam = 1,
                                                  scoreType = "std"), numeric(1))
write_tsv(data.frame(pathway = names(matched), ES = num(es), stringsAsFactors = FALSE),
          "calc_gsea_stat.tsv.gz")

# --- goseq (Wallenius, supplied PWF = p_detect) -----------------------------------------------
genes <- read_tsv("ora_genes.tsv.gz", c("character", "numeric", "integer"))
cats <- read_tsv("ora_categories.tsv.gz", c("character", "character"))
pwf <- data.frame(DEgenes = genes$detected, bias.data = 1, pwf = genes$p_detect,
                  row.names = genes$gene)
gene2cat <- data.frame(gene = cats$gene, category = cats$category, stringsAsFactors = FALSE)

set.seed(GOLDEN_SEED)
go <- capture(goseq(pwf, gene2cat = gene2cat, method = "Wallenius",
                    use_genes_without_cat = TRUE))
go <- go[match(unique(cats$category), go$category), ]
write_tsv(data.frame(
  category = go$category, over_represented_pvalue = num(go$over_represented_pvalue),
  under_represented_pvalue = num(go$under_represented_pvalue),
  numDEInCat = as.character(go$numDEInCat), numInCat = as.character(go$numInCat),
  stringsAsFactors = FALSE
), "goseq.tsv.gz")

# --- VERSIONS.txt -----------------------------------------------------------------------------
lock_sha <- if (file.exists(lock_file)) {
  strsplit(system2("sha256sum", lock_file, stdout = TRUE), " ")[[1]][1]
} else {
  "absent"
}
out <- c(
  paste("fgsea", as.character(packageVersion("fgsea"))),
  paste("goseq", as.character(packageVersion("goseq"))),
  paste("BiasedUrn", as.character(packageVersion("BiasedUrn"))),
  paste("env lock sha256", lock_sha),
  paste("GOLDEN_SEED", GOLDEN_SEED),
  "",
  "## warnings",
  if (length(warned)) paste("-", unique(warned)) else "none",
  "",
  "## sessionInfo()",
  capture.output(print(sessionInfo()))
)
writeLines(out, versions_file)
