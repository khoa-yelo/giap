#!/usr/bin/env -S Rscript --vanilla
# fgseaMultilevel for GIPPA's function step (M8-17a; ADR-0038; GIPPA_SPEC.md §7.7.4, R5).
# Started only as a subprocess by gippa.io.fgsea_r in the gippa-fgsea env (workflow/envs/fgsea.yml).
# Installs and downloads nothing.
#   gippa_fgsea.R --stats stats.tsv --sets sets.tsv --jobs jobs.tsv --sample-size 101 --eps 0 \
#                --n-perm-simple 1000 --seed <int> --out results.tsv --versions versions.txt
# Exit codes: 0 ok, 2 usage error, 3 invalid input, 1 any R error (R's default).
# Per job: set.seed(seed), then fgseaMultilevel(..., scoreType = "std", nproc = 1); rows in the
# input set order. results.tsv and versions.txt are written to *.tmp and renamed, never partial.

usage_error <- function(msg) {
  cat("gippa_fgsea.R: usage error: ", msg, "\n", file = stderr(), sep = "")
  quit(status = 2, save = "no")
}

invalid_input <- function(msg) {
  cat("gippa_fgsea.R: invalid input: ", msg, "\n", file = stderr(), sep = "")
  quit(status = 3, save = "no")
}

# ---------------------------------------------------------------- arguments (exit 2)

OPTIONS <- c("--stats", "--sets", "--jobs", "--sample-size", "--eps", "--n-perm-simple",
             "--seed", "--out", "--versions")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) %% 2 != 0) usage_error("options must come as --name value pairs")
keys <- argv[c(TRUE, FALSE)]
vals <- argv[c(FALSE, TRUE)]
unknown <- setdiff(keys, OPTIONS)
if (length(unknown)) usage_error(paste("unknown argument", unknown[1]))
if (anyDuplicated(keys)) usage_error(paste("repeated argument", keys[anyDuplicated(keys)]))
missing_opts <- setdiff(OPTIONS, keys)
if (length(missing_opts)) usage_error(paste("missing argument", missing_opts[1]))
opt <- setNames(as.list(vals), keys)

as_int <- function(name, lo) {
  x <- opt[[name]]
  if (!grepl("^-?[0-9]+$", x)) usage_error(paste(name, "is not an integer:", x))
  v <- suppressWarnings(as.integer(x))
  if (is.na(v) || v < lo) usage_error(paste(name, "is out of range:", x))
  v
}
seed <- as_int("--seed", -.Machine$integer.max)
sample_size <- as_int("--sample-size", 1L)
n_perm_simple <- as_int("--n-perm-simple", 1L)
eps <- suppressWarnings(as.numeric(opt[["--eps"]]))
if (is.na(eps) || !is.finite(eps) || eps < 0) usage_error(paste("--eps is invalid:", opt[["--eps"]]))

# ---------------------------------------------------------------- inputs (exit 3)

read_tab <- function(path, cols) {
  if (!file.exists(path)) invalid_input(paste("no such file", path))
  df <- read.delim(path, colClasses = "character", quote = "", comment.char = "",
                   na.strings = character(0), check.names = FALSE, encoding = "UTF-8")
  if (!identical(names(df), cols)) {
    invalid_input(paste0(path, ": columns must be ", paste(cols, collapse = ", ")))
  }
  df
}

st <- read_tab(opt[["--stats"]], c("tu_id", "score"))
score <- suppressWarnings(as.numeric(st$score))
if (any(!is.finite(score))) {
  invalid_input(paste("non-finite or NA score for", st$tu_id[!is.finite(score)][1]))
}
if (anyDuplicated(st$tu_id)) invalid_input(paste("duplicate tu_id", st$tu_id[anyDuplicated(st$tu_id)]))
stats <- setNames(score, st$tu_id)

jb <- read_tab(opt[["--jobs"]], c("job_id", "min_size", "max_size"))
if (anyDuplicated(jb$job_id)) invalid_input(paste("duplicate job", jb$job_id[anyDuplicated(jb$job_id)]))
int_col <- function(x, name) {
  ok <- grepl("^[0-9]+$", x)
  if (!all(ok)) invalid_input(paste(name, "is not an integer:", x[!ok][1]))
  as.integer(x)
}
min_size <- int_col(jb$min_size, "min_size")
max_size <- int_col(jb$max_size, "max_size")
if (any(min_size < 1)) invalid_input("min_size < 1")
if (any(max_size < min_size)) invalid_input("max_size < min_size")

gs <- read_tab(opt[["--sets"]], c("job_id", "set_id", "member"))
stray <- setdiff(unique(gs$job_id), jb$job_id)
if (length(stray)) invalid_input(paste("sets.tsv job", stray[1], "is not in jobs.tsv"))
unranked <- setdiff(gs$member, names(stats))
if (length(unranked)) invalid_input(paste("set member", unranked[1], "is not in stats"))

# ---------------------------------------------------------------- fgsea (R errors exit 1)

suppressPackageStartupMessages(library(fgsea))

num <- function(x) {
  vapply(x, function(v) if (is.na(v)) "NA" else format(v, digits = 17), character(1),
         USE.NAMES = FALSE)
}

warned <- character()
warned_job <- character()
out <- "job_id\tpathway\tsize\tES\tNES\tpval\tlog2err\tleadingEdge"
for (i in seq_len(nrow(jb))) {
  job <- jb$job_id[i]
  rows <- gs[gs$job_id == job, , drop = FALSE]
  pathways <- split(rows$member, factor(rows$set_id, levels = unique(rows$set_id)))
  sizes <- vapply(pathways, function(p) length(unique(p)), integer(1))
  if (!any(sizes >= min_size[i] & sizes <= max_size[i])) next  # no rows, not an error
  set.seed(seed)
  fg <- withCallingHandlers(
    fgseaMultilevel(pathways, stats, sampleSize = sample_size, minSize = min_size[i],
                    maxSize = max_size[i], eps = eps, nPermSimple = n_perm_simple,
                    scoreType = "std", nproc = 1),
    warning = function(w) {
      msg <- gsub("[\r\n]+", " ", conditionMessage(w))
      if (!any(warned == msg & warned_job == job)) {
        warned <<- c(warned, msg)
        warned_job <<- c(warned_job, job)
      }
      invokeRestart("muffleWarning")
    }
  )
  fg <- fg[match(intersect(names(pathways), fg$pathway), fg$pathway), ]
  if (nrow(fg) == 0) next
  edge <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
  out <- c(out, paste(job, fg$pathway, as.character(fg$size), num(fg$ES), num(fg$NES),
                      num(fg$pval), num(fg$log2err), edge, sep = "\t"))
}

versions <- c(
  paste0("R\t", R.version.string),
  paste0("fgsea\t", as.character(packageVersion("fgsea"))),
  paste0("data.table\t", as.character(packageVersion("data.table"))),
  paste0("BiocParallel\t", as.character(packageVersion("BiocParallel"))),
  paste0("seed\t", seed),
  paste0("jobs\t", nrow(jb)),
  if (length(warned)) paste0("warning\t", warned_job, "\t", warned)
)

write_atomic <- function(lines, path) {
  tmp <- paste0(path, ".tmp")
  con <- file(tmp, "wb")
  writeLines(lines, con, sep = "\n", useBytes = TRUE)
  close(con)
  if (!file.rename(tmp, path)) stop("cannot rename ", tmp, " to ", path)
}
write_atomic(versions, opt[["--versions"]])
write_atomic(out, opt[["--out"]])
