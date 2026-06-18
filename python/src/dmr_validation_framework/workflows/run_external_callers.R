#!/usr/bin/env Rscript

arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name)
  idx <- which(args == key)
  if (length(idx) == 0 || idx[length(idx)] >= length(args)) {
    return(default)
  }
  args[idx[length(idx)] + 1]
}

arg_values <- function(name) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name)
  idx <- which(args == key)
  values <- character()
  for (i in idx) {
    if (i < length(args)) {
      values <- c(values, args[i + 1])
    }
  }
  values
}

split_csv <- function(value) {
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

root <- normalizePath(arg_value("root"), mustWork = FALSE)
out_root <- normalizePath(arg_value("out-root", file.path(root, "external_callers")), mustWork = FALSE)
contexts <- toupper(split_csv(arg_value("contexts", "CG,CHG,CHH")))
callers <- split_csv(arg_value("callers", "DSS,methylKit,dmrseq,BSmooth,comb-p,metilene"))
min_coverage <- as.integer(arg_value("min-coverage", "5"))
min_sites <- as.integer(arg_value("min-sites", "3"))
max_gap <- as.integer(arg_value("max-gap", "1000"))
q_threshold <- as.numeric(arg_value("q-threshold", "0.10"))
window_size <- as.integer(arg_value("window-size", "1000"))
step_size <- as.integer(arg_value("step-size", "1000"))
# methylKit effect-size threshold in percentage points (0-100). NOTE: this is
# caller-specific and not directly comparable with DSS/metilene delta (fraction)
# or BSmooth areaStat; methylKit::getMethylDiff filters on |meth.diff| >= this.
methylkit_min_diff <- as.numeric(arg_value("methylkit-min-diff", "10"))
if (is.na(methylkit_min_diff) || methylkit_min_diff < 0) {
  methylkit_min_diff <- 10
}
max_sites_per_context <- as.integer(arg_value("max-sites-per-context", "0"))
read_chunk_size <- as.integer(arg_value("read-chunk-size", "1000000"))
if (is.na(read_chunk_size) || read_chunk_size <= 0) {
  read_chunk_size <- 1000000
}
dmrseq_permutations <- as.integer(arg_value("dmrseq-permutations", "0"))
if (is.na(dmrseq_permutations) || dmrseq_permutations < 0) {
  dmrseq_permutations <- 0
}
metilene_bin <- arg_value("metilene", "")

extra_libs <- arg_values("r-lib")
if (length(extra_libs) > 0) {
  .libPaths(c(extra_libs, Sys.getenv("R_LIBS_USER"), .libPaths()))
} else {
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

status_rows <- list()
add_status <- function(caller, context, status, message = "") {
  status_rows[[length(status_rows) + 1]] <<- data.frame(
    caller = caller,
    context = context,
    status = status,
    message = message,
    stringsAsFactors = FALSE
  )
}

write_status <- function() {
  status <- if (length(status_rows) > 0) {
    do.call(rbind, status_rows)
  } else {
    data.frame(caller = character(), context = character(), status = character(), message = character())
  }
  utils::write.table(
    status,
    file.path(out_root, "external_caller_run_status.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

require_pkg <- function(pkg, caller, context = "") {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    add_status(caller, context, "missing_package", paste0("Missing R package: ", pkg))
  }
  ok
}

if (!require_pkg("data.table", "preflight")) {
  write_status()
  stop("data.table is required")
}
library(data.table)
library(parallel)
detectCores <- parallel::detectCores
if (requireNamespace("matrixStats", quietly = TRUE)) {
  rowVars <- matrixStats::rowVars
}

manifest_path <- file.path(root, "metadata", "sample_manifest.tsv")
if (!file.exists(manifest_path)) {
  write_status()
  stop(paste("Missing sample manifest:", manifest_path))
}
manifest <- data.table::fread(manifest_path)
required_manifest <- c("sample_id", "condition", "file")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest) > 0) {
  write_status()
  stop(paste("sample_manifest.tsv misses columns:", paste(missing_manifest, collapse = ",")))
}
manifest[, raw_path := file.path(root, "raw_cx", file)]
manifest[file.exists(file), raw_path := file]
if (!all(file.exists(manifest$raw_path))) {
  write_status()
  stop(paste("Missing raw CX files:", paste(manifest$raw_path[!file.exists(manifest$raw_path)], collapse = "; ")))
}

sample_ids <- as.character(manifest$sample_id)
control_ids <- as.character(manifest[condition == "control", sample_id])
treatment_ids <- as.character(manifest[condition == "treatment", sample_id])
treatment_vector <- as.integer(manifest$condition == "treatment")

cx_cols <- c("chrom", "pos", "strand", "meth", "unmeth", "context", "trinuc")

read_cx_chunk <- function(con) {
  tryCatch(
    utils::read.table(
      con,
      sep = "\t",
      header = FALSE,
      col.names = cx_cols,
      colClasses = c("character", "integer", "character", "integer", "integer", "character", "character"),
      quote = "",
      comment.char = "",
      nrows = read_chunk_size
    ),
    error = function(e) {
      if (grepl("no lines available", conditionMessage(e), ignore.case = TRUE)) {
        return(NULL)
      }
      stop(e)
    }
  )
}

read_cx_counts <- function(path, context) {
  ctx <- context
  message("[read] ", basename(path), " context=", ctx)
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  chunks <- list()
  chunks_read <- 0L
  rows_read <- 0L
  rows_kept <- 0L

  while (TRUE) {
    x <- read_cx_chunk(con)
    if (is.null(x) || nrow(x) == 0) {
      break
    }
    chunks_read <- chunks_read + 1L
    rows_read <- rows_read + nrow(x)
    x <- data.table::as.data.table(x)
    x <- x[x[["context"]] == ctx]
    if (nrow(x) > 0) {
      bad_contexts <- setdiff(unique(x[["context"]]), ctx)
      if (length(bad_contexts) > 0) {
        stop(
          "Context filter failed for ",
          basename(path),
          ": requested ",
          ctx,
          ", observed ",
          paste(bad_contexts, collapse = ",")
        )
      }
      x[, N := meth + unmeth]
      x <- x[N >= min_coverage]
    }
    if (nrow(x) > 0) {
      x[, X := meth]
      x <- x[, .(chr = chrom, pos = as.integer(pos), N = as.integer(N), X = as.integer(X))]
      rows_kept <- rows_kept + nrow(x)
      chunks[[length(chunks) + 1]] <- x
      if (max_sites_per_context > 0) {
        keep <- data.table::rbindlist(chunks, use.names = TRUE)
        data.table::setorder(keep, -N)
        if (nrow(keep) > max_sites_per_context) {
          keep <- keep[seq_len(max_sites_per_context)]
        }
        chunks <- list(keep)
      }
    }
    if (chunks_read %% 10 == 0) {
      message("[read] ", basename(path), " chunks=", chunks_read, " rows=", rows_read, " kept=", rows_kept)
    }
  }

  if (length(chunks) == 0) {
    x <- data.table::data.table(chr = character(), pos = integer(), N = integer(), X = integer())
  } else {
    x <- data.table::rbindlist(chunks, use.names = TRUE)
  }
  if (max_sites_per_context > 0 && nrow(x) > max_sites_per_context) {
    data.table::setorder(x, -N)
    x <- x[seq_len(max_sites_per_context)]
  }
  data.table::setorder(x, chr, pos)
  message("[read] ", basename(path), " kept ", nrow(x), " ", ctx, " sites")
  as.data.frame(x)
}

read_context_counts <- function(context) {
  out <- vector("list", length(sample_ids))
  names(out) <- sample_ids
  for (i in seq_along(sample_ids)) {
    out[[i]] <- read_cx_counts(manifest$raw_path[i], context)
  }
  out
}

safe_col <- function(df, names, default = NA) {
  for (name in names) {
    if (name %in% colnames(df)) {
      return(df[[name]])
    }
  }
  rep(default, nrow(df))
}

write_normalized <- function(df, out_path, caller, context) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(df) || nrow(df) == 0) {
    empty <- data.frame(
      dmr_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      context = character(),
      delta = numeric(),
      p_value = numeric(),
      q_value = numeric(),
      n_sites = integer(),
      source_caller = character()
    )
    data.table::fwrite(empty, out_path, sep = "\t")
    return()
  }
  df$context <- context
  df$source_caller <- caller
  data.table::fwrite(as.data.table(df), out_path, sep = "\t")
}

merge_bsseq_counts <- function(dat_list) {
  merged <- NULL
  for (sample_id in names(dat_list)) {
    x <- data.table::as.data.table(dat_list[[sample_id]])
    data.table::setnames(x, c("N", "X"), c(paste0("Cov_", sample_id), paste0("M_", sample_id)))
    if (is.null(merged)) {
      merged <- x
    } else {
      merged <- merge(merged, x, by = c("chr", "pos"), all = TRUE)
    }
  }
  setorder(merged, chr, pos)
  merged
}

make_bsseq <- function(dat_list) {
  if (!require_pkg("bsseq", "bsseq")) {
    return(NULL)
  }
  merged <- merge_bsseq_counts(dat_list)
  m_cols <- paste0("M_", sample_ids)
  cov_cols <- paste0("Cov_", sample_ids)
  m <- as.matrix(merged[, ..m_cols])
  cov <- as.matrix(merged[, ..cov_cols])
  colnames(m) <- sample_ids
  colnames(cov) <- sample_ids
  m[is.na(m)] <- 0
  cov[is.na(cov)] <- 0
  bs <- bsseq::BSseq(chr = merged$chr, pos = merged$pos, M = m, Cov = cov, sampleNames = sample_ids)
  Biobase::pData(bs)$condition <- factor(manifest$condition, levels = c("control", "treatment"))
  bs
}

filter_bsseq_for_dmrseq <- function(bs, context) {
  cov <- bsseq::getCoverage(bs, type = "Cov")
  control_cols <- colnames(cov) %in% control_ids
  treatment_cols <- colnames(cov) %in% treatment_ids
  if (!any(control_cols) || !any(treatment_cols)) {
    add_status("dmrseq", context, "skipped", "could not map control/treatment samples in BSseq coverage matrix")
    return(NULL)
  }
  keep <- rowSums(cov[, control_cols, drop = FALSE] > 0, na.rm = TRUE) > 0 &
    rowSums(cov[, treatment_cols, drop = FALSE] > 0, na.rm = TRUE) > 0
  n_removed <- length(keep) - sum(keep)
  if (sum(keep) < min_sites) {
    add_status(
      "dmrseq",
      context,
      "skipped",
      paste0("only ", sum(keep), " loci remain after condition-coverage filter; removed ", n_removed)
    )
    return(NULL)
  }
  if (n_removed > 0) {
    add_status(
      "dmrseq",
      context,
      "filtered",
      paste0("removed ", n_removed, " loci with zero coverage in all samples of at least one condition")
    )
    write_status()
  }
  bs[keep, ]
}

run_dss <- function(context, dat_list) {
  if (!require_pkg("DSS", "DSS", context)) {
    return(NULL)
  }
  suppressPackageStartupMessages(library(DSS))
  add_status("DSS", context, "started", "running DSS::DMLtest and DSS::callDMR")
  write_status()
  out_dir <- file.path(out_root, "dss")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bs <- DSS::makeBSseqData(dat_list, sample_ids)
  dml <- DSS::DMLtest(bs, group1 = control_ids, group2 = treatment_ids, smoothing = TRUE)
  data.table::fwrite(as.data.table(dml), file.path(out_dir, paste0("dss_dml_sites_", context, ".tsv")), sep = "\t")
  dmr <- DSS::callDMR(dml, p.threshold = q_threshold, minCG = min_sites, dis.merge = max_gap)
  dmr <- as.data.frame(dmr)
  out <- data.frame(
    dmr_id = paste0("DSS_", context, "_", seq_len(nrow(dmr))),
    chrom = safe_col(dmr, c("chr", "chrom")),
    start = safe_col(dmr, c("start")),
    end = safe_col(dmr, c("end")),
    delta = safe_col(dmr, c("diff.Methy", "diff", "delta")),
    p_value = safe_col(dmr, c("pval", "pvalue", "p_value")),
    q_value = safe_col(dmr, c("fdr", "FDR", "q_value", "qvalue")),
    n_sites = safe_col(dmr, c("nCG", "n_cytosines", "n_sites"))
  )
  write_normalized(out, file.path(out_dir, paste0("dss_dmrs_", context, ".tsv")), "DSS", context)
  add_status("DSS", context, "ok", paste0(nrow(out), " DMRs"))
  dml
}

fisher_pvalue <- function(pvalues) {
  pvalues <- pvalues[!is.na(pvalues)]
  if (length(pvalues) == 0) {
    return(NA_real_)
  }
  pvalues <- pmin(pmax(as.numeric(pvalues), 1e-300), 1)
  stats::pchisq(-2 * sum(log(pvalues)), df = 2 * length(pvalues), lower.tail = FALSE)
}

# Kost & McKay (2002) approximation of Cov(-2 ln p_i, -2 ln p_j) as a function
# of the correlation r between the two underlying test statistics. Only the
# positive-correlation branch matters: positive dependence inflates Var(X) and
# makes Fisher's method anti-conservative.
kost_cov <- function(r) {
  r * (3.263 + 0.710 * r + 0.027 * r * r)
}

# Genome-wide lag-1 autocorrelation of the per-site test statistics, used as a
# single mean correlation r for Brown's correction. Estimated once from all
# ordered sites; clamped to [0, 0.99] because only positive dependence is
# relevant here and r -> 1 is numerically unstable.
estimate_mean_correlation <- function(stat) {
  s <- as.numeric(stat)
  s <- s[is.finite(s)]
  if (length(s) < 3) {
    return(0)
  }
  r <- suppressWarnings(stats::cor(s[-length(s)], s[-1]))
  if (!is.finite(r)) {
    return(0)
  }
  max(0, min(r, 0.99))
}

# Brown's correction: approximate the dependent Fisher statistic X = -2 sum ln p
# by a scaled chi-square c * chi^2_f, matching the first two moments. With a
# common pairwise correlation r: E[X] = 2k, Var[X] = 4k + 2*C(k,2)*kost_cov(r),
# c = Var/(2 E), f = 2 E^2 / Var, and P = P(chi^2_f >= X / c). For r = 0 (or
# k = 1) this collapses back to the classical chi^2_{2k} Fisher p-value.
brown_pvalue <- function(pvalues, r_mean) {
  pvalues <- pvalues[!is.na(pvalues)]
  k <- length(pvalues)
  if (k == 0) {
    return(NA_real_)
  }
  pvalues <- pmin(pmax(as.numeric(pvalues), 1e-300), 1)
  x_stat <- -2 * sum(log(pvalues))
  if (k == 1 || is.na(r_mean) || r_mean <= 0) {
    return(stats::pchisq(x_stat, df = 2 * k, lower.tail = FALSE))
  }
  expected <- 2 * k
  variance <- 4 * k + 2 * choose(k, 2) * kost_cov(r_mean)
  c_scale <- variance / (2 * expected)
  f_dof <- 2 * expected^2 / variance
  stats::pchisq(x_stat / c_scale, df = f_dof, lower.tail = FALSE)
}

bh_q <- function(pvalues) {
  stats::p.adjust(pvalues, method = "BH")
}

run_combp_like <- function(context, dml) {
  if (is.null(dml) || !("comb-p" %in% callers)) {
    return()
  }
  add_status("comb-p", context, "started", "merging DSS per-site p-values")
  write_status()
  out_dir <- file.path(out_root, "combp")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  x <- as.data.table(dml)
  p_col <- intersect(c("pval", "pvalue", "p_value"), names(x))[1]
  delta_col <- intersect(c("diff", "diff.Methy", "delta"), names(x))[1]
  if (is.na(p_col) || !("chr" %in% names(x)) || !("pos" %in% names(x))) {
    add_status("comb-p", context, "skipped", "DSS DML table lacks chr/pos/p-value columns")
    return()
  }
  # Estimate the genome-wide correlation of the per-site statistics BEFORE
  # filtering to significant sites, so Brown's correction reflects the true
  # spatial dependence of the test, not just the selected tail.
  stat_col <- intersect(c("stat", "tstat", "statistic"), names(x))[1]
  full_setorder <- data.table::copy(x)
  data.table::setorder(full_setorder, chr, pos)
  r_mean <- if (!is.na(stat_col)) estimate_mean_correlation(full_setorder[[stat_col]]) else 0
  x[, p_for_merge := as.numeric(get(p_col))]
  x <- x[!is.na(p_for_merge) & p_for_merge <= q_threshold]
  if (nrow(x) == 0) {
    write_normalized(NULL, file.path(out_dir, paste0("combp_dmrs_", context, ".tsv")), "comb-p", context)
    add_status("comb-p", context, "ok", "0 regions")
    return()
  }
  setorder(x, chr, pos)
  rows <- list()
  for (chrom in unique(x$chr)) {
    sub <- x[chr == chrom]
    group_id <- cumsum(c(TRUE, diff(sub$pos) > max_gap))
    split_sub <- split(sub, group_id)
    for (g in split_sub) {
      if (nrow(g) < min_sites) {
        next
      }
      delta <- if (!is.na(delta_col)) mean(as.numeric(g[[delta_col]]), na.rm = TRUE) else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        chrom = chrom,
        start = min(g$pos),
        end = max(g$pos) + 1,
        delta = delta,
        p_value = brown_pvalue(g$p_for_merge, r_mean),
        n_sites = nrow(g)
      )
    }
  }
  out <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
  if (nrow(out) > 0) {
    out$q_value <- bh_q(out$p_value)
    out <- out[out$q_value <= q_threshold, , drop = FALSE]
    out$dmr_id <- paste0("combp_", context, "_", seq_len(nrow(out)))
  }
  write_normalized(out, file.path(out_dir, paste0("combp_dmrs_", context, ".tsv")), "comb-p", context)
  add_status(
    "comb-p",
    context,
    "ok",
    paste0(nrow(out), " regions; Brown-corrected Fisher (mean r=", round(r_mean, 3), ")")
  )
}

run_dmrseq <- function(context, bs) {
  if (!("dmrseq" %in% callers)) {
    return()
  }
  if (is.null(bs) || !require_pkg("dmrseq", "dmrseq", context)) {
    return()
  }
  add_status("dmrseq", context, "started", "running dmrseq::dmrseq")
  write_status()
  out_dir <- file.path(out_root, "dmrseq")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bs <- filter_bsseq_for_dmrseq(bs, context)
  if (is.null(bs)) {
    write_normalized(NULL, file.path(out_dir, paste0("dmrseq_dmrs_", context, ".tsv")), "dmrseq", context)
    return()
  }
  dmrseq_args <- list(bs = bs, testCovariate = "condition", cutoff = q_threshold)
  if (dmrseq_permutations > 0) {
    dmrseq_args$maxPerms <- dmrseq_permutations
  }
  dmr <- do.call(dmrseq::dmrseq, dmrseq_args)
  dmr <- as.data.frame(dmr)
  out <- data.frame(
    dmr_id = paste0("dmrseq_", context, "_", seq_len(nrow(dmr))),
    chrom = safe_col(dmr, c("chr", "chrom")),
    start = safe_col(dmr, c("start")),
    end = safe_col(dmr, c("end")),
    delta = safe_col(dmr, c("beta", "delta")),
    p_value = safe_col(dmr, c("pval", "p_value")),
    q_value = safe_col(dmr, c("qval", "q_value")),
    n_sites = safe_col(dmr, c("L", "n", "n_sites"))
  )
  write_normalized(out, file.path(out_dir, paste0("dmrseq_dmrs_", context, ".tsv")), "dmrseq", context)
  add_status("dmrseq", context, "ok", paste0(nrow(out), " DMRs"))
}

run_bsmooth <- function(context, bs) {
  if (!("BSmooth" %in% callers)) {
    return()
  }
  if (is.null(bs)) {
    return()
  }
  add_status("BSmooth", context, "started", "running bsseq::BSmooth and dmrFinder")
  write_status()
  out_dir <- file.path(out_root, "bsmooth")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  smoothed <- bsseq::BSmooth(bs)
  tstat <- bsseq::BSmooth.tstat(smoothed, group1 = control_ids, group2 = treatment_ids)
  dmr <- as.data.frame(bsseq::dmrFinder(tstat))
  out <- data.frame(
    dmr_id = paste0("BSmooth_", context, "_", seq_len(nrow(dmr))),
    chrom = safe_col(dmr, c("chr", "chrom")),
    start = safe_col(dmr, c("start")),
    end = safe_col(dmr, c("end")),
    delta = safe_col(dmr, c("areaStat", "meanDiff", "delta")),
    p_value = NA_real_,
    q_value = NA_real_,
    n_sites = safe_col(dmr, c("n", "n_sites"))
  )
  write_normalized(out, file.path(out_dir, paste0("bsmooth_dmrs_", context, ".tsv")), "BSmooth", context)
  add_status("BSmooth", context, "ok", paste0(nrow(out), " DMRs"))
}

write_methylkit_coverage <- function(context, dat_list) {
  temp_dir <- file.path(out_root, "methylkit_input", context)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character(length(sample_ids))
  for (i in seq_along(sample_ids)) {
    x <- as.data.table(dat_list[[sample_ids[i]]])
    out <- data.table(
      chr = x$chr,
      start = x$pos,
      end = x$pos,
      methylation_percent = ifelse(x$N > 0, 100 * x$X / x$N, 0),
      count_methylated = x$X,
      count_unmethylated = x$N - x$X
    )
    paths[i] <- file.path(temp_dir, paste0(sample_ids[i], ".bismark_coverage.tsv"))
    fwrite(out, paths[i], sep = "\t", col.names = FALSE)
  }
  paths
}

run_methylkit <- function(context, dat_list) {
  if (!("methylKit" %in% callers)) {
    return()
  }
  if (!require_pkg("methylKit", "methylKit", context)) {
    return()
  }
  add_status("methylKit", context, "started", "running methylKit regional differential methylation")
  write_status()
  out_dir <- file.path(out_root, "methylkit")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  coverage_paths <- write_methylkit_coverage(context, dat_list)
  raw <- methylKit::methRead(
    location = as.list(coverage_paths),
    sample.id = as.list(sample_ids),
    assembly = "SL3.0",
    treatment = treatment_vector,
    context = context,
    pipeline = "bismarkCoverage"
  )
  if (exists("tileMethylCounts", where = asNamespace("methylKit"), mode = "function")) {
    raw <- methylKit::tileMethylCounts(raw, win.size = window_size, step.size = step_size)
  }
  united <- methylKit::unite(raw, destrand = FALSE)
  diff <- methylKit::calculateDiffMeth(united)
  dmr <- as.data.frame(methylKit::getMethylDiff(diff, difference = methylkit_min_diff, qvalue = q_threshold))
  out <- data.frame(
    dmr_id = paste0("methylKit_", context, "_", seq_len(nrow(dmr))),
    chrom = safe_col(dmr, c("chr", "chrom")),
    start = safe_col(dmr, c("start")),
    end = safe_col(dmr, c("end")),
    delta = safe_col(dmr, c("meth.diff", "meth_diff")) / 100,
    p_value = safe_col(dmr, c("pvalue", "p_value")),
    q_value = safe_col(dmr, c("qvalue", "q_value")),
    n_sites = NA_integer_
  )
  write_normalized(out, file.path(out_dir, paste0("methylkit_dmrs_", context, ".tsv")), "methylKit", context)
  add_status("methylKit", context, "ok", paste0(nrow(out), " regions"))
}

write_metilene_input <- function(context, dat_list) {
  if (!("metilene" %in% callers)) {
    return()
  }
  add_status("metilene", context, "started", "writing metilene input matrix")
  write_status()
  out_dir <- file.path(out_root, "metilene_input")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  merged <- merge_bsseq_counts(dat_list)
  out <- data.table(chr = merged$chr, pos = merged$pos)
  for (sample_id in sample_ids) {
    m <- merged[[paste0("M_", sample_id)]]
    cov <- merged[[paste0("Cov_", sample_id)]]
    sample_condition <- as.character(manifest[["condition"]][match(sample_id, manifest[["sample_id"]])])
    metilene_sample_id <- paste(sample_condition, sample_id, sep = "_")
    out[[metilene_sample_id]] <- ifelse(!is.na(cov) & cov > 0, m / cov, NA)
  }
  control_cols <- paste("control", control_ids, sep = "_")
  treatment_cols <- paste("treatment", treatment_ids, sep = "_")
  control_cols <- intersect(control_cols, names(out))
  treatment_cols <- intersect(treatment_cols, names(out))
  if (length(control_cols) == 0 || length(treatment_cols) == 0) {
    write_normalized(NULL, file.path(out_root, "metilene", paste0("metilene_dmrs_", context, ".tsv")), "metilene", context)
    add_status("metilene", context, "skipped", "missing control or treatment columns in metilene matrix")
    write_status()
    return()
  }
  out[, n_control_values := rowSums(!is.na(.SD)), .SDcols = control_cols]
  out[, n_treatment_values := rowSums(!is.na(.SD)), .SDcols = treatment_cols]
  out <- out[n_control_values > 0 & n_treatment_values > 0]
  out[, c("n_control_values", "n_treatment_values") := NULL]
  data.table::setorder(out, chr, pos)
  input_path <- file.path(out_dir, paste0("metilene_", context, ".tsv"))
  if (nrow(out) < min_sites) {
    data.table::fwrite(out, input_path, sep = "\t", na = "NA")
    write_normalized(NULL, file.path(out_root, "metilene", paste0("metilene_dmrs_", context, ".tsv")), "metilene", context)
    add_status("metilene", context, "skipped", paste0("only ", nrow(out), " usable sites after group coverage filtering"))
    write_status()
    return()
  }
  fwrite(out, input_path, sep = "\t", na = "NA")
  run_metilene(context, input_path)
}

normalize_metilene_output <- function(raw_path, normalized_path, context) {
  if (!file.exists(raw_path) || file.info(raw_path)$size == 0) {
    write_normalized(NULL, normalized_path, "metilene", context)
    return(0L)
  }
  raw <- data.table::fread(raw_path, header = FALSE)
  if (nrow(raw) == 0) {
    write_normalized(NULL, normalized_path, "metilene", context)
    return(0L)
  }
  expected <- c("chrom", "start", "end", "q_value", "delta", "n_sites", "p_value", "p2_value", "m1", "m2")
  data.table::setnames(raw, seq_len(min(length(expected), ncol(raw))), expected[seq_len(min(length(expected), ncol(raw)))])
  out <- data.frame(
    dmr_id = paste0("metilene_", context, "_", seq_len(nrow(raw))),
    chrom = safe_col(raw, c("chrom")),
    start = safe_col(raw, c("start")),
    end = safe_col(raw, c("end")),
    delta = safe_col(raw, c("delta")),
    p_value = safe_col(raw, c("p_value")),
    q_value = safe_col(raw, c("q_value")),
    n_sites = safe_col(raw, c("n_sites"))
  )
  write_normalized(out, normalized_path, "metilene", context)
  nrow(out)
}

run_metilene <- function(context, input_path) {
  out_dir <- file.path(out_root, "metilene")
  input_dir <- file.path(out_root, "metilene_input")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  raw_path <- file.path(out_dir, paste0("metilene_dmrs_", context, ".raw.tsv"))
  normalized_path <- file.path(out_dir, paste0("metilene_dmrs_", context, ".tsv"))
  stderr_path <- file.path(out_dir, paste0("metilene_", context, ".stderr.log"))
  exe <- if (nchar(metilene_bin) > 0) metilene_bin else Sys.which("metilene")
  metilene_args <- c(
    "-a", "control",
    "-b", "treatment",
    "-M", as.character(max_gap),
    "-m", as.character(min_sites),
    "-c", "2",
    "-X", "1",
    "-Y", "1",
    input_path
  )
  command <- paste(
    ifelse(nchar(exe) > 0, exe, "metilene"),
    paste(shQuote(metilene_args), collapse = " "),
    ">",
    shQuote(raw_path)
  )
  writeLines(command, file.path(input_dir, paste0("run_metilene_", context, ".cmd.txt")))
  if (nchar(exe) == 0) {
    write_normalized(NULL, normalized_path, "metilene", context)
    add_status("metilene", context, "input_ready", "metilene input matrix written; binary not found on PATH and --metilene was not supplied")
    write_status()
    return()
  }
  exit_code <- system2(exe, args = metilene_args, stdout = raw_path, stderr = stderr_path)
  if (!identical(as.integer(exit_code), 0L)) {
    add_status("metilene", context, "error", paste0("metilene exited with code ", exit_code, "; see ", stderr_path))
    write_status()
    return()
  }
  n_regions <- normalize_metilene_output(raw_path, normalized_path, context)
  add_status("metilene", context, "ok", paste0(n_regions, " DMRs"))
  write_status()
}

if ("DMRcate" %in% callers) {
  add_status("DMRcate", "", "not_implemented", "DMRcate execution is not implemented in this runner yet; import/harmonization is supported.")
}

for (ctx in contexts) {
  message("=== context ", ctx, " ===")
  dat_list <- read_context_counts(ctx)
  dml <- NULL
  if ("DSS" %in% callers || "comb-p" %in% callers) {
    dml <- tryCatch(run_dss(ctx, dat_list), error = function(e) {
      add_status("DSS", ctx, "error", conditionMessage(e))
      NULL
    })
  }
  tryCatch(run_combp_like(ctx, dml), error = function(e) add_status("comb-p", ctx, "error", conditionMessage(e)))
  bs <- NULL
  if ("dmrseq" %in% callers || "BSmooth" %in% callers) {
    bs <- tryCatch(make_bsseq(dat_list), error = function(e) {
      add_status("bsseq", ctx, "error", conditionMessage(e))
      NULL
    })
  }
  tryCatch(run_dmrseq(ctx, bs), error = function(e) add_status("dmrseq", ctx, "error", conditionMessage(e)))
  tryCatch(run_bsmooth(ctx, bs), error = function(e) add_status("BSmooth", ctx, "error", conditionMessage(e)))
  tryCatch(run_methylkit(ctx, dat_list), error = function(e) add_status("methylKit", ctx, "error", conditionMessage(e)))
  tryCatch(write_metilene_input(ctx, dat_list), error = function(e) add_status("metilene", ctx, "error", conditionMessage(e)))
  write_status()
}

write_status()
