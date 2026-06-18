#!/usr/bin/env Rscript
# Build ComplexUpset plots of DMR caller-support intersections.
#
# Input: a membership TSV with one row per consensus region, 0/1 columns per
# caller (the "sets"), and optional `region_id` and `context` columns.
# Output: upset_caller_support_overall.{pdf,png} and one plot per context.

arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name)
  idx <- which(args == key)
  if (length(idx) == 0 || idx[length(idx)] >= length(args)) {
    return(default)
  }
  args[idx[length(idx)] + 1]
}

membership_path <- arg_value("membership")
out_dir <- arg_value("out-dir", ".")
colors_path <- arg_value("colors")
min_size <- as.integer(arg_value("min-size", "1"))
if (is.na(min_size) || min_size < 0) {
  min_size <- 1
}

# Optional caller -> hex color map (matches the Python report palette).
caller_colors <- list()
if (!is.null(colors_path) && file.exists(colors_path)) {
  color_table <- utils::read.delim(colors_path, stringsAsFactors = FALSE)
  if (all(c("caller", "color") %in% names(color_table))) {
    for (i in seq_len(nrow(color_table))) {
      caller_colors[[as.character(color_table$caller[i])]] <- as.character(color_table$color[i])
    }
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(out_dir, "upset_status.txt")
write_status <- function(text) {
  cat(text, "\n", file = status_path, append = FALSE)
}

if (is.null(membership_path) || !file.exists(membership_path)) {
  write_status(paste0("membership table not found: ", membership_path))
  quit(status = 3)
}

ok <- requireNamespace("ggplot2", quietly = TRUE) &&
  requireNamespace("ComplexUpset", quietly = TRUE)
if (!ok) {
  write_status("missing R packages: install ggplot2 and ComplexUpset")
  quit(status = 4)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexUpset)
})

df <- utils::read.delim(membership_path, check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- intersect(c("region_id", "context"), names(df))
set_cols_all <- setdiff(names(df), meta_cols)
if (length(set_cols_all) == 0) {
  write_status("no caller set columns found in membership table")
  quit(status = 5)
}
df[set_cols_all] <- lapply(df[set_cols_all], function(x) as.logical(as.integer(x)))

render_upset <- function(data, label) {
  set_cols <- set_cols_all[vapply(set_cols_all, function(col) any(data[[col]], na.rm = TRUE), logical(1))]
  if (length(set_cols) < 1 || nrow(data) == 0) {
    return(FALSE)
  }
  # One colored query per caller: colors the set-size bar and the matrix points
  # for that caller, using the shared report palette.
  queries <- list()
  for (col in set_cols) {
    if (!is.null(caller_colors[[col]])) {
      queries[[length(queries) + 1]] <- ComplexUpset::upset_query(
        set = col,
        fill = caller_colors[[col]],
        color = caller_colors[[col]]
      )
    }
  }
  # Render the intersection-size counts as vertical labels placed *above* every
  # bar (bar_number_threshold = Inf), with extra headroom on the y axis so the
  # numbers never overlap the bars or each other, even for the long tail of
  # small intersections.
  size_annotation <- ComplexUpset::intersection_size(
    bar_number_threshold = Inf,
    text = list(angle = 90, vjust = 0.5, hjust = -0.1, size = 2.6),
    text_colors = c(on_background = "black", on_bar = "black")
  ) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.3)))
  plot <- ComplexUpset::upset(
    data,
    intersect = set_cols,
    name = "DMR caller support",
    min_size = min_size,
    sort_intersections_by = "cardinality",
    queries = queries,
    base_annotations = list("Intersection size" = size_annotation),
    set_sizes = ComplexUpset::upset_set_size()
  ) + ggplot2::labs(title = paste0("DMR caller support intersections: ", label))
  base <- file.path(out_dir, paste0("upset_caller_support_", label))
  # Widen the figure with the number of intersection bars actually shown so the
  # columns (and their vertical count labels) stay legible instead of crowding.
  codes <- apply(data[set_cols], 1, function(r) paste(as.integer(r), collapse = ""))
  n_intersections <- sum(table(codes) >= min_size)
  width <- max(7, min(26, n_intersections * 0.34), length(set_cols) * 1.6)
  tryCatch(
    {
      ggplot2::ggsave(paste0(base, ".pdf"), plot, width = width, height = 6, limitsize = FALSE)
      ggplot2::ggsave(paste0(base, ".png"), plot, width = width, height = 6, dpi = 200, limitsize = FALSE)
      TRUE
    },
    error = function(e) {
      cat("failed to render ", label, ": ", conditionMessage(e), "\n", file = status_path, append = TRUE)
      FALSE
    }
  )
}

built <- character()
if (render_upset(df, "overall")) {
  built <- c(built, "overall")
}
if ("context" %in% names(df)) {
  for (ctx in sort(unique(df$context))) {
    ctx_label <- gsub("[^A-Za-z0-9_]+", "_", as.character(ctx))
    if (render_upset(df[df$context == ctx, , drop = FALSE], ctx_label)) {
      built <- c(built, ctx_label)
    }
  }
}

write_status(paste0("built: ", paste(built, collapse = ", ")))
quit(status = 0)
