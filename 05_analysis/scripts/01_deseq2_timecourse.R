#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_pkgs <- c("DESeq2", "apeglm")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    sprintf(
      "Missing required R packages: %s\nInstall them before running.",
      paste(missing_pkgs, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
})

control_genes_default <- c(
  "FBgn0000158", # bam
  "FBgn0000146", # aub
  "FBgn0032473", # kmg
  "FBgn0004372", # aly
  "FBgn0260400"  # elav
)

parse_args <- function(args) {
  opts <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop(sprintf("Invalid argument format: %s", arg), call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    val <- sub("^--[^=]+=(.*)$", "\\1", arg)
    opts[[key]] <- val
  }
  opts
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

res_to_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  df$Geneid <- rownames(df)
  cols <- c("Geneid", setdiff(colnames(df), "Geneid"))
  df <- df[, cols, drop = FALSE]
  if ("padj" %in% colnames(df)) {
    df <- df[order(df$padj, na.last = TRUE), , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, mustWork = FALSE)
  repeat {
    has_samples <- file.exists(file.path(cur, "00_admin", "metadata", "samples.tsv"))
    has_counts <- file.exists(file.path(cur, "02_work", "06_counts", "final", "gene_counts.matrix.tsv"))
    if (has_samples && has_counts) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      break
    }
    cur <- parent
  }
  normalizePath(start_dir, mustWork = FALSE)
}

extract_gtf_from_featurecounts_header <- function(featurecounts_path) {
  if (!file.exists(featurecounts_path)) {
    return(NA_character_)
  }
  first <- readLines(featurecounts_path, n = 1, warn = FALSE)
  if (length(first) == 0) {
    return(NA_character_)
  }
  m <- regexec('"-a" "([^"]+)"', first)
  reg <- regmatches(first, m)[[1]]
  if (length(reg) >= 2) {
    return(reg[2])
  }
  NA_character_
}

infer_gtf_path <- function(project_root, featurecounts_path, gtf_arg = NULL) {
  candidates <- character(0)

  if (!is.null(gtf_arg) && nzchar(gtf_arg)) {
    candidates <- c(candidates, normalizePath(gtf_arg, mustWork = FALSE))
  }

  from_header <- extract_gtf_from_featurecounts_header(featurecounts_path)
  if (!is.na(from_header)) {
    candidates <- c(candidates, normalizePath(from_header, mustWork = FALSE))
    candidates <- c(candidates, file.path(project_root, "02_work", "04_ref", "reference", basename(from_header)))
  }

  local_ref_dir <- file.path(project_root, "02_work", "04_ref", "reference")
  if (dir.exists(local_ref_dir)) {
    local_gtf <- list.files(local_ref_dir, pattern = "\\.gtf(\\.gz)?$", full.names = TRUE)
    if (length(local_gtf) > 0) {
      candidates <- c(candidates, local_gtf)
    }
  }

  if (length(candidates) == 0) {
    return(NA_character_)
  }

  candidates <- unique(candidates)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[1])
  }

  candidates[1]
}

extract_rrna_gene_ids_from_gtf <- function(gtf_path) {
  if (is.na(gtf_path) || !nzchar(gtf_path) || !file.exists(gtf_path)) {
    return(character(0))
  }

  con <- if (grepl("\\.gz$", gtf_path)) gzfile(gtf_path, "rt") else file(gtf_path, "rt")
  on.exit(close(con), add = TRUE)

  lines <- readLines(con, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (length(lines) == 0) {
    return(character(0))
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)
  attr <- vapply(fields, function(x) if (length(x) >= 9) x[9] else "", character(1))

  gene_id <- ifelse(
    grepl('gene_id "', attr),
    sub('.*gene_id "([^"]+)".*', '\\1', attr),
    NA_character_
  )

  biotype <- rep("", length(attr))
  for (k in c("gene_biotype", "gene_type", "transcript_biotype", "transcript_type")) {
    idx <- biotype == "" & grepl(paste0(k, ' "'), attr)
    if (any(idx)) {
      biotype[idx] <- sub(paste0('.*', k, ' "([^"]+)".*'), '\\1', attr[idx])
    }
  }

  gene_name <- ifelse(
    grepl('gene_name "', attr),
    sub('.*gene_name "([^"]+)".*', '\\1', attr),
    ""
  )

  rrna_idx <- !is.na(gene_id) & (
    grepl("rrna", biotype, ignore.case = TRUE) |
      grepl("rrna", gene_name, ignore.case = TRUE)
  )

  unique(gene_id[rrna_idx])
}

is_mito_chr <- function(chr_vec) {
  sapply(strsplit(as.character(chr_vec), ";", fixed = TRUE), function(tokens) {
    toks <- trimws(tokens)
    any(grepl("mitochond|^MT$|^M$|^chrM$|^chrMT$", toks, ignore.case = TRUE))
  })
}

compute_sample_qc_metrics <- function(featurecounts_path, fc_summary_path, coldata, outdir, gtf_path = NA_character_) {
  if (!file.exists(featurecounts_path)) {
    stop(sprintf("featureCounts table not found: %s", featurecounts_path), call. = FALSE)
  }
  if (!file.exists(fc_summary_path)) {
    stop(sprintf("featureCounts summary not found: %s", fc_summary_path), call. = FALSE)
  }

  fc <- read.delim(featurecounts_path, comment.char = "#", check.names = FALSE)
  if (ncol(fc) < 7) {
    stop("featureCounts table has unexpected format.", call. = FALSE)
  }

  sample_cols <- colnames(fc)[7:ncol(fc)]
  sample_ids_fc <- sub("\\.sorted\\.bam$", "", basename(sample_cols))

  fc_counts <- as.matrix(fc[, sample_cols, drop = FALSE])
  storage.mode(fc_counts) <- "numeric"
  colnames(fc_counts) <- sample_ids_fc

  sm <- read.delim(fc_summary_path, check.names = FALSE)
  status <- sm[[1]]
  sm_mat <- as.matrix(sm[, -1, drop = FALSE])
  storage.mode(sm_mat) <- "numeric"
  rownames(sm_mat) <- status
  colnames(sm_mat) <- sub("\\.sorted\\.bam$", "", basename(colnames(sm_mat)))

  required_samples <- rownames(coldata)
  missing_fc <- setdiff(required_samples, colnames(fc_counts))
  missing_sm <- setdiff(required_samples, colnames(sm_mat))
  if (length(missing_fc) > 0 || length(missing_sm) > 0) {
    stop(
      sprintf(
        "Sample mismatch in featureCounts inputs. Missing in counts: [%s]; missing in summary: [%s]",
        paste(missing_fc, collapse = ","),
        paste(missing_sm, collapse = ",")
      ),
      call. = FALSE
    )
  }

  samples <- required_samples
  assigned_reads <- sm_mat["Assigned", samples]
  total_reads <- colSums(sm_mat[, samples, drop = FALSE], na.rm = TRUE)
  assigned_pct <- ifelse(total_reads > 0, 100 * assigned_reads / total_reads, NA_real_)

  mito_mask <- is_mito_chr(fc$Chr)
  mito_reads <- if (any(mito_mask)) {
    colSums(fc_counts[mito_mask, samples, drop = FALSE], na.rm = TRUE)
  } else {
    setNames(rep(0, length(samples)), samples)
  }

  rrna_genes <- extract_rrna_gene_ids_from_gtf(gtf_path)
  if (length(rrna_genes) > 0) {
    rrna_mask <- fc$Geneid %in% rrna_genes
    rrna_method <- "GTF biotype/gene_name (rRNA)"
  } else {
    rrna_mask <- grepl("rrna", fc$Geneid, ignore.case = TRUE)
    rrna_method <- "Geneid regex fallback (rrna)"
  }

  rrna_reads <- if (any(rrna_mask)) {
    colSums(fc_counts[rrna_mask, samples, drop = FALSE], na.rm = TRUE)
  } else {
    setNames(rep(0, length(samples)), samples)
  }

  rrna_pct_total <- ifelse(total_reads > 0, 100 * rrna_reads / total_reads, NA_real_)
  mito_pct_total <- ifelse(total_reads > 0, 100 * mito_reads / total_reads, NA_real_)
  genes_detected <- colSums(fc_counts[, samples, drop = FALSE] > 0, na.rm = TRUE)

  out <- data.frame(
    sample = samples,
    timepoint = as.character(coldata[samples, "timepoint"]),
    lane = as.character(coldata[samples, "lane"]),
    total_reads_featurecounts = as.numeric(total_reads),
    assigned_reads = as.numeric(assigned_reads),
    assigned_pct = round(as.numeric(assigned_pct), 4),
    rrna_reads = as.numeric(rrna_reads),
    rrna_pct_total_reads = round(as.numeric(rrna_pct_total), 4),
    mito_reads = as.numeric(mito_reads),
    mito_pct_total_reads = round(as.numeric(mito_pct_total), 4),
    genes_detected = as.integer(genes_detected),
    stringsAsFactors = FALSE
  )

  write_tsv(out, file.path(outdir, "sample_qc_metrics.tsv"))
  notes <- data.frame(
    key = c("featurecounts_table", "featurecounts_summary", "gtf_used", "rrna_method", "mito_definition"),
    value = c(
      featurecounts_path,
      fc_summary_path,
      ifelse(is.na(gtf_path), "NA", gtf_path),
      rrna_method,
      "Chr contains mitochond/MT/M/chrM/chrMT"
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(notes, file.path(outdir, "sample_qc_metrics_notes.tsv"))

  invisible(out)
}

plot_volcano <- function(df, title, png_path, pdf_path, padj_thr = 0.05, lfc_thr = 1) {
  if (!all(c("log2FoldChange", "padj") %in% colnames(df))) {
    return(invisible(FALSE))
  }

  x <- df$log2FoldChange
  y <- -log10(pmax(df$padj, .Machine$double.xmin))
  y[is.na(y)] <- 0
  sig <- !is.na(df$padj) & df$padj < padj_thr & !is.na(x) & abs(x) >= lfc_thr

  x_range <- suppressWarnings(quantile(abs(x[is.finite(x)]), 0.99, na.rm = TRUE))
  if (!is.finite(x_range) || x_range < 2) {
    x_range <- 2
  }
  y_max <- suppressWarnings(max(y, na.rm = TRUE))
  if (!is.finite(y_max) || y_max < 5) {
    y_max <- 5
  }

  draw_plot <- function() {
    plot(
      x,
      y,
      pch = 16,
      cex = 0.5,
      col = ifelse(sig, "firebrick3", "grey75"),
      xlab = "log2FoldChange",
      ylab = "-log10(padj)",
      main = title,
      xlim = c(-x_range, x_range),
      ylim = c(0, y_max * 1.05)
    )
    abline(v = c(-lfc_thr, lfc_thr), lty = 2)
    abline(h = -log10(padj_thr), lty = 2)
  }

  png(filename = png_path, width = 1800, height = 1400, res = 220)
  draw_plot()
  dev.off()

  pdf(file = pdf_path, width = 8, height = 6)
  draw_plot()
  dev.off()

  invisible(TRUE)
}

plot_pca <- function(vsd, coldata, prefix) {
  mat <- assay(vsd)
  pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_expl <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  pca_df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  pca_df$timepoint <- coldata[pca_df$sample, "timepoint"]
  pca_df$lane <- coldata[pca_df$sample, "lane"]

  draw_scatter <- function(color_by, cols, title_suffix, out_stem) {
    group_vals <- as.character(pca_df[[color_by]])
    point_cols <- cols[group_vals]

    draw <- function() {
      plot(
        pca_df$PC1,
        pca_df$PC2,
        pch = 19,
        col = point_cols,
        xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
        ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
        main = sprintf("PCA (%s)", title_suffix)
      )
      text(pca_df$PC1, pca_df$PC2, labels = pca_df$sample, pos = 3, cex = 0.65)
      legend("topright", legend = names(cols), col = cols, pch = 19, bty = "n")
    }

    png(file.path(prefix, paste0(out_stem, ".png")), width = 1800, height = 1400, res = 220)
    draw()
    dev.off()

    pdf(file.path(prefix, paste0(out_stem, ".pdf")), width = 8, height = 6)
    draw()
    dev.off()
  }

  tp_levels <- levels(coldata$timepoint)
  tp_palette <- c("0H" = "#1b9e77", "2H" = "#d95f02", "4H" = "#7570b3", "6H" = "#e7298a")
  tp_cols <- tp_palette[tp_levels]
  names(tp_cols) <- tp_levels

  lane_levels <- levels(coldata$lane)
  lane_palette <- c("L7" = "steelblue3", "L8" = "firebrick3", "NA" = "grey50")
  lane_cols <- lane_palette[lane_levels]
  names(lane_cols) <- lane_levels

  draw_scatter("timepoint", tp_cols, "colored by timepoint", "pca_timepoint")
  draw_scatter("lane", lane_cols, "colored by lane", "pca_lane")

  invisible(pca_df)
}

compute_pc1_loadings <- function(pca_obj) {
  load <- pca_obj$rotation[, 1]
  contrib <- (load^2 / sum(load^2)) * 100
  load_df <- data.frame(
    Geneid = names(load),
    PC1_loading = as.numeric(load),
    abs_loading = abs(as.numeric(load)),
    contribution_pct_PC1 = as.numeric(contrib),
    stringsAsFactors = FALSE
  )
  load_df[order(load_df$abs_loading, decreasing = TRUE), , drop = FALSE]
}

extract_go_term_genes <- function(gene_id_field) {
  if (is.na(gene_id_field) || !nzchar(gene_id_field)) {
    return(character(0))
  }
  toks <- trimws(unlist(strsplit(gene_id_field, "/", fixed = TRUE)))
  toks[nzchar(toks)]
}

run_quick_go <- function(top_loading_df, universe_genes, out_prefix,
                         plot_title = "Quick GO BP (top PC1 loading genes)") {
  status_file <- paste0(out_prefix, "_status.txt")

  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
    writeLines(
      "GO quick skipped: clusterProfiler and/or org.Dm.eg.db not installed.",
      con = status_file
    )
    return(invisible(NULL))
  }

  if (is.null(top_loading_df) || nrow(top_loading_df) == 0) {
    writeLines("GO quick skipped: empty loading table.", con = status_file)
    return(invisible(NULL))
  }

  genes <- unique(top_loading_df$Geneid[grepl("^FBgn", top_loading_df$Geneid)])
  universe <- unique(universe_genes[grepl("^FBgn", universe_genes)])

  if (length(genes) < 5 || length(universe) < 20) {
    writeLines(
      sprintf("GO quick skipped: not enough FlyBase IDs (genes=%d, universe=%d).", length(genes), length(universe)),
      con = status_file
    )
    return(invisible(NULL))
  }

  ego <- clusterProfiler::enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.Dm.eg.db::org.Dm.eg.db,
    keyType = "FLYBASE",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = FALSE
  )

  go_df <- as.data.frame(ego)
  if (nrow(go_df) == 0) {
    writeLines("GO quick returned no significant terms.", con = status_file)
    write_tsv(data.frame(), paste0(out_prefix, ".tsv"))
    return(invisible(NULL))
  }

  top_loading_df <- top_loading_df[!duplicated(top_loading_df$Geneid), , drop = FALSE]
  contrib_map <- setNames(top_loading_df$contribution_pct_PC1, top_loading_df$Geneid)
  total_top_contrib <- sum(top_loading_df$contribution_pct_PC1, na.rm = TRUE)

  term_genes <- lapply(go_df$geneID, extract_go_term_genes)
  go_df$pc1_contribution_pct <- vapply(
    term_genes,
    function(gs) sum(contrib_map[intersect(gs, names(contrib_map))], na.rm = TRUE),
    numeric(1)
  )
  go_df$pc1_contribution_pct_of_top50 <- ifelse(
    total_top_contrib > 0,
    100 * go_df$pc1_contribution_pct / total_top_contrib,
    NA_real_
  )
  go_df$n_top50_genes_in_term <- vapply(
    term_genes,
    function(gs) length(intersect(gs, names(contrib_map))),
    integer(1)
  )

  go_df <- go_df[order(go_df$p.adjust, -go_df$pc1_contribution_pct), , drop = FALSE]
  write_tsv(go_df, paste0(out_prefix, ".tsv"))

  top_go <- head(go_df, 15)
  labels <- substr(top_go$Description, 1, 80)
  vals <- top_go$pc1_contribution_pct

  draw <- function() {
    op <- par(mar = c(5, 15, 4, 2))
    on.exit(par(op), add = TRUE)
    barplot(
      rev(vals),
      names.arg = rev(labels),
      horiz = TRUE,
      las = 1,
      col = "steelblue3",
      border = NA,
      xlab = "Contribution to PC1 (%)",
      main = plot_title
    )
  }

  png(filename = paste0(out_prefix, ".png"), width = 2000, height = 1400, res = 220)
  draw()
  dev.off()

  pdf(file = paste0(out_prefix, ".pdf"), width = 10, height = 7)
  draw()
  dev.off()

  writeLines(
    sprintf("GO quick completed with %d significant terms.", nrow(go_df)),
    con = status_file
  )

  invisible(go_df)
}

plot_pca_global_de <- function(vsd, coldata, global_genes, prefix) {
  genes <- unique(global_genes)
  genes <- genes[genes %in% rownames(vsd)]
  if (length(genes) < 2) {
    return(invisible(NULL))
  }

  mat <- assay(vsd)[genes, , drop = FALSE]
  pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_expl <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  pca_df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
  pca_df$timepoint <- coldata[pca_df$sample, "timepoint"]

  tp_levels <- levels(coldata$timepoint)
  tp_palette <- c("0H" = "#1b9e77", "2H" = "#d95f02", "4H" = "#7570b3", "6H" = "#e7298a")
  tp_cols <- tp_palette[tp_levels]
  point_cols <- tp_cols[as.character(pca_df$timepoint)]

  draw <- function() {
    plot(
      pca_df$PC1,
      pca_df$PC2,
      pch = 19,
      col = point_cols,
      xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
      ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
      main = sprintf("PCA global DE genes (n=%d)", length(genes))
    )
    text(pca_df$PC1, pca_df$PC2, labels = pca_df$sample, pos = 3, cex = 0.65)
    legend("topright", legend = names(tp_cols), col = tp_cols, pch = 19, bty = "n")
  }

  png(file.path(prefix, "pca_global_de_timepoint.png"), width = 1800, height = 1400, res = 220)
  draw()
  dev.off()

  pdf(file.path(prefix, "pca_global_de_timepoint.pdf"), width = 8, height = 6)
  draw()
  dev.off()

  load_df <- compute_pc1_loadings(pca)
  top50 <- head(load_df, 50)

  write_tsv(load_df, file.path(prefix, "pc1_all_loadings_contrib_global_de_all_samples.tsv"))
  write_tsv(top50, file.path(prefix, "pc1_top50_loadings_contrib_global_de_all_samples.tsv"))

  run_quick_go(
    top_loading_df = top50,
    universe_genes = rownames(mat),
    out_prefix = file.path(prefix, "go_quick_pc1_top50_global_de_all_samples"),
    plot_title = "Quick GO BP (PC1 rare, all samples)"
  )

  invisible(list(load_df = load_df, top50 = top50))
}

plot_pca_global_de_excluded_pc1 <- function(vsd, coldata, global_genes, prefix,
                                            exclude_samples = c("HSF4_2H_3", "HSF4_4H_2", "HSF4_4H_3")) {
  genes <- unique(global_genes)
  genes <- genes[genes %in% rownames(vsd)]

  status_file <- file.path(prefix, "pca_global_de_excl_2H3_4H2_4H3_status.txt")
  if (length(genes) < 2) {
    writeLines("Skipped excluded-sample global-DE PCA: <2 global DE genes.", con = status_file)
    return(invisible(NULL))
  }

  keep_samples <- setdiff(colnames(vsd), exclude_samples)
  keep_samples <- keep_samples[keep_samples %in% rownames(coldata)]

  if (length(keep_samples) < 3) {
    writeLines("Skipped excluded-sample global-DE PCA: <3 samples after exclusion.", con = status_file)
    return(invisible(NULL))
  }

  mat <- assay(vsd)[genes, keep_samples, drop = FALSE]
  pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_expl <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  pca_df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  pca_df$timepoint <- coldata[pca_df$sample, "timepoint"]

  tp_levels <- levels(coldata$timepoint)
  tp_palette <- c("0H" = "#1b9e77", "2H" = "#d95f02", "4H" = "#7570b3", "6H" = "#e7298a")
  tp_cols <- tp_palette[tp_levels]
  point_cols <- tp_cols[as.character(pca_df$timepoint)]

  draw <- function() {
    plot(
      pca_df$PC1,
      pca_df$PC2,
      pch = 19,
      col = point_cols,
      xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
      ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
      main = sprintf("PCA global DE genes excl: 2H_3,4H_2,4H_3 (n=%d genes)", length(genes))
    )
    text(pca_df$PC1, pca_df$PC2, labels = pca_df$sample, pos = 3, cex = 0.7)
    legend("topright", legend = names(tp_cols), col = tp_cols, pch = 19, bty = "n")
  }

  png(file.path(prefix, "pca_global_de_timepoint_excl_2H3_4H2_4H3.png"), width = 1800, height = 1400, res = 220)
  draw()
  dev.off()

  pdf(file.path(prefix, "pca_global_de_timepoint_excl_2H3_4H2_4H3.pdf"), width = 8, height = 6)
  draw()
  dev.off()

  load_df <- compute_pc1_loadings(pca)
  top50 <- head(load_df, 50)

  write_tsv(load_df, file.path(prefix, "pc1_all_loadings_contrib_global_de_excl_2H3_4H2_4H3.tsv"))
  write_tsv(top50, file.path(prefix, "pc1_top50_loadings_contrib_global_de_excl_2H3_4H2_4H3.tsv"))

  run_quick_go(
    top_loading_df = top50,
    universe_genes = rownames(mat),
    out_prefix = file.path(prefix, "go_quick_pc1_top50_global_de_excl_2H3_4H2_4H3"),
    plot_title = "Quick GO BP (PC1 time, outliers removed)"
  )

  writeLines(
    sprintf("Excluded-sample PCA completed with %d genes and %d samples.", nrow(mat), ncol(mat)),
    con = status_file
  )

  invisible(list(load_df = load_df, top50 = top50))
}

is_gzipped_file <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  sig <- readBin(con, what = "raw", n = 2)
  length(sig) == 2 && as.integer(sig[1]) == 31 && as.integer(sig[2]) == 139
}

read_testis_specificity_db <- function(testis_db_path, target_genes = NULL) {
  out <- data.frame(
    Geneid = character(0),
    testis_mated_4day_tissue_specificity_index = numeric(0),
    stringsAsFactors = FALSE
  )

  if (is.na(testis_db_path) || !nzchar(testis_db_path) || !file.exists(testis_db_path)) {
    return(out)
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    return(out)
  }

  raw <- readxl::read_excel(testis_db_path, skip = 1)
  if (ncol(raw) < 2) {
    return(out)
  }

  nm <- trimws(names(raw))
  nm_low <- tolower(nm)
  gene_idx <- which(nm_low == "gene_id")[1]
  spec_idx <- which(grepl("testis", nm_low) & grepl("specificity", nm_low))[1]

  if (is.na(gene_idx) || is.na(spec_idx)) {
    return(out)
  }

  out <- data.frame(
    Geneid = as.character(raw[[gene_idx]]),
    testis_mated_4day_tissue_specificity_index = suppressWarnings(as.numeric(raw[[spec_idx]])),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$Geneid) & nzchar(out$Geneid), , drop = FALSE]
  out <- out[order(is.na(out$testis_mated_4day_tissue_specificity_index)), , drop = FALSE]
  out <- out[!duplicated(out$Geneid), , drop = FALSE]

  if (!is.null(target_genes)) {
    out <- out[out$Geneid %in% target_genes, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

read_scrna_db_for_genes <- function(scrna_db_path, target_genes, chunk_size = 100000L) {
  out <- data.frame(
    Geneid = character(0),
    Pub_ID = character(0),
    Cluster_Name = character(0),
    Mean_Expression = numeric(0),
    stringsAsFactors = FALSE
  )

  target_genes <- unique(target_genes[!is.na(target_genes) & nzchar(target_genes)])
  if (length(target_genes) == 0) {
    return(out)
  }
  if (is.na(scrna_db_path) || !nzchar(scrna_db_path) || !file.exists(scrna_db_path)) {
    return(out)
  }

  con <- if (is_gzipped_file(scrna_db_path)) gzfile(scrna_db_path, "rt") else file(scrna_db_path, "rt")
  on.exit(close(con), add = TRUE)

  header_line <- character(0)
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) {
      return(out)
    }
    if (startsWith(line, "##")) {
      next
    }
    header_line <- line
    break
  }

  header <- sub("^#", "", header_line)
  cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  required <- c("Pub_ID", "Cluster_Name", "Gene_ID", "Mean_Expression")
  idx <- match(required, cols)
  while (any(is.na(idx))) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) {
      return(out)
    }
    if (startsWith(line, "##")) {
      next
    }
    header <- sub("^#", "", line)
    cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
    idx <- match(required, cols)
  }

  chunks <- list()
  n_chunks <- 0L
  max_idx <- max(idx)
  repeat {
    chunk <- tryCatch(
      read.delim(
        con,
        header = FALSE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE,
        fill = TRUE,
        nrows = chunk_size
      ),
      error = function(e) {
        if (grepl("no lines available in input", conditionMessage(e), fixed = TRUE)) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (is.null(chunk) || nrow(chunk) == 0) {
      break
    }
    if (ncol(chunk) < max_idx) {
      next
    }

    gene_ids <- as.character(chunk[[idx[3]]])
    keep <- !is.na(gene_ids) & gene_ids %in% target_genes
    if (!any(keep)) {
      next
    }

    kept <- chunk[keep, idx, drop = FALSE]
    colnames(kept) <- required
    kept$Geneid <- as.character(kept$Gene_ID)
    kept$Mean_Expression <- suppressWarnings(as.numeric(kept$Mean_Expression))

    chunks[[length(chunks) + 1L]] <- kept[, c("Geneid", "Pub_ID", "Cluster_Name", "Mean_Expression"), drop = FALSE]
    n_chunks <- n_chunks + 1L
    if (n_chunks %% 10L == 0L) {
      gc(verbose = FALSE)
    }
  }

  if (length(chunks) == 0) {
    return(out)
  }
  out <- do.call(rbind, chunks)
  rownames(out) <- NULL
  out
}

annotate_pc_loading_with_expression <- function(load_df, output_prefix, testis_db_path,
                                                scrna_db_path, control_genes,
                                                testis_prefetched = NULL, scrna_prefetched = NULL) {
  if (is.null(load_df) || nrow(load_df) == 0) {
    return(invisible(NULL))
  }

  load_df <- load_df[!duplicated(load_df$Geneid), , drop = FALSE]
  load_df$control_rank <- match(load_df$Geneid, control_genes)

  missing_controls <- setdiff(control_genes, load_df$Geneid)
  if (length(missing_controls) > 0) {
    add_df <- data.frame(
      Geneid = missing_controls,
      PC1_loading = NA_real_,
      abs_loading = NA_real_,
      contribution_pct_PC1 = NA_real_,
      control_rank = match(missing_controls, control_genes),
      stringsAsFactors = FALSE
    )
    load_df <- rbind(load_df, add_df)
  }

  load_df$is_control_gene <- ifelse(is.na(load_df$control_rank), "no", "yes")
  load_df$sort_abs <- ifelse(is.na(load_df$abs_loading), -Inf, load_df$abs_loading)
  load_df <- load_df[order(is.na(load_df$control_rank), load_df$control_rank, -load_df$sort_abs), , drop = FALSE]
  load_df$sort_abs <- NULL

  testis_df <- if (is.null(testis_prefetched)) {
    read_testis_specificity_db(testis_db_path, target_genes = load_df$Geneid)
  } else {
    testis_prefetched[testis_prefetched$Geneid %in% load_df$Geneid, , drop = FALSE]
  }
  load_df$testis_mated_4day_tissue_specificity_index <- NA_real_
  if (nrow(testis_df) > 0) {
    m <- match(load_df$Geneid, testis_df$Geneid)
    load_df$testis_mated_4day_tissue_specificity_index <- testis_df$testis_mated_4day_tissue_specificity_index[m]
  }

  scrna_df <- if (is.null(scrna_prefetched)) {
    read_scrna_db_for_genes(scrna_db_path, target_genes = load_df$Geneid)
  } else {
    scrna_prefetched[scrna_prefetched$Geneid %in% load_df$Geneid, , drop = FALSE]
  }
  load_df$scrna_records <- 0L
  load_df$scrna_unique_pub_ids <- NA_character_
  load_df$scrna_unique_clusters <- NA_character_
  load_df$scrna_mean_expression_mean <- NA_real_
  load_df$scrna_mean_expression_max <- NA_real_

  if (nrow(scrna_df) > 0) {
    spl <- split(scrna_df, scrna_df$Geneid)
    scrna_summary <- do.call(
      rbind,
      lapply(names(spl), function(g) {
        d <- spl[[g]]
        data.frame(
          Geneid = g,
          scrna_records = nrow(d),
          scrna_unique_pub_ids = paste(sort(unique(as.character(d$Pub_ID))), collapse = ";"),
          scrna_unique_clusters = paste(sort(unique(as.character(d$Cluster_Name))), collapse = ";"),
          scrna_mean_expression_mean = suppressWarnings(mean(d$Mean_Expression, na.rm = TRUE)),
          scrna_mean_expression_max = suppressWarnings(max(d$Mean_Expression, na.rm = TRUE)),
          stringsAsFactors = FALSE
        )
      })
    )

    m <- match(load_df$Geneid, scrna_summary$Geneid)
    load_df$scrna_records <- ifelse(is.na(m), 0L, as.integer(scrna_summary$scrna_records[m]))
    load_df$scrna_unique_pub_ids <- scrna_summary$scrna_unique_pub_ids[m]
    load_df$scrna_unique_clusters <- scrna_summary$scrna_unique_clusters[m]
    load_df$scrna_mean_expression_mean <- scrna_summary$scrna_mean_expression_mean[m]
    load_df$scrna_mean_expression_max <- scrna_summary$scrna_mean_expression_max[m]
  }

  load_df$order_index <- seq_len(nrow(load_df))
  write_tsv(load_df[, setdiff(colnames(load_df), "order_index"), drop = FALSE], paste0(output_prefix, "_with_testis_scrna.tsv"))

  long_df <- merge(
    load_df[, c("Geneid", "order_index", "PC1_loading", "abs_loading", "contribution_pct_PC1",
                "is_control_gene", "testis_mated_4day_tissue_specificity_index"), drop = FALSE],
    scrna_df,
    by = "Geneid",
    all.x = TRUE,
    sort = FALSE
  )
  ord_expr <- ifelse(is.na(long_df$Mean_Expression), -Inf, long_df$Mean_Expression)
  long_df <- long_df[order(long_df$order_index, -ord_expr), , drop = FALSE]
  long_df$order_index <- NULL
  write_tsv(long_df, paste0(output_prefix, "_with_testis_scrna_long.tsv"))

  status_lines <- c(
    sprintf("testis_db_path=%s", testis_db_path),
    sprintf("scrna_db_path=%s", scrna_db_path),
    sprintf("n_genes_annotated=%d", nrow(load_df)),
    sprintf("n_scrna_rows_matched=%d", nrow(scrna_df))
  )
  writeLines(status_lines, con = paste0(output_prefix, "_annotation_status.txt"))

  invisible(load_df)
}

plot_heatmap <- function(vsd, genes, coldata, title, png_path, pdf_path) {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(vsd)]
  if (length(genes) < 2) {
    return(invisible(FALSE))
  }

  mat <- assay(vsd)[genes, , drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  ord <- order(coldata$timepoint)
  mat <- mat[, ord, drop = FALSE]
  cd <- coldata[colnames(mat), , drop = FALSE]

  tp_palette <- c("0H" = "#1b9e77", "2H" = "#d95f02", "4H" = "#7570b3", "6H" = "#e7298a")
  side_cols <- tp_palette[as.character(cd$timepoint)]

  draw <- function() {
    heatmap(
      mat,
      Rowv = NA,
      Colv = NA,
      scale = "none",
      col = colorRampPalette(c("navy", "white", "firebrick3"))(101),
      ColSideColors = side_cols,
      margins = c(9, 8),
      cexCol = 0.8,
      cexRow = 0.5,
      main = title,
      labRow = rownames(mat),
      labCol = colnames(mat)
    )
    legend("topright", legend = names(tp_palette), fill = tp_palette, border = NA, bty = "n", cex = 0.8)
  }

  png(filename = png_path, width = 1900, height = 1500, res = 220)
  draw()
  dev.off()

  pdf(file = pdf_path, width = 10, height = 8)
  draw()
  dev.off()

  invisible(TRUE)
}

run_analysis <- function(count_mat, coldata, outdir, apply_filter = TRUE, padj_thr = 0.05, lfc_thr = 1,
                         annotate_external = FALSE, testis_db_path = NA_character_,
                         scrna_db_path = NA_character_, control_genes = character(0)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ timepoint
  )

  n_before <- nrow(dds)
  if (apply_filter) {
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
  }
  n_after <- nrow(dds)

  if (n_after < 10) {
    stop("Too few genes left after filtering; cannot run DE.", call. = FALSE)
  }

  lrt_dds <- DESeq(dds, test = "LRT", reduced = ~1, quiet = TRUE)
  global_res <- results(lrt_dds)
  global_df <- res_to_df(global_res)
  write_tsv(global_df, file.path(outdir, "global_lrt_results.tsv"))

  wald_dds <- DESeq(dds, test = "Wald", quiet = TRUE)

  vsd <- vst(wald_dds, blind = FALSE)
  plot_pca(vsd, coldata = as.data.frame(colData(wald_dds)), prefix = outdir)

  global_sig <- subset(global_df, !is.na(padj) & padj < padj_thr)
  if (nrow(global_sig) > 0) {
    write_tsv(global_sig, file.path(outdir, "sig_global_lrt.tsv"))
  }

  pca_rare <- plot_pca_global_de(
    vsd,
    coldata = as.data.frame(colData(wald_dds)),
    global_genes = global_sig$Geneid,
    prefix = outdir
  )

  pca_time <- plot_pca_global_de_excluded_pc1(
    vsd,
    coldata = as.data.frame(colData(wald_dds)),
    global_genes = global_sig$Geneid,
    prefix = outdir,
    exclude_samples = c("HSF4_2H_3", "HSF4_4H_2", "HSF4_4H_3")
  )

  if (annotate_external) {
    prefetch_genes <- unique(c(
      control_genes,
      if (is.list(pca_rare) && !is.null(pca_rare$top50)) pca_rare$top50$Geneid else character(0),
      if (is.list(pca_time) && !is.null(pca_time$top50)) pca_time$top50$Geneid else character(0)
    ))
    testis_prefetched <- read_testis_specificity_db(testis_db_path, target_genes = prefetch_genes)
    scrna_prefetched <- read_scrna_db_for_genes(scrna_db_path, target_genes = prefetch_genes)

    annotate_pc_loading_with_expression(
      load_df = if (is.list(pca_rare)) pca_rare$top50 else NULL,
      output_prefix = file.path(outdir, "pc1_top50_loadings_contrib_global_de_all_samples"),
      testis_db_path = testis_db_path,
      scrna_db_path = scrna_db_path,
      control_genes = control_genes,
      testis_prefetched = testis_prefetched,
      scrna_prefetched = scrna_prefetched
    )
    annotate_pc_loading_with_expression(
      load_df = if (is.list(pca_time)) pca_time$top50 else NULL,
      output_prefix = file.path(outdir, "pc1_top50_loadings_contrib_global_de_excl_2H3_4H2_4H3"),
      testis_db_path = testis_db_path,
      scrna_db_path = scrna_db_path,
      control_genes = control_genes,
      testis_prefetched = testis_prefetched,
      scrna_prefetched = scrna_prefetched
    )
  }

  global_top <- head(global_sig$Geneid, 50)
  if (length(global_top) < 20) {
    global_top <- head(global_df$Geneid[!is.na(global_df$padj)], 50)
  }
  plot_heatmap(
    vsd,
    genes = global_top,
    coldata = as.data.frame(colData(wald_dds)),
    title = "Top global time-effect genes (LRT)",
    png_path = file.path(outdir, "heatmap_top_genes_global.png"),
    pdf_path = file.path(outdir, "heatmap_top_genes_global.pdf")
  )

  res_names <- resultsNames(wald_dds)
  contrasts <- c("2H", "4H", "6H")
  summary_rows <- list()

  for (tp in contrasts) {
    coef_name <- paste0("timepoint_", tp, "_vs_0H")
    if (!(coef_name %in% res_names)) {
      stop(sprintf("Expected coefficient not found: %s\nAvailable: %s", coef_name, paste(res_names, collapse = ", ")), call. = FALSE)
    }

    shr <- lfcShrink(wald_dds, coef = coef_name, type = "apeglm")
    cdf <- res_to_df(shr)

    contrast_name <- paste0(tp, "_vs_0H")
    write_tsv(cdf, file.path(outdir, paste0("contrast_", contrast_name, ".tsv")))

    sig <- subset(cdf, !is.na(padj) & padj < padj_thr & !is.na(log2FoldChange) & abs(log2FoldChange) >= lfc_thr)
    write_tsv(sig, file.path(outdir, paste0("sig_contrast_", contrast_name, ".tsv")))

    plot_volcano(
      cdf,
      title = paste0("Volcano ", contrast_name),
      png_path = file.path(outdir, paste0("volcano_", contrast_name, ".png")),
      pdf_path = file.path(outdir, paste0("volcano_", contrast_name, ".pdf")),
      padj_thr = padj_thr,
      lfc_thr = lfc_thr
    )

    top_genes <- head(sig$Geneid, 50)
    if (length(top_genes) < 20) {
      top_genes <- head(cdf$Geneid[!is.na(cdf$padj)], 50)
    }

    plot_heatmap(
      vsd,
      genes = top_genes,
      coldata = as.data.frame(colData(wald_dds)),
      title = paste0("Top genes ", contrast_name),
      png_path = file.path(outdir, paste0("heatmap_top_genes_", contrast_name, ".png")),
      pdf_path = file.path(outdir, paste0("heatmap_top_genes_", contrast_name, ".pdf"))
    )

    summary_rows[[contrast_name]] <- data.frame(
      analysis = basename(outdir),
      contrast = contrast_name,
      n_sig = nrow(sig),
      stringsAsFactors = FALSE
    )
  }

  metrics_df <- data.frame(
    metric = c("genes_before_filter", "genes_after_filter", "global_sig_lrt_padj_lt_0.05"),
    value = c(n_before, n_after, nrow(global_sig)),
    stringsAsFactors = FALSE
  )
  write_tsv(metrics_df, file.path(outdir, "run_metrics.tsv"))

  summary_df <- do.call(rbind, summary_rows)
  write_tsv(summary_df, file.path(outdir, "contrast_sig_summary.tsv"))

  list(
    n_before = n_before,
    n_after = n_after,
    global_sig = nrow(global_sig),
    contrast_sig = summary_df
  )
}

project_root_default <- find_project_root(getwd())

opts <- parse_args(commandArgs(trailingOnly = TRUE))
counts_path <- normalizePath(if (!is.null(opts$counts)) opts$counts else file.path(project_root_default, "02_work", "06_counts", "final", "gene_counts.matrix.tsv"), mustWork = FALSE)
samples_path <- normalizePath(if (!is.null(opts$samples)) opts$samples else file.path(project_root_default, "00_admin", "metadata", "samples.tsv"), mustWork = FALSE)
outdir <- normalizePath(if (!is.null(opts$outdir)) opts$outdir else file.path(project_root_default, "05_analysis", "de_timecourse"), mustWork = FALSE)
featurecounts_path <- normalizePath(if (!is.null(opts$featurecounts)) opts$featurecounts else file.path(project_root_default, "02_work", "06_counts", "final", "gene_counts.tsv"), mustWork = FALSE)
fc_summary_path <- normalizePath(if (!is.null(opts$fcsummary)) opts$fcsummary else file.path(project_root_default, "02_work", "06_counts", "final", "gene_counts.tsv.summary"), mustWork = FALSE)
testis_db_path <- normalizePath(if (!is.null(opts$testisdb)) opts$testisdb else file.path(project_root_default, "05_analysis", "12864_2018_5085_MOESM2_ESM.xlsx"), mustWork = FALSE)
scrna_db_path <- normalizePath(if (!is.null(opts$scrnadb)) opts$scrnadb else file.path(project_root_default, "05_analysis", "scRNA-Seq_gene_expression_fb_2025_05.tsv.gz.tsv"), mustWork = FALSE)

if (!file.exists(counts_path)) {
  stop(sprintf("Counts file not found: %s", counts_path), call. = FALSE)
}
if (!file.exists(samples_path)) {
  stop(sprintf("Samples metadata file not found: %s", samples_path), call. = FALSE)
}
if (!file.exists(featurecounts_path)) {
  stop(sprintf("featureCounts table not found: %s", featurecounts_path), call. = FALSE)
}
if (!file.exists(fc_summary_path)) {
  stop(sprintf("featureCounts summary not found: %s", fc_summary_path), call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "sensitivity_no_filter"), recursive = TRUE, showWarnings = FALSE)

counts_df <- read.delim(counts_path, check.names = FALSE)
if (ncol(counts_df) < 3) {
  stop("Counts matrix must contain Geneid + sample columns.", call. = FALSE)
}
if (!identical(colnames(counts_df)[1], "Geneid")) {
  stop("First column of counts matrix must be named 'Geneid'.", call. = FALSE)
}

count_mat <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(count_mat) <- counts_df$Geneid

if (any(is.na(count_mat))) {
  stop("Counts matrix contains NA values.", call. = FALSE)
}

storage.mode(count_mat) <- "numeric"
if (any(count_mat < 0)) {
  stop("Counts matrix contains negative values.", call. = FALSE)
}
if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop("Counts matrix must contain integer-like values.", call. = FALSE)
}
storage.mode(count_mat) <- "integer"

sample_names <- colnames(count_mat)

sample_df <- data.frame(sample = sample_names, stringsAsFactors = FALSE)
sample_df$timepoint <- sub("^HSF4_([0-9]+H)_.*$", "\\1", sample_df$sample)
if (any(sample_df$timepoint == sample_df$sample)) {
  stop("Failed to parse timepoint from one or more sample names.", call. = FALSE)
}

samples_meta <- read.delim(samples_path, check.names = FALSE)
if (!all(c("sample_id", "R1") %in% colnames(samples_meta))) {
  stop("samples.tsv must contain columns 'sample_id' and 'R1'.", call. = FALSE)
}

samples_meta$lane <- ifelse(grepl("_L7_", samples_meta$R1), "L7", ifelse(grepl("_L8_", samples_meta$R1), "L8", "NA"))
lane_map <- setNames(samples_meta$lane, samples_meta$sample_id)
sample_df$lane <- lane_map[sample_df$sample]
sample_df$lane[is.na(sample_df$lane)] <- "NA"

sample_df$timepoint <- factor(sample_df$timepoint, levels = c("0H", "2H", "4H", "6H"))
if (any(is.na(sample_df$timepoint))) {
  stop("Detected samples with timepoints outside expected levels: 0H,2H,4H,6H", call. = FALSE)
}
sample_df$lane <- factor(sample_df$lane, levels = c("L7", "L8", "NA"))

rownames(sample_df) <- sample_df$sample
sample_df <- sample_df[colnames(count_mat), , drop = FALSE]

if (!all(rownames(sample_df) == colnames(count_mat))) {
  stop("Sample metadata is not aligned with count matrix columns.", call. = FALSE)
}

if (ncol(count_mat) != 12) {
  warning(sprintf("Expected 12 samples; found %d", ncol(count_mat)))
}

time_counts <- table(sample_df$timepoint)
write_tsv(
  data.frame(timepoint = names(time_counts), n_samples = as.integer(time_counts), stringsAsFactors = FALSE),
  file.path(outdir, "sample_distribution.tsv")
)
write_tsv(
  data.frame(sample = rownames(sample_df), timepoint = sample_df$timepoint, lane = sample_df$lane, stringsAsFactors = FALSE),
  file.path(outdir, "coldata.tsv")
)

gtf_path <- infer_gtf_path(
  project_root = project_root_default,
  featurecounts_path = featurecounts_path,
  gtf_arg = opts$gtf
)

compute_sample_qc_metrics(
  featurecounts_path = featurecounts_path,
  fc_summary_path = fc_summary_path,
  coldata = sample_df,
  outdir = outdir,
  gtf_path = gtf_path
)

main_summary <- run_analysis(
  count_mat = count_mat,
  coldata = sample_df,
  outdir = outdir,
  apply_filter = TRUE,
  padj_thr = 0.05,
  lfc_thr = 1,
  annotate_external = TRUE,
  testis_db_path = testis_db_path,
  scrna_db_path = scrna_db_path,
  control_genes = control_genes_default
)

nofilter_dir <- file.path(outdir, "sensitivity_no_filter")
nofilter_summary <- run_analysis(
  count_mat = count_mat,
  coldata = sample_df,
  outdir = nofilter_dir,
  apply_filter = FALSE,
  padj_thr = 0.05,
  lfc_thr = 1,
  annotate_external = FALSE,
  testis_db_path = testis_db_path,
  scrna_db_path = scrna_db_path,
  control_genes = control_genes_default
)

sensitivity_summary <- data.frame(
  metric = c(
    "genes_before_filter_main",
    "genes_after_filter_main",
    "global_sig_main",
    "genes_before_filter_nofilter",
    "genes_after_filter_nofilter",
    "global_sig_nofilter"
  ),
  value = c(
    main_summary$n_before,
    main_summary$n_after,
    main_summary$global_sig,
    nofilter_summary$n_before,
    nofilter_summary$n_after,
    nofilter_summary$global_sig
  ),
  stringsAsFactors = FALSE
)
write_tsv(sensitivity_summary, file.path(outdir, "sensitivity_summary.tsv"))

report_path <- file.path(outdir, "analysis_report.txt")
con <- file(report_path, open = "wt")
writeLines(c(
  "Differential Expression Analysis Report",
  "====================================",
  "",
  sprintf("Counts file: %s", counts_path),
  sprintf("Samples file: %s", samples_path),
  sprintf("featureCounts table: %s", featurecounts_path),
  sprintf("featureCounts summary: %s", fc_summary_path),
  sprintf("GTF used for rRNA annotation: %s", ifelse(is.na(gtf_path), "NA", gtf_path)),
  sprintf("Testis specificity DB: %s", testis_db_path),
  sprintf("scRNA expression DB: %s", scrna_db_path),
  sprintf("Output directory: %s", outdir),
  "",
  "Design and tests:",
  "- Primary model: ~ timepoint (categorical; reference 0H)",
  "- Global effect: DESeq2 LRT (reduced = ~1)",
  "- Pairwise contrasts: 2H vs 0H, 4H vs 0H, 6H vs 0H (Wald + lfcShrink apeglm)",
  "- Significance threshold: padj < 0.05 and |log2FoldChange| >= 1",
  "",
  "Sample checks:",
  sprintf("- n samples: %d", ncol(count_mat)),
  sprintf("- timepoint counts: %s", paste(names(time_counts), as.integer(time_counts), sep = "=", collapse = ", ")),
  "",
  "Confounding warning:",
  "- lane and timepoint are partially confounded in this dataset (0H/2H in L8; 4H/6H in L7).",
  "- Cross-lane contrasts (e.g., 4H vs 0H, 6H vs 0H) may reflect both biology and lane effects.",
  "",
  "Main (filtered) summary:",
  sprintf("- genes before filter: %d", main_summary$n_before),
  sprintf("- genes after filter: %d", main_summary$n_after),
  sprintf("- global LRT significant genes (padj<0.05): %d", main_summary$global_sig),
  "",
  "No-filter sensitivity summary:",
  sprintf("- genes before filter: %d", nofilter_summary$n_before),
  sprintf("- genes after filter: %d", nofilter_summary$n_after),
  sprintf("- global LRT significant genes (padj<0.05): %d", nofilter_summary$global_sig),
  "",
  "Primary outputs:",
  "- sample_qc_metrics.tsv (assigned%, rRNA%, mito%, genes_detected per sample)",
  "- global_lrt_results.tsv",
  "- contrast_2H_vs_0H.tsv, contrast_4H_vs_0H.tsv, contrast_6H_vs_0H.tsv",
  "- sig_contrast_*.tsv and sig_global_lrt.tsv",
  "- pca_timepoint.(png/pdf), pca_lane.(png/pdf)",
  "- pca_global_de_timepoint.(png/pdf)",
  "- pca_global_de_timepoint_excl_2H3_4H2_4H3.(png/pdf)",
  "- pc1_all_loadings_contrib_global_de_all_samples.tsv",
  "- pc1_top50_loadings_contrib_global_de_all_samples.tsv",
  "- go_quick_pc1_top50_global_de_all_samples.(tsv/png/pdf or status.txt)",
  "- pc1_top50_loadings_contrib_global_de_all_samples_with_testis_scrna.tsv",
  "- pc1_top50_loadings_contrib_global_de_all_samples_with_testis_scrna_long.tsv",
  "- pc1_top50_loadings_contrib_global_de_excl_2H3_4H2_4H3.tsv",
  "- go_quick_pc1_top50_global_de_excl_2H3_4H2_4H3.(tsv/png/pdf or status.txt)",
  "- pc1_top50_loadings_contrib_global_de_excl_2H3_4H2_4H3_with_testis_scrna.tsv",
  "- pc1_top50_loadings_contrib_global_de_excl_2H3_4H2_4H3_with_testis_scrna_long.tsv",
  "- volcano_*.png/pdf",
  "- heatmap_top_genes_*.png/pdf"
), con)
close(con)

message("[INFO] DE analysis completed successfully.")
message("[INFO] Output directory: ", outdir)
