suppressPackageStartupMessages({
  library(vegan)
})

pairwise_permanova_base <- function(dist_obj, group, permutations = 999) {
  group <- droplevels(as.factor(group))
  levs <- levels(group)
  out <- list()

  if (length(levs) < 2) {
    return(data.frame())
  }

  idx <- 1
  for (i in 1:(length(levs) - 1)) {
    for (j in (i + 1):length(levs)) {
      keep <- group %in% c(levs[i], levs[j])
      g_sub <- droplevels(group[keep])

      d_mat <- as.matrix(dist_obj)
      d_sub <- as.dist(d_mat[keep, keep])

      ad <- adonis2(d_sub ~ g_sub, permutations = permutations)

      out[[idx]] <- data.frame(
        group1 = levs[i],
        group2 = levs[j],
        F = unname(ad$F[1]),
        R2 = unname(ad$R2[1]),
        P = unname(ad$`Pr(>F)`[1]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  res <- do.call(rbind, out)
  res$P_adj_BH <- p.adjust(res$P, method = "BH")
  res
}

run_one_analysis <- function(mat_file,
                             meta_file,
                             data_type = c("ARG", "VIR", "REP"),
                             group_var = c("Country", "Source"),
                             dist_method = "bray",
                             binary_jaccard = FALSE,
                             permutations = 999) {

  data_type <- match.arg(data_type)
  group_var <- match.arg(group_var)

  cat("\n==============================\n")
  cat("Running:", data_type, "~", group_var, "\n")
  cat("==============================\n")

  # 读矩阵
  mat_df <- read.delim(mat_file, check.names = FALSE, stringsAsFactors = FALSE)

  if (!("SampleID" %in% colnames(mat_df))) {
    colnames(mat_df)[1] <- "SampleID"
  }

  mat_df$SampleID <- as.character(mat_df$SampleID)

  mat <- as.matrix(mat_df[, setdiff(colnames(mat_df), "SampleID"), drop = FALSE])
  rownames(mat) <- mat_df$SampleID
  storage.mode(mat) <- "numeric"

  # 读元数据
  meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  meta$sample_id <- as.character(meta$sample_id)

  # 统一分组变量
  if (group_var == "Country") {
    meta$Group <- trimws(as.character(meta$region))
  } else {
    meta$Group <- trimws(as.character(meta$source))
    meta$Group[tolower(meta$Group) == "clinical"] <- "human"
  }

  # 去 unknown / NA / 空值
  bad <- is.na(meta$Group) |
         meta$Group == "" |
         meta$Group == "NA" |
         tolower(meta$Group) == "unknown"
  meta <- meta[!bad, , drop = FALSE]

  # 对齐样本
  common_ids <- intersect(rownames(mat), meta$sample_id)
  mat <- mat[common_ids, , drop = FALSE]
  meta <- meta[match(common_ids, meta$sample_id), , drop = FALSE]

  if (!identical(rownames(mat), meta$sample_id)) {
    stop("Sample order mismatch after matching.")
  }

  # 去掉全0行
  keep_nonzero <- rowSums(mat, na.rm = TRUE) > 0
  mat <- mat[keep_nonzero, , drop = FALSE]
  meta <- meta[keep_nonzero, , drop = FALSE]

  # 只保留 top3 组，和你现在作图保持一致
  grp_tab <- sort(table(meta$Group), decreasing = TRUE)
  top_groups <- names(grp_tab)[1:min(3, length(grp_tab))]
  keep_top <- meta$Group %in% top_groups
  mat <- mat[keep_top, , drop = FALSE]
  meta <- meta[keep_top, , drop = FALSE]
  grp <- droplevels(factor(meta$Group, levels = top_groups))

  if (nrow(mat) < 3) {
    stop(paste("Too few samples for", data_type, "~", group_var))
  }
  if (nlevels(grp) < 2) {
    stop(paste("Too few groups for", data_type, "~", group_var))
  }

  cat("Samples kept:", nrow(mat), "\n")
  cat("Features:", ncol(mat), "\n")
  cat("Groups:\n")
  print(table(grp))

  # 距离
  if (dist_method == "jaccard") {
    d <- vegdist(mat, method = "jaccard", binary = binary_jaccard)
  } else {
    d <- vegdist(mat, method = dist_method)
  }

  # PERMANOVA
  ad <- adonis2(d ~ grp, permutations = permutations)

  # betadisper
  bd <- betadisper(d, grp)
  bd_perm <- permutest(bd, permutations = permutations)

  # pairwise PERMANOVA
  pw <- pairwise_permanova_base(d, grp, permutations = permutations)

  # 输出文件前缀
  prefix <- paste0(data_type, "_", group_var)

  # 全局结果表
  global_res <- data.frame(
    dataset = data_type,
    grouping = group_var,
    distance = ifelse(dist_method == "jaccard", "Jaccard", "Bray-Curtis"),
    n_samples = nrow(mat),
    n_features = ncol(mat),
    n_groups = nlevels(grp),
    groups = paste(levels(grp), collapse = ";"),
    permanova_F = unname(ad$F[1]),
    permanova_R2 = unname(ad$R2[1]),
    permanova_P = unname(ad$`Pr(>F)`[1]),
    betadisper_F = unname(bd_perm$tab[1, "F"]),
    betadisper_P = unname(bd_perm$tab[1, "Pr(>F)"]),
    stringsAsFactors = FALSE
  )

  write.table(global_res,
              file = paste0(prefix, "_global_stats.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (nrow(pw) > 0) {
    write.table(pw,
                file = paste0(prefix, "_pairwise_PERMANOVA.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # 更详细的文本输出
  sink(paste0(prefix, "_details.txt"))
  cat("Dataset:", data_type, "\n")
  cat("Grouping:", group_var, "\n")
  cat("Distance:", ifelse(dist_method == "jaccard", "Jaccard", "Bray-Curtis"), "\n\n")

  cat("Sample counts per group:\n")
  print(table(grp))

  cat("\nPERMANOVA (adonis2):\n")
  print(ad)

  cat("\nHomogeneity of dispersion (betadisper ANOVA):\n")
  print(anova(bd))

  cat("\nHomogeneity of dispersion (permutest):\n")
  print(bd_perm)

  if (nrow(pw) > 0) {
    cat("\nPairwise PERMANOVA:\n")
    print(pw)
  }
  sink()

  return(global_res)
}

all_global <- rbind(
  run_one_analysis("ARG_matrix_01_01_for_PCoA.tsv",      "final_dataset_for_test.tsv", "ARG", "Country"),
  run_one_analysis("ARG_matrix_01_01_for_PCoA.tsv",      "final_dataset_for_test.tsv", "ARG", "Source"),
  run_one_analysis("VIR_matrix_full_01_01_for_PCoA.tsv", "final_dataset_for_test.tsv", "VIR", "Country"),
  run_one_analysis("VIR_matrix_full_01_01_for_PCoA.tsv", "final_dataset_for_test.tsv", "VIR", "Source"),
  run_one_analysis("REP_matrix_01_01_for_PCoA.tsv",      "final_dataset_for_test.tsv", "REP", "Country"),
  run_one_analysis("REP_matrix_01_01_for_PCoA.tsv",      "final_dataset_for_test.tsv", "REP", "Source")
)

write.table(all_global,
            file = "PCoA_all_global_stats_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nAll analyses finished.\n")
cat("Summary file: PCoA_all_global_stats_summary.tsv\n")
