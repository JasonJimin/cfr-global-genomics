suppressPackageStartupMessages({
  library(vegan)
})

# =========================
# 参数
# =========================
meta_file    <- "final_dataset_for_test.tsv"
species_file <- "cfrV.csv"

# 与你当前PCoA保持一致：只分析每个分组变量中样本数最多的前3组
use_top3_groups <- TRUE

# 为避免 species 水平太碎导致模型过拟合，
# 将样本数 < 3 的 species 合并为 "Other"
min_species_n <- 3

distance_method <- "bray"
nperm <- 999

# =========================
# 工具函数
# =========================
read_matrix_file <- function(mat_file) {
  x <- read.delim(mat_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"SampleID" %in% colnames(x)) {
    colnames(x)[1] <- "SampleID"
  }
  x$SampleID <- as.character(x$SampleID)

  mat <- as.matrix(x[, setdiff(colnames(x), "SampleID"), drop = FALSE])
  rownames(mat) <- x$SampleID
  storage.mode(mat) <- "numeric"
  mat
}

read_species_file <- function(species_file) {
  sp <- read.csv(species_file, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(sp) <- trimws(colnames(sp))

  need <- c("strain_id", "Scientific_name")
  miss <- setdiff(need, colnames(sp))
  if (length(miss) > 0) {
    stop("species file缺少列: ", paste(miss, collapse = ", "))
  }

  sp <- sp[, c("strain_id", "Scientific_name"), drop = FALSE]
  sp$strain_id <- trimws(as.character(sp$strain_id))
  sp$Scientific_name <- trimws(as.character(sp$Scientific_name))
  sp <- sp[sp$strain_id != "" & sp$Scientific_name != "", , drop = FALSE]
  sp <- sp[!duplicated(sp$strain_id), , drop = FALSE]
  sp
}

prepare_metadata <- function(meta_file, species_file, group_var) {
  meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  meta$sample_id <- trimws(as.character(meta$sample_id))

  sp <- read_species_file(species_file)

  meta <- merge(
    meta,
    sp,
    by.x = "sample_id",
    by.y = "strain_id",
    all.x = TRUE,
    sort = FALSE
  )

  if (group_var == "Country") {
    meta$Group <- trimws(as.character(meta$region))
  } else if (group_var == "Source") {
    meta$Group <- trimws(as.character(meta$source))
    meta$Group[tolower(meta$Group) == "clinical"] <- "human"
  } else {
    stop("group_var must be Country or Source")
  }

  # 去 unknown / NA / 空值
  bad_group <- is.na(meta$Group) |
               meta$Group == "" |
               meta$Group == "NA" |
               tolower(meta$Group) == "unknown"

  bad_species <- is.na(meta$Scientific_name) |
                 meta$Scientific_name == "" |
                 meta$Scientific_name == "NA" |
                 tolower(meta$Scientific_name) == "unknown"

  meta <- meta[!(bad_group | bad_species), , drop = FALSE]

  meta
}

collapse_rare_species <- function(species_vec, min_n = 3) {
  species_vec <- as.character(species_vec)
  tb <- table(species_vec)
  out <- ifelse(tb[species_vec] < min_n, "Other", species_vec)
  factor(out)
}

get_adonis_row <- function(ad, term_name) {
  rn <- rownames(ad)
  if (!(term_name %in% rn)) return(NULL)
  ad[term_name, , drop = FALSE]
}

pairwise_permanova_adjusted <- function(dist_obj, grp, sp, permutations = 999) {
  grp <- droplevels(as.factor(grp))
  sp  <- droplevels(as.factor(sp))
  levs <- levels(grp)

  if (length(levs) < 2) {
    return(data.frame())
  }

  out <- list()
  k <- 1

  for (i in 1:(length(levs)-1)) {
    for (j in (i+1):length(levs)) {

      keep <- grp %in% c(levs[i], levs[j])
      g_sub <- droplevels(grp[keep])
      s_sub <- droplevels(sp[keep])

      dmat <- as.matrix(dist_obj)
      dsub <- as.dist(dmat[keep, keep])

      # 如果子集里species只有一个水平，就退化成普通 pairwise adonis2
      if (nlevels(s_sub) >= 2) {
        ad_sub <- adonis2(dsub ~ g_sub + s_sub, by = "margin", permutations = permutations)
        row_sub <- get_adonis_row(ad_sub, "g_sub")
        model_type <- "group + species (margin)"
      } else {
        ad_sub <- adonis2(dsub ~ g_sub, permutations = permutations)
        row_sub <- get_adonis_row(ad_sub, "g_sub")
        model_type <- "group only (species invariant)"
      }

      out[[k]] <- data.frame(
        group1 = levs[i],
        group2 = levs[j],
        model = model_type,
        F = unname(row_sub$F[1]),
        R2 = unname(row_sub$R2[1]),
        P = unname(row_sub$`Pr(>F)`[1]),
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  res <- do.call(rbind, out)
  res$P_adj_BH <- p.adjust(res$P, method = "BH")
  res
}

run_one <- function(mat_file, data_type, group_var) {
  cat("\n=============================\n")
  cat("Running:", data_type, "~", group_var, "\n")
  cat("=============================\n")

  mat  <- read_matrix_file(mat_file)
  meta <- prepare_metadata(meta_file, species_file, group_var)

  common_ids <- intersect(rownames(mat), meta$sample_id)
  mat  <- mat[common_ids, , drop = FALSE]
  meta <- meta[match(common_ids, meta$sample_id), , drop = FALSE]

  if (!identical(rownames(mat), meta$sample_id)) {
    stop("样本顺序未对齐: ", data_type, " ~ ", group_var)
  }

  # 去全0行
  keep_nonzero <- rowSums(mat, na.rm = TRUE) > 0
  mat  <- mat[keep_nonzero, , drop = FALSE]
  meta <- meta[keep_nonzero, , drop = FALSE]

  # 保持与现有PCoA图一致，只取top3组
  if (use_top3_groups) {
    tb <- sort(table(meta$Group), decreasing = TRUE)
    top_groups <- names(tb)[1:min(3, length(tb))]
    keep_top <- meta$Group %in% top_groups
    mat  <- mat[keep_top, , drop = FALSE]
    meta <- meta[keep_top, , drop = FALSE]
  }

  grp <- droplevels(factor(meta$Group))
  sp_raw <- factor(meta$Scientific_name)
  sp <- collapse_rare_species(meta$Scientific_name, min_species_n)

  if (nrow(mat) < 3) stop("样本数太少: ", data_type, " ~ ", group_var)
  if (nlevels(grp) < 2) stop("组别不足: ", data_type, " ~ ", group_var)

  cat("Samples kept:", nrow(mat), "\n")
  cat("Features:", ncol(mat), "\n")
  cat("Group counts:\n")
  print(table(grp))
  cat("Raw species levels:", nlevels(sp_raw), "\n")
  cat("Collapsed species levels:", nlevels(sp), "\n")

  d <- vegdist(mat, method = distance_method)

  # 1) 原始 PERMANOVA
  ad_raw <- adonis2(d ~ grp, permutations = nperm)
  row_raw <- get_adonis_row(ad_raw, "grp")

  # 2) 控制 species 后的 PERMANOVA
  # 用 by="margin" 提取 group 的独立效应
  if (nlevels(sp) >= 2) {
    ad_adj <- adonis2(d ~ grp + sp, by = "margin", permutations = nperm)
    row_grp_adj <- get_adonis_row(ad_adj, "grp")
    row_sp_adj  <- get_adonis_row(ad_adj, "sp")

    # 再给一个顺序模型，直观看“species先进入后，group还剩多少解释力”
    ad_seq <- adonis2(d ~ sp + grp, by = "terms", permutations = nperm)
    row_grp_seq <- get_adonis_row(ad_seq, "grp")
  } else {
    ad_adj <- NULL
    row_grp_adj <- NULL
    row_sp_adj <- NULL
    ad_seq <- NULL
    row_grp_seq <- NULL
  }

  # 3) dispersion
  bd <- betadisper(d, grp)
  bd_perm <- permutest(bd, permutations = nperm)

  # 4) pairwise，控制species
  pw <- pairwise_permanova_adjusted(d, grp, sp, permutations = nperm)

  prefix <- paste0(data_type, "_", group_var, "_species_adjusted")

  global_res <- data.frame(
    dataset = data_type,
    grouping = group_var,
    distance = "Bray-Curtis",
    n_samples = nrow(mat),
    n_features = ncol(mat),
    groups = paste(levels(grp), collapse = ";"),
    n_group_levels = nlevels(grp),
    n_species_raw = nlevels(sp_raw),
    n_species_collapsed = nlevels(sp),

    raw_PERMANOVA_F = unname(row_raw$F[1]),
    raw_PERMANOVA_R2 = unname(row_raw$R2[1]),
    raw_PERMANOVA_P = unname(row_raw$`Pr(>F)`[1]),

    adj_group_F = if (!is.null(row_grp_adj)) unname(row_grp_adj$F[1]) else NA,
    adj_group_R2 = if (!is.null(row_grp_adj)) unname(row_grp_adj$R2[1]) else NA,
    adj_group_P = if (!is.null(row_grp_adj)) unname(row_grp_adj$`Pr(>F)`[1]) else NA,

    adj_species_F = if (!is.null(row_sp_adj)) unname(row_sp_adj$F[1]) else NA,
    adj_species_R2 = if (!is.null(row_sp_adj)) unname(row_sp_adj$R2[1]) else NA,
    adj_species_P = if (!is.null(row_sp_adj)) unname(row_sp_adj$`Pr(>F)`[1]) else NA,

    seq_group_after_species_F = if (!is.null(row_grp_seq)) unname(row_grp_seq$F[1]) else NA,
    seq_group_after_species_R2 = if (!is.null(row_grp_seq)) unname(row_grp_seq$R2[1]) else NA,
    seq_group_after_species_P = if (!is.null(row_grp_seq)) unname(row_grp_seq$`Pr(>F)`[1]) else NA,

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

  # species频次
  sp_count <- sort(table(sp_raw), decreasing = TRUE)
  write.table(
    data.frame(Species = names(sp_count), n = as.integer(sp_count)),
    file = paste0(prefix, "_species_counts.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # 详细结果
  sink(paste0(prefix, "_details.txt"))
  cat("Dataset:", data_type, "\n")
  cat("Grouping:", group_var, "\n")
  cat("Distance:", distance_method, "\n")
  cat("use_top3_groups:", use_top3_groups, "\n")
  cat("min_species_n:", min_species_n, "\n\n")

  cat("Group counts:\n")
  print(table(grp))

  cat("\nRaw species counts (top 20):\n")
  print(head(sort(table(sp_raw), decreasing = TRUE), 20))

  cat("\nCollapsed species counts:\n")
  print(sort(table(sp), decreasing = TRUE))

  cat("\nRaw PERMANOVA (group only):\n")
  print(ad_raw)

  cat("\nAdjusted PERMANOVA (group + species, by='margin'):\n")
  if (!is.null(ad_adj)) {
    print(ad_adj)
  } else {
    cat("Species只有一个水平，未运行adjusted model。\n")
  }

  cat("\nSequential PERMANOVA (species first, then group; by='terms'):\n")
  if (!is.null(ad_seq)) {
    print(ad_seq)
  } else {
    cat("Species只有一个水平，未运行sequential model。\n")
  }

  cat("\nBetadisper ANOVA:\n")
  print(anova(bd))

  cat("\nBetadisper permutation test:\n")
  print(bd_perm)

  cat("\nPairwise species-adjusted PERMANOVA:\n")
  if (nrow(pw) > 0) print(pw)

  sink()

  global_res
}

all_res <- rbind(
  run_one("ARG_matrix_01_01_for_PCoA.tsv",      "ARG", "Country"),
  run_one("ARG_matrix_01_01_for_PCoA.tsv",      "ARG", "Source"),
  run_one("VIR_matrix_full_01_01_for_PCoA.tsv", "VIR", "Country"),
  run_one("VIR_matrix_full_01_01_for_PCoA.tsv", "VIR", "Source"),
  run_one("REP_matrix_01_01_for_PCoA.tsv",      "REP", "Country"),
  run_one("REP_matrix_01_01_for_PCoA.tsv",      "REP", "Source")
)

write.table(all_res,
            file = "PCoA_species_adjusted_all_global_stats_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nAll species-adjusted analyses finished.\n")
cat("Main summary: PCoA_species_adjusted_all_global_stats_summary.tsv\n")
