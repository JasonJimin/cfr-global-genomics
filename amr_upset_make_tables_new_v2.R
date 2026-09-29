# ============================================================
# amr_upset_make_tables_new_v2.R
# 修复点：自动跳过 ResFinder.tsv 开头的 “## ResFinder” 等段落标记行
# ============================================================

## 让 R 优先使用 conda 环境中的 R 包（如在 conda 中运行）
.libPaths(unique(c(file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library"), .libPaths())))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(readr); library(purrr); library(tibble)
})

basedir <- getwd()

# ---- 兼容命名：优先 meta.tsv；否则自动找你生成的 meta_as_metadata_with_cfr.tsv ----
meta_candidates <- c("meta.tsv", "meta_as_metadata_with_cfr.tsv", "metadata.tsv", "meta_as_metadata.tsv")
meta_file <- meta_candidates[file.exists(file.path(basedir, meta_candidates))][1]
if (is.na(meta_file)) stop("找不到元信息文件：请提供 meta.tsv（或 meta_as_metadata_with_cfr.tsv / metadata.tsv）")
meta_file <- file.path(basedir, meta_file)

# ---- RGI 预测表型文件（两列：sample_id<TAB>Predicted_phenotypes；通常无表头） ----
rgi_pheno_file <- file.path(basedir, "rgi_predicted_phenotypes.tsv")
if (!file.exists(rgi_pheno_file)) stop("找不到 rgi_predicted_phenotypes.tsv")

# ---- ResFinder 文件：优先 ResFinder.tsv；兼容旧名 all_resfinder.tsv ----
res_candidates <- c("ResFinder.tsv", "all_resfinder.tsv")
resf_file <- res_candidates[file.exists(file.path(basedir, res_candidates))][1]
if (is.na(resf_file)) {
  warning("未找到 ResFinder.tsv / all_resfinder.tsv，将只使用 RGI 预测表型。")
  resf_file <- NA_character_
} else {
  resf_file <- file.path(basedir, resf_file)
}

outdir <- file.path(basedir, "amr_upset_figs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Using files:")
message("  meta: ", meta_file)
message("  rgi : ", rgi_pheno_file)
message("  res : ", ifelse(is.na(resf_file), "(none)", resf_file))
message("Output dir: ", outdir)

## ===== 归一化到这些“药物类别”（按需增删） =====
canon_classes <- c(
  "carbapenem","cephalosporin","penicillin","beta-lactamase inhibitor",
  "fluoroquinolone","aminoglycoside","tetracycline","polymyxin",
  "macrolide","sulfonamide","trimethoprim","phenicol","glycopeptide",
  "fosfomycin","rifamycin","nitroimidazole",
  "lincosamide","streptogramin","oxazolidinone","pleuromutilin",
  "nitrofuran","aminocoumarin","nucleoside","peptide",
  "glycylcycline","disinfecting agents and antiseptics","monobactam"
)

## ===== 工具函数 =====
split_multi <- function(x){
  if (is.na(x) || is.null(x)) return(character())
  unlist(str_split(x, "\\s*[,;|/]\\s*")) %>% str_squish() %>% discard(~.x=="")
}

canonize <- function(v){
  lv <- tolower(v)
  lv <- str_replace_all(lv, c(
    "beta.?-?lactamase inhibitor.*"="beta-lactamase inhibitor",
    ".*polymyxin.*"="polymyxin",
    ".*carbapenem.*"="carbapenem",
    ".*cephalosporin.*"="cephalosporin",
    ".*monobactam.*"="monobactam",
    ".*glycylcycline.*"="glycylcycline",
    ".*quinolone.*|.*fluoroquinolone.*"="fluoroquinolone",
    ".*aminoglycoside.*"="aminoglycoside",
    ".*tetracycline.*"="tetracycline",
    ".*macrolide.*"="macrolide",
    ".*sulfonamide.*"="sulfonamide",
    ".*trimethoprim.*|.*diaminopyrimidine.*"="trimethoprim",
    ".*phenicol.*"="phenicol",
    ".*glycopeptide.*"="glycopeptide",
    ".*penicillin.*"="penicillin",
    ".*fosfomycin.*|.*phosphonic acid.*"="fosfomycin",
    ".*rifampin.*|.*rifamycin.*"="rifamycin",
    ".*metronidazole.*|.*nitroimidazole.*"="nitroimidazole",
    ".*nitrofuran.*"="nitrofuran",
    ".*lincosamide.*"="lincosamide",
    ".*streptogramin.*"="streptogramin",
    ".*oxazolidinone.*"="oxazolidinone",
    ".*pleuromutilin.*"="pleuromutilin",
    ".*aminocoumarin.*"="aminocoumarin",
    ".*nucleoside.*"="nucleoside",
    ".*peptide antibiotic.*"="peptide",
    ".*disinfecting agents and antiseptics.*"="disinfecting agents and antiseptics"
  ))
  lv[lv %in% canon_classes]
}

# 自动定位“表头行”并 skip 到那里（用于 ResFinder 这种开头有 ## 的情况）
read_tsv_after_header <- function(path, header_regex){
  lines <- readLines(path, warn = FALSE)
  idx <- which(grepl(header_regex, lines, ignore.case = TRUE))[1]
  if (is.na(idx)) {
    stop("无法在文件中找到表头行：", basename(path), "\n请检查是否是 TSV，且包含表头。")
  }
  # 从表头行开始读
  df <- readr::read_delim(
    file = path,
    delim = "\t",
    skip = idx - 1,
    col_names = TRUE,
    show_col_types = FALSE,
    progress = FALSE,
    quote = ""
  )
  # 如果第一列名带 #（例如 #FILE），去掉 #
  names(df)[1] <- sub("^#", "", names(df)[1])
  df
}

## ===== 读 meta =====
meta <- readr::read_delim(
  meta_file, delim="\t", col_names=TRUE,
  show_col_types=FALSE, progress=FALSE, quote=""
)
stopifnot("sample_id" %in% names(meta), "Year" %in% names(meta))
meta <- meta %>%
  mutate(Year = suppressWarnings(as.integer(as.character(.data$Year))))

## ===== 读 RGI 预测表型 =====
first_line <- readLines(rgi_pheno_file, n=1, warn=FALSE)
has_header <- grepl("^\\s*sample_id\\t", first_line, ignore.case = TRUE) ||
              grepl("^\\s*sample\\t", first_line, ignore.case = TRUE)

if (has_header) {
  rgi <- readr::read_delim(
    rgi_pheno_file, delim="\t", col_names=TRUE,
    show_col_types=FALSE, progress=FALSE, quote=""
  )
  if (!("sample_id" %in% names(rgi))) names(rgi)[1] <- "sample_id"
  if (!("Predicted_phenotypes" %in% names(rgi))) names(rgi)[2] <- "Predicted_phenotypes"
} else {
  rgi <- readr::read_delim(
    rgi_pheno_file, delim="\t", col_names=FALSE,
    show_col_types=FALSE, progress=FALSE, quote=""
  )
  if (ncol(rgi) < 2) stop("RGI 文件列数不足：需要至少两列 (sample_id, Predicted_phenotypes)")
  rgi <- rgi[, 1:2]
  names(rgi) <- c("sample_id", "Predicted_phenotypes")
}

rgi_long <- rgi %>%
  transmute(sample_id = .data$sample_id, ph_raw = .data$Predicted_phenotypes) %>%
  mutate(term = purrr::map(ph_raw, split_multi)) %>%
  unnest_longer(term, keep_empty=FALSE) %>%
  mutate(drug_class = sapply(term, function(x){
    y <- canonize(x)
    if (length(y)==0) NA_character_ else y
  })) %>%
  filter(!is.na(drug_class)) %>%
  distinct(sample_id, drug_class)

## ===== 读 ResFinder（可选） =====
res_long <- tibble(sample_id=character(), drug_class=character())

if (!is.na(resf_file) && file.exists(resf_file)) {
  # 关键：ResFinder.tsv 你这个文件第一行是 “## ResFinder”
  # 所以我们定位到表头 “Sample\tResistance gene...” 再开始读
  resf <- read_tsv_after_header(resf_file, header_regex = "^\\s*Sample\\t")

  # 必要列检查
  if (!("Sample" %in% names(resf))) {
    stop("ResFinder 表中没有 Sample 列。读到的列名：\n", paste(names(resf), collapse=", "))
  }
  if (!("Phenotype" %in% names(resf))) {
    stop("ResFinder 表中没有 Phenotype 列。读到的列名：\n", paste(names(resf), collapse=", "))
  }

  res_long <- resf %>%
    transmute(sample_id = .data$Sample, ph_raw = .data$Phenotype) %>%
    mutate(term = purrr::map(ph_raw, split_multi)) %>%
    unnest_longer(term, keep_empty=FALSE) %>%
    mutate(drug_class = sapply(term, function(x){
      y <- canonize(x)
      if (length(y)==0) NA_character_ else y
    })) %>%
    filter(!is.na(drug_class)) %>%
    distinct(sample_id, drug_class)
}

## ===== 合并 + 对齐元信息 =====
amr_long <- bind_rows(rgi_long, res_long) %>% distinct()

if (nrow(amr_long) == 0) {
  stop("未解析到任何药物类别；请检查 RGI/ResFinder 内容，以及 canon_classes 的匹配规则。")
}

# 用 meta 对齐（只保留 meta 里存在的样本）
amr_long2 <- amr_long %>% inner_join(meta, by="sample_id")

if (nrow(amr_long2) == 0) {
  stop("AMR 结果与 meta 的 sample_id 完全对不上。请确认：meta.tsv 的 sample_id 与 ResFinder/RGI 的样本名一致。")
}

## ===== 样本 × 类别 0/1 矩阵 =====
bin_df <- amr_long2 %>%
  mutate(value = 1L) %>%
  select(sample_id, drug_class, value) %>%
  pivot_wider(names_from = drug_class, values_from = value, values_fill = 0) %>%
  arrange(sample_id)

## ===== 每个样本的 Profile（组合） + Year 计数 =====
up_in <- bin_df %>% column_to_rownames("sample_id") %>% as.data.frame()

profile_df <- up_in %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to="drug_class", values_to="present") %>%
  filter(present == 1) %>%
  group_by(sample_id) %>%
  summarise(Profile = paste(sort(drug_class), collapse="+"), .groups="drop") %>%
  right_join(meta, by="sample_id") %>%
  mutate(Profile = tidyr::replace_na(Profile, "None"))

year_prof <- profile_df %>%
  count(Year, Profile, name="n") %>%
  arrange(Year, desc(n))

## ===== 写出 =====
write.table(profile_df, file.path(outdir, "per_sample_profile.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(year_prof, file.path(outdir, "year_profile_counts.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(amr_long2 %>% select(sample_id, drug_class) %>% distinct(),
            file.path(outdir, "_classes_used.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(bin_df, file.path(outdir, "binary_matrix.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

message("Done. Tables in: ", outdir)
