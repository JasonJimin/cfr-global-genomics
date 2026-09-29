#!/usr/bin/env python3
import os
import re
import pandas as pd

ROOT = "/public/home/220322096/FYM/cfr-SRR3/2499"
CLUSTER_FILE = f"{ROOT}/snv_cluster/cfr_clustered_manifest_all.tsv"
REVISED_TABLE = f"{ROOT}/123_with_revised_source_categories.csv"
OUTDIR = f"{ROOT}/snv_cluster/fig7_inputs"
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------
# Helper functions
# ----------------------------
def pick_first_existing(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def collapse_source(x):
    s = str(x).strip().lower()

    if s in {"unknown", "", "not available", "missing", "nan"}:
        return "Unknown"

    # environment-related
    if re.search(r'wastewater|sewage|effluent|river|lake|water|aquatic|marine|sediment|soil|sludge', s):
        return "Aquatic/soil environment"
    if re.search(r'hospital|ward|icu|sink|drain|surface|environmental swab|equipment', s):
        return "Healthcare environment"
    if re.search(r'food|meat|chicken meat|pork|beef|retail', s):
        return "Food"
    if re.search(r'feed|fodder', s):
        return "Feed"

    # animal-related
    if re.search(r'pig|swine|porcine', s):
        return "Pig"
    if re.search(r'chicken|broiler|hen|poultry|duck|goose', s):
        return "Poultry"
    if re.search(r'cattle|cow|bovine|calf', s):
        return "Cattle"
    if re.search(r'sheep|goat|ovine|caprine', s):
        return "Small ruminants"
    if re.search(r'dog|canine|cat|feline|pet', s):
        return "Companion animals"
    if re.search(r'bird|gull|wild bird|avian', s):
        return "Wild birds"
    if re.search(r'feces|faeces|stool|gut|intestinal|cecal|caecal|rectal', s):
        return "Gut/fecal material"

    return "Other"

# ----------------------------
# A. cluster-level source overlap
# ----------------------------
cl = pd.read_csv(CLUSTER_FILE, sep="\t", dtype=str).fillna("")

source_col_A = pick_first_existing(cl, ["Isolation_source_cat_revised", "Isolation_source_cat"])
if source_col_A is None:
    raise SystemExit("No source category column found in cluster file.")

cl[source_col_A] = cl[source_col_A].str.strip().str.lower()
cl = cl[cl["cfr_cluster_id"].notna() & (cl["cfr_cluster_id"] != "")].copy()
cl = cl[cl[source_col_A].isin(["human", "animal", "environment"])].copy()

def overlap_label(vals):
    s = sorted(set(vals))
    return " + ".join(s) if s else "unknown"

cluster_overlap = (
    cl.groupby("cfr_cluster_id", as_index=False)
      .agg(
          group_id=("group_id", "first"),
          cfr_variant_types=("cfr_variant_types", "first"),
          n_isolates=("strain_id", "size"),
          n_human=(source_col_A, lambda x: (x == "human").sum()),
          n_animal=(source_col_A, lambda x: (x == "animal").sum()),
          n_environment=(source_col_A, lambda x: (x == "environment").sum()),
          source_membership=(source_col_A, overlap_label)
      )
)

cluster_overlap["membership_code"] = cluster_overlap.apply(
    lambda r: "".join([
        "H" if r["n_human"] > 0 else "",
        "A" if r["n_animal"] > 0 else "",
        "E" if r["n_environment"] > 0 else ""
    ]),
    axis=1
)

cluster_overlap = cluster_overlap.sort_values(
    ["membership_code", "n_isolates", "cfr_cluster_id"],
    ascending=[True, False, True]
)

cluster_overlap.to_csv(f"{OUTDIR}/fig7A_cluster_source_overlap.tsv", sep="\t", index=False)

overlap_counts = (
    cluster_overlap.groupby(["membership_code", "source_membership"], as_index=False)
    .agg(
        n_clusters=("cfr_cluster_id", "nunique"),
        n_isolates=("n_isolates", "sum")
    )
    .sort_values(["n_clusters", "source_membership"], ascending=[False, True])
)

overlap_counts.to_csv(f"{OUTDIR}/fig7A_cluster_source_overlap_counts.tsv", sep="\t", index=False)

variant_overlap_counts = (
    cluster_overlap.groupby(["cfr_variant_types", "membership_code"], as_index=False)
    .agg(n_clusters=("cfr_cluster_id", "nunique"))
    .sort_values(["cfr_variant_types", "n_clusters"], ascending=[True, False])
)

variant_overlap_counts.to_csv(f"{OUTDIR}/fig7A_variant_specific_overlap_counts.tsv", sep="\t", index=False)

summary_lines = []
summary_lines.append(("total_clusters_with_source", cluster_overlap["cfr_cluster_id"].nunique()))
for code in ["H", "A", "E", "HA", "HE", "AE", "HAE"]:
    summary_lines.append((f"clusters_{code}", int((cluster_overlap["membership_code"] == code).sum())))

pd.DataFrame(summary_lines, columns=["metric", "value"]).to_csv(
    f"{OUTDIR}/fig7_summary_numbers.tsv", sep="\t", index=False
)

# ----------------------------
# B. specific non-human sources
# Use the full revised table, not cluster manifest
# ----------------------------
if not os.path.exists(REVISED_TABLE):
    raise SystemExit(f"Missing revised source table: {REVISED_TABLE}")

full = pd.read_csv(REVISED_TABLE, dtype=str).fillna("")

source_col_B = pick_first_existing(full, ["Isolation_source_cat_revised", "Isolation_source_cat"])
orig_col_B = pick_first_existing(
    full,
    [
        "Original_isolation_source_from_456",
        "Original_isolation_source_from_456_y",
        "Original_isolation_source_from_456_x"
    ]
)

if source_col_B is None:
    raise SystemExit("No source category column found in revised table.")
if orig_col_B is None:
    raise SystemExit("No original source column found in revised table.")

full[source_col_B] = full[source_col_B].str.strip().str.lower()
full[orig_col_B] = full[orig_col_B].fillna("").astype(str).str.strip()

nonhuman = full[full[source_col_B].isin(["animal", "environment"])].copy()
nonhuman["original_source_raw"] = nonhuman[orig_col_B].replace("", "unknown")
nonhuman["source_group_for_plot"] = nonhuman["original_source_raw"].map(collapse_source)

raw_counts = (
    nonhuman.groupby([source_col_B, "original_source_raw"], as_index=False)
    .agg(n_isolates=("strain_id", "size"))
    .sort_values(["n_isolates", "original_source_raw"], ascending=[False, True])
)

raw_counts.to_csv(f"{OUTDIR}/fig7B_nonhuman_original_source_raw_counts.tsv", sep="\t", index=False)

group_counts = (
    nonhuman.groupby([source_col_B, "source_group_for_plot"], as_index=False)
    .agg(n_isolates=("strain_id", "size"))
    .sort_values([source_col_B, "n_isolates"], ascending=[True, False])
)

group_counts.to_csv(f"{OUTDIR}/fig7B_nonhuman_source_group_counts.tsv", sep="\t", index=False)

print("Wrote corrected Fig.7 input tables to:", OUTDIR)
