#!/usr/bin/env python3
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Circle, Patch

ROOT = "/public/home/220322096/FYM/cfr-SRR3/2499"
FIG7_DIR = f"{ROOT}/snv_cluster/fig7_inputs"
OUTDIR = f"{ROOT}/snv_cluster/fig7_plots"
os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------
# Global style
# -----------------------------
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 10
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["xtick.major.width"] = 0.8
plt.rcParams["ytick.major.width"] = 0.8

summary_file = f"{FIG7_DIR}/fig7_summary_numbers.tsv"
revised_table = f"{FIG7_DIR}/123_with_revised_source_categories.csv"

if not os.path.exists(summary_file):
    raise SystemExit(f"Missing: {summary_file}")
if not os.path.exists(revised_table):
    revised_table = f"{ROOT}/123_with_revised_source_categories.csv"
    if not os.path.exists(revised_table):
        raise SystemExit("Cannot find revised source table.")

# -----------------------------
# Read Fig. 7A summary
# -----------------------------
summary = pd.read_csv(summary_file, sep="\t", dtype=str).fillna("")
summary_dict = dict(zip(summary["metric"], summary["value"]))

H   = int(summary_dict.get("clusters_H", 0))
A   = int(summary_dict.get("clusters_A", 0))
E   = int(summary_dict.get("clusters_E", 0))
HA  = int(summary_dict.get("clusters_HA", 0))
HE  = int(summary_dict.get("clusters_HE", 0))
AE  = int(summary_dict.get("clusters_AE", 0))
HAE = int(summary_dict.get("clusters_HAE", 0))

# -----------------------------
# Read revised source table
# -----------------------------
df = pd.read_csv(revised_table, dtype=str).fillna("")

source_col = "Isolation_source_cat_revised" if "Isolation_source_cat_revised" in df.columns else "Isolation_source_cat"

orig_candidates = [
    "Original_isolation_source_from_456",
    "Original_isolation_source_from_456_y",
    "Original_isolation_source_from_456_x",
]
orig_col = None
for c in orig_candidates:
    if c in df.columns:
        orig_col = c
        break
if orig_col is None:
    raise SystemExit("Cannot find original isolation source column.")

df[source_col] = df[source_col].astype(str).str.strip().str.lower()
df[orig_col] = df[orig_col].astype(str).str.strip()

# -----------------------------
# Reviewer-friendly labels for Fig. 7B
# B panel = non-human specific sources only
# Food remains under Environmental
# No "Other animals" / "Other environment"
# -----------------------------
nonhuman = df[df[source_col].isin(["animal", "environment"])].copy()

def classify_nonhuman(row):
    cat = str(row[source_col]).strip().lower()
    s = str(row[orig_col]).strip().lower()

    # ---------- ANIMAL ----------
    if cat == "animal":
        if re.search(r'pig|swine|porcine', s):
            return ("Animal", "Pig")
        if re.search(r'cattle|cow|bovine|calf', s):
            return ("Animal", "Cattle")
        if re.search(r'chicken|broiler|hen|poultry|duck|goose', s):
            return ("Animal", "Poultry")
        if re.search(r'bird|gull|wild bird|avian', s):
            return ("Animal", "Wild birds")
        if re.search(r'feces|faeces|stool|gut|intestinal|cecal|caecal|rectal|cloacal|cecum', s):
            return ("Animal", "Gut/fecal material")

        # 这些原来容易落进 Other animals 的，尽量并入更明确的环境/生物类
        if re.search(r'fly|mactra veneriformis|shellfish|clam|oyster|mussel', s):
            return ("Environmental", "Aquatic animals/plants/environment")

        # 极少数动物组织/拭子，但不想再单独拉一个 residual animal bar
        if re.search(r'tonsil|swab', s):
            return ("Animal", "Gut/fecal material")

        # 其余难分类 animal 标签，归到环境 Unknown，避免出现 Other animals
        return ("Environmental", "Unknown environment")

    # ---------- ENVIRONMENT ----------
    if s in {"", "unknown", "nan", "not available", "missing"}:
        return ("Environmental", "Unknown environment")

    # Food stays under environment
    if re.search(
        r'food|retail|meat|beef|pork|chicken meat|duck meat|goose meat|milk|egg|'
        r'peanut butter|canned cucumber|chili|carrot juice|infant rice cereal|'
        r'vegetable|vegetables-potato|ham|honey|choy sum|comminuted chicken',
        s
    ):
        return ("Environmental", "Food")

    if re.search(r'wastewater|sewage|effluent|river|lake|water|aquatic|marine|sediment|soil|sludge|canal|sea', s):
        return ("Environmental", "Aquatic animals/plants/environment")

    if re.search(
        r'hospital|ward|icu|sink|drain|surface|environmental swab|equipment|'
        r'bed remote|bedrail|floor|sewer|healthcare system|clinical environment',
        s
    ):
        return ("Environmental", "Healthcare environment")

    if re.search(r'feed|fodder', s):
        return ("Environmental", "Feed")

    if re.search(r'manure|litter|farm|barn|pen|house dust|dust sample|slurry', s):
        return ("Environmental", "Agricultural environment")

    # 广义或模糊环境标签统一放 Unknown environment
    if re.search(r'environment$|^environmental$|environmental sample|in-house environment|surveillance', s):
        return ("Environmental", "Unknown environment")

    # 其余都并入 Unknown environment，避免 residual bar 命名含糊
    return ("Environmental", "Unknown environment")

classified = nonhuman.apply(classify_nonhuman, axis=1, result_type="expand")
classified.columns = ["plot_block", "source_group_final"]
nonhuman = pd.concat([nonhuman.reset_index(drop=True), classified.reset_index(drop=True)], axis=1)

b_counts = (
    nonhuman.groupby(["plot_block", "source_group_final"], as_index=False)
    .agg(n_isolates=("strain_id", "size"))
)

animal_order = [
    "Pig",
    "Cattle",
    "Wild birds",
    "Poultry",
    "Gut/fecal material",
]
env_order = [
    "Aquatic animals/plants/environment",
    "Food",
    "Healthcare environment",
    "Feed",
    "Agricultural environment",
    "Unknown environment",
]

animal_df = b_counts[b_counts["plot_block"] == "Animal"].copy()
env_df = b_counts[b_counts["plot_block"] == "Environmental"].copy()

animal_df["source_group_final"] = pd.Categorical(animal_df["source_group_final"], categories=animal_order, ordered=True)
env_df["source_group_final"] = pd.Categorical(env_df["source_group_final"], categories=env_order, ordered=True)

animal_df = animal_df.groupby("source_group_final", as_index=False)["n_isolates"].sum()
env_df = env_df.groupby("source_group_final", as_index=False)["n_isolates"].sum()

animal_df = animal_df.dropna(subset=["source_group_final"])
env_df = env_df.dropna(subset=["source_group_final"])

animal_df = animal_df[animal_df["n_isolates"] > 0].sort_values("source_group_final")
env_df = env_df[env_df["n_isolates"] > 0].sort_values("source_group_final")

clean_counts = pd.concat([
    animal_df.assign(source_block="Animal"),
    env_df.assign(source_block="Environmental"),
], ignore_index=True)
clean_counts.to_csv(f"{OUTDIR}/fig7B_source_group_counts_clean.tsv", sep="\t", index=False)

# -----------------------------
# Nature-like muted palette
# -----------------------------
human_fill   = "#C88996"   # muted rose
animal_fill  = "#7D96B3"   # muted slate blue
env_fill     = "#8DB6A0"   # muted sage green

animal_bar   = "#5F84A8"
env_bar      = "#6FA187"

# -----------------------------
# Helper: equal-size Venn
# -----------------------------
def draw_equal_venn(ax, H, A, E, HA, HE, AE, HAE):
    r = 1.42
    c_animal = (-0.95, 0.45)
    c_env    = ( 0.95, 0.45)
    c_human  = ( 0.00,-0.72)

    circles = [
        Circle(c_animal, r, facecolor=animal_fill, edgecolor="black", lw=0.8, alpha=0.58),
        Circle(c_env,    r, facecolor=env_fill,    edgecolor="black", lw=0.8, alpha=0.58),
        Circle(c_human,  r, facecolor=human_fill,  edgecolor="black", lw=0.8, alpha=0.58),
    ]
    for c in circles:
        ax.add_patch(c)

    ax.text(c_animal[0], c_animal[1] + r + 0.18, "Animal", ha="center", va="bottom", fontsize=11)
    ax.text(c_env[0],    c_env[1]    + r + 0.18, "Environmental", ha="center", va="bottom", fontsize=11)
    ax.text(c_human[0],  c_human[1]  - r - 0.18, "Human", ha="center", va="top", fontsize=11)

    # fixed positions so 0 also shows
    ax.text(-1.28,  0.42, str(A),   ha="center", va="center", fontsize=12)
    ax.text( 1.28,  0.42, str(E),   ha="center", va="center", fontsize=12)
    ax.text( 0.00, -1.42, str(H),   ha="center", va="center", fontsize=12)

    ax.text( 0.00,  0.86, str(AE),  ha="center", va="center", fontsize=10)
    ax.text(-0.62, -0.28, str(HA),  ha="center", va="center", fontsize=10)
    ax.text( 0.62, -0.28, str(HE),  ha="center", va="center", fontsize=10)
    ax.text( 0.00,  0.02, str(HAE), ha="center", va="center", fontsize=10)

    ax.set_xlim(-2.75, 2.75)
    ax.set_ylim(-2.55, 2.35)
    ax.set_aspect("equal")
    ax.axis("off")

# -----------------------------
# Build rows for panel B
# -----------------------------
plot_rows = []
for _, r in animal_df.iterrows():
    plot_rows.append(("Animal", str(r["source_group_final"]), int(r["n_isolates"])))

plot_rows.append(("gap", "", 0))

for _, r in env_df.iterrows():
    plot_rows.append(("Environmental", str(r["source_group_final"]), int(r["n_isolates"])))

plot_df = pd.DataFrame(plot_rows, columns=["block", "label", "n_isolates"])
plot_df["y"] = list(range(len(plot_df)))[::-1]

legend_handles = [
    Patch(facecolor=animal_bar, edgecolor="black", label="Animal"),
    Patch(facecolor=env_bar, edgecolor="black", label="Environmental"),
]

# -----------------------------
# Combined figure
# -----------------------------
fig = plt.figure(figsize=(13.0, 5.9))
gs = gridspec.GridSpec(1, 2, width_ratios=[1.02, 1.46], wspace=0.62)

# Panel A
ax1 = fig.add_subplot(gs[0, 0])
draw_equal_venn(ax1, H, A, E, HA, HE, AE, HAE)
ax1.set_title("A", loc="left", fontweight="bold", fontsize=13, pad=2)

# Panel B
ax2 = fig.add_subplot(gs[0, 1])

for _, r in plot_df.iterrows():
    if r["block"] == "gap":
        continue
    color = animal_bar if r["block"] == "Animal" else env_bar
    ax2.barh(r["y"], r["n_isolates"], height=0.72, color=color, edgecolor="black", linewidth=0.6)
    ax2.text(r["n_isolates"] + 2, r["y"], str(r["n_isolates"]), va="center", ha="left", fontsize=9)

ax2.set_yticks(plot_df[plot_df["block"] != "gap"]["y"])
ax2.set_yticklabels(plot_df[plot_df["block"] != "gap"]["label"], fontsize=10)
ax2.set_xlabel("Number of isolates", fontsize=11)
ax2.set_title("B", loc="left", fontweight="bold", fontsize=13, pad=2)
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.45)
ax2.set_axisbelow(True)
ax2.tick_params(axis="y", pad=4)

if not animal_df.empty:
    y_animal_top = plot_df[plot_df["block"] == "Animal"]["y"].max()
    ax2.text(-0.10, y_animal_top + 0.8, "Animal", transform=ax2.get_yaxis_transform(),
             fontsize=10.5, fontweight="bold", ha="right", va="center")

if not env_df.empty:
    y_env_top = plot_df[plot_df["block"] == "Environmental"]["y"].max()
    ax2.text(-0.10, y_env_top + 0.8, "Environmental", transform=ax2.get_yaxis_transform(),
             fontsize=10.5, fontweight="bold", ha="right", va="center")

ax2.legend(handles=legend_handles, title="Category", frameon=False,
           loc="upper right", fontsize=9, title_fontsize=9)

for spine in ["top", "right"]:
    ax2.spines[spine].set_visible(False)

fig.savefig(f"{OUTDIR}/Fig7_combined_final.pdf", bbox_inches="tight")
fig.savefig(f"{OUTDIR}/Fig7_combined_final.png", dpi=600, bbox_inches="tight")

# -----------------------------
# Panel A only
# -----------------------------
figA, axA = plt.subplots(figsize=(5.2, 5.0))
draw_equal_venn(axA, H, A, E, HA, HE, AE, HAE)
axA.set_title("A", loc="left", fontweight="bold", fontsize=13, pad=2)
figA.savefig(f"{OUTDIR}/Fig7A_venn_final.pdf", bbox_inches="tight")
figA.savefig(f"{OUTDIR}/Fig7A_venn_final.png", dpi=600, bbox_inches="tight")

# -----------------------------
# Panel B only
# -----------------------------
figB, axB = plt.subplots(figsize=(7.8, 5.8))

for _, r in plot_df.iterrows():
    if r["block"] == "gap":
        continue
    color = animal_bar if r["block"] == "Animal" else env_bar
    axB.barh(r["y"], r["n_isolates"], height=0.72, color=color, edgecolor="black", linewidth=0.6)
    axB.text(r["n_isolates"] + 2, r["y"], str(r["n_isolates"]), va="center", ha="left", fontsize=9)

axB.set_yticks(plot_df[plot_df["block"] != "gap"]["y"])
axB.set_yticklabels(plot_df[plot_df["block"] != "gap"]["label"], fontsize=10)
axB.set_xlabel("Number of isolates", fontsize=11)
axB.xaxis.set_major_locator(MaxNLocator(integer=True))
axB.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.45)
axB.set_axisbelow(True)
axB.tick_params(axis="y", pad=4)

if not animal_df.empty:
    y_animal_top = plot_df[plot_df["block"] == "Animal"]["y"].max()
    axB.text(-0.10, y_animal_top + 0.8, "Animal", transform=axB.get_yaxis_transform(),
             fontsize=10.5, fontweight="bold", ha="right", va="center")

if not env_df.empty:
    y_env_top = plot_df[plot_df["block"] == "Environmental"]["y"].max()
    axB.text(-0.10, y_env_top + 0.8, "Environmental", transform=axB.get_yaxis_transform(),
             fontsize=10.5, fontweight="bold", ha="right", va="center")

axB.legend(handles=legend_handles, title="Category", frameon=False,
           loc="upper right", fontsize=9, title_fontsize=9)

for spine in ["top", "right"]:
    axB.spines[spine].set_visible(False)

figB.savefig(f"{OUTDIR}/Fig7B_sources_final.pdf", bbox_inches="tight")
figB.savefig(f"{OUTDIR}/Fig7B_sources_final.png", dpi=600, bbox_inches="tight")

print("Saved:")
print(f"{OUTDIR}/Fig7_combined_final.pdf")
print(f"{OUTDIR}/Fig7_combined_final.png")
print(f"{OUTDIR}/Fig7A_venn_final.pdf")
print(f"{OUTDIR}/Fig7A_venn_final.png")
print(f"{OUTDIR}/Fig7B_sources_final.pdf")
print(f"{OUTDIR}/Fig7B_sources_final.png")
print(f"{OUTDIR}/fig7B_source_group_counts_clean.tsv")