#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
from collections import Counter, defaultdict, OrderedDict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

VARIANT_ORDER = ["cfr", "cfr(B)", "cfr(C)", "cfr(D)", "cfr(E)"]
CATEGORY_ORDER = [
    "Integron/IS/ICE",
    "Integron/IS",
    "Integron/ICE",
    "Integron",
    "IS/ICE",
    "IS",
    "ICE",
    "Unclassified",
]
COLORS = OrderedDict([
    ("Integron/IS/ICE", "#c44e52"),
    ("Integron/IS", "#dd8452"),
    ("Integron/ICE", "#ccb974"),
    ("Integron", "#e0c596"),
    ("IS/ICE", "#55a868"),
    ("IS", "#9fbf84"),
    ("ICE", "#4c9f9f"),
    ("Unclassified", "#9a9a9a"),
])

def norm_variant(v: str) -> str:
    return {"Cfr(D)": "cfr(D)", "Cfr(E)": "cfr(E)"}.get(str(v), str(v))

def norm_sample(s: str) -> str:
    s = str(s)
    s = re.sub(r"^mob_recon_", "", s)
    s = re.sub(r"__(chromosome|plasmid.*)$", "", s)
    return s

def parse_header_map(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            seq_id, header = parts
            # expected style:
            # mob_recon_ERR10014996__plasmid_AE822|cfr|contig=60|gene=123-456(+)|window=...|source_len=...
            toks = header.split("|")
            sample = toks[0]
            variant = toks[1] if len(toks) > 1 else None
            m_contig = re.search(r"(?:^|\|)contig=([^|]+)", header)
            m_gene = re.search(r"(?:^|\|)gene=(\d+)-(\d+)\(([+-])\)", header)

            contig = m_contig.group(1) if m_contig else None
            cfr_start = int(m_gene.group(1)) if m_gene else None
            cfr_end = int(m_gene.group(2)) if m_gene else None
            strand = m_gene.group(3) if m_gene else None

            rows.append({
                "Seq_ID": seq_id,
                "Sample": sample,
                "Sample_norm": norm_sample(sample),
                "cfr_variant": norm_variant(variant),
                "Contig_id": str(contig) if contig is not None else None,
                "cfr_start": cfr_start,
                "cfr_end": cfr_end,
                "Strand": strand,
                "Original_header": header,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No records parsed from header map.")
    required = ["Seq_ID", "Sample", "Sample_norm", "cfr_variant", "Contig_id", "cfr_start", "cfr_end"]
    miss = [c for c in required if c not in df.columns]
    if miss:
        raise RuntimeError(f"Header map missing parsed columns: {miss}")
    return df

def parse_integron_summary(path: str) -> pd.DataFrame:
    header = None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            if line.startswith("ID_replicon\t"):
                header = line.split("\t")
                continue
            if header is None:
                continue
            vals = line.split("\t")
            if len(vals) < len(header):
                continue
            row = dict(zip(header, vals))
            seq_id = row["ID_replicon"]
            calin = int(row.get("CALIN", 0))
            complete = int(row.get("complete", 0))
            ino = int(row.get("In0", 0))
            present = any(x > 0 for x in [calin, complete, ino])
            if complete > 0:
                subtype = "complete integron"
            elif calin > 0:
                subtype = "CALIN"
            elif ino > 0:
                subtype = "In0"
            else:
                subtype = "none"
            rows.append({
                "Seq_ID": seq_id,
                "integron_present": present,
                "integron_subtype": subtype,
                "integron_CALIN": calin,
                "integron_complete": complete,
                "integron_In0": ino,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No records parsed from integron summary.")
    return df

def detect_is_from_attr(attr_text: str) -> bool:
    s = attr_text.lower()
    keywords = [
        "transposase",
        "insertion sequence",
        "transposon",
    ]
    if any(k in s for k in keywords):
        return True
    # common transposase / IS labels
    if re.search(r"\bis[0-9a-z][0-9a-z\-_/]*\b", s):
        return True
    if re.search(r"\btnp[a-z0-9]*\b", s):
        return True
    if re.search(r"\btn[0-9][0-9a-z\-_/]*\b", s):
        return True
    return False

def parse_prokka_gff(path: str) -> pd.DataFrame:
    seen_seq_ids = set()
    rows = []
    matched_products = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seq_id, source, feature, start, end, score, strand, phase, attrs = parts
            if feature != "CDS":
                continue
            if detect_is_from_attr(attrs):
                matched_products[seq_id].append(attrs)
                seen_seq_ids.add(seq_id)
    for seq_id in seen_seq_ids:
        rows.append({"Seq_ID": seq_id, "IS_present": True, "IS_hit_count": len(matched_products[seq_id]),
                     "IS_examples": " || ".join(matched_products[seq_id][:5])})
    df = pd.DataFrame(rows)
    if df.empty:
        # allow empty; later fill False
        return pd.DataFrame(columns=["Seq_ID", "IS_present", "IS_hit_count", "IS_examples"])
    return df

def parse_ice_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for c in ["Sample", "cfr_variant", "Contig_id"]:
        df[c] = df[c].astype(str)
    df["Sample_norm"] = df["Sample"].map(norm_sample)
    df["cfr_variant"] = df["cfr_variant"].map(norm_variant)
    # coerce numeric cols
    for c in ["cfr_start", "cfr_end"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    df["ICE_present"] = df["ICE_same_contig"].astype(str).str.lower().eq("yes")
    return df

def assign_category(row) -> str:
    parts = []
    if bool(row.get("integron_present", False)):
        parts.append("Integron")
    if bool(row.get("IS_present", False)):
        parts.append("IS")
    if bool(row.get("ICE_present", False)):
        parts.append("ICE")
    return "/".join(parts) if parts else "Unclassified"

def make_plot(wide_df: pd.DataFrame, out_png: str, out_pdf: str, n_total: int):
    fig, ax = plt.subplots(figsize=(13.5, 5.4))
    y_labels = [v for v in VARIANT_ORDER if v in set(wide_df["cfr_variant"])]
    data = wide_df.set_index("cfr_variant").reindex(y_labels).fillna(0)

    totals = data[CATEGORY_ORDER].sum(axis=1)
    left = [0.0] * len(data)
    legend_labels = []

    for cat in CATEGORY_ORDER:
        vals = data[cat].tolist()
        props = [(v / t if t else 0) for v, t in zip(vals, totals)]
        if sum(vals) == 0:
            # keep zero classes in legend with 0 count, but don't draw visible bars
            legend_labels.append((cat, f"{cat} (0)"))
            continue
        bars = ax.barh(
            range(len(data)),
            props,
            left=left,
            color=COLORS[cat],
            edgecolor="white",
            linewidth=0.8,
            label=f"{cat} ({sum(vals)})",
        )
        legend_labels.append((cat, f"{cat} ({sum(vals)})"))
        for i, (bar, v, p) in enumerate(zip(bars, vals, props)):
            if v <= 0:
                continue
            if p >= 0.035 or v >= 20:
                ax.text(
                    left[i] + p / 2,
                    bar.get_y() + bar.get_height() / 2,
                    str(int(v)),
                    ha="center",
                    va="center",
                    fontsize=10,
                )
        left = [l + p for l, p in zip(left, props)]

    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data.index.tolist(), fontsize=12)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of qualified cfr-flanking sequences", fontsize=12)
    ax.set_ylabel("cfr variant", fontsize=12)
    ax.set_title(f"Associated mobile elements by cfr variant (n={n_total})", fontsize=14)

    # Build legend in fixed order, including zeros
    handles = []
    labels = []
    for cat in CATEGORY_ORDER:
        from matplotlib.patches import Patch
        handles.append(Patch(facecolor=COLORS[cat], edgecolor="white"))
        total_cat = int(data[cat].sum()) if cat in data.columns else 0
        labels.append(f"{cat} ({total_cat})")
    ax.legend(handles, labels, title="Associated mobile element", bbox_to_anchor=(1.02, 0.5),
              loc="center left", frameon=False, fontsize=11, title_fontsize=12)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Make corrected cfr variant stacked MGE plot with fixed ICE sample matching")
    ap.add_argument("--header-map", required=True)
    ap.add_argument("--integron-summary", required=True)
    ap.add_argument("--prokka-gff", required=True)
    ap.add_argument("--ice-table", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    header_df = parse_header_map(args.header_map)
    int_df = parse_integron_summary(args.integron_summary)
    is_df = parse_prokka_gff(args.prokka_gff)
    ice_df = parse_ice_table(args.ice_table)

    merged = header_df.merge(int_df, on="Seq_ID", how="left")
    merged = merged.merge(is_df, on="Seq_ID", how="left")

    # exact ICE match on normalized sample + variant + contig + coordinates
    ice_key_cols = ["Sample_norm", "cfr_variant", "Contig_id", "cfr_start", "cfr_end"]
    merged = merged.merge(
        ice_df[ice_key_cols + ["ICE_present", "ICE_same_contig", "ICE_overlap", "Min_distance_bp", "ICE_ref", "ICE_hit_start", "ICE_hit_end", "ICE_aln_len", "ICE_identity"]],
        on=ice_key_cols,
        how="left",
    )

    # fill missing flags
    merged["integron_present"] = merged["integron_present"].fillna(False).astype(bool)
    merged["integron_subtype"] = merged["integron_subtype"].fillna("none")
    merged["IS_present"] = merged["IS_present"].fillna(False).astype(bool)
    merged["IS_hit_count"] = merged["IS_hit_count"].fillna(0).astype(int)
    merged["IS_examples"] = merged["IS_examples"].fillna("")
    merged["ICE_present"] = merged["ICE_present"].fillna(False).astype(bool)
    merged["ICE_same_contig"] = merged["ICE_same_contig"].fillna("No")
    merged["ICE_overlap"] = merged["ICE_overlap"].fillna("No")

    merged["MGE_category"] = merged.apply(assign_category, axis=1)

    # QC
    dup_n = merged.duplicated(subset=["Sample", "cfr_variant", "Contig_id", "cfr_start", "cfr_end"]).sum()
    exact_ice_yes = int(merged["ICE_present"].sum())

    count_df = (
        merged.groupby(["cfr_variant", "MGE_category"])
        .size()
        .reset_index(name="n")
    )
    count_df["cfr_variant"] = pd.Categorical(count_df["cfr_variant"], categories=VARIANT_ORDER, ordered=True)
    count_df["MGE_category"] = pd.Categorical(count_df["MGE_category"], categories=CATEGORY_ORDER, ordered=True)
    count_df = count_df.sort_values(["cfr_variant", "MGE_category"])

    wide_df = (
        merged.groupby(["cfr_variant", "MGE_category"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    for cat in CATEGORY_ORDER:
        if cat not in wide_df.columns:
            wide_df[cat] = 0
    wide_df["cfr_variant"] = pd.Categorical(wide_df["cfr_variant"], categories=VARIANT_ORDER, ordered=True)
    wide_df = wide_df.sort_values("cfr_variant")
    wide_df["total"] = wide_df[CATEGORY_ORDER].sum(axis=1)

    out_prefix = args.out_prefix
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    merged.to_csv(out_prefix + ".per_sequence.tsv", sep="\t", index=False)
    count_df.to_csv(out_prefix + ".count_table.tsv", sep="\t", index=False)
    wide_df.to_csv(out_prefix + ".wide_table.tsv", sep="\t", index=False)

    summary_lines = []
    summary_lines.append(f"Qualified sequences: {len(merged)}")
    summary_lines.append(f"Duplicate biological keys (Sample+cfr_variant+Contig_id+cfr_start+cfr_end): {dup_n}")
    summary_lines.append(f"Exact ICE same-contig matches after fixed sample normalization: {exact_ice_yes}")
    summary_lines.append("")
    summary_lines.append("By cfr variant:")
    for v in VARIANT_ORDER:
        n = int((merged["cfr_variant"] == v).sum())
        if n > 0:
            summary_lines.append(f"  {v}\t{n}")
    summary_lines.append("")
    summary_lines.append("MGE categories:")
    total = len(merged)
    cat_counts = merged["MGE_category"].value_counts()
    for cat in CATEGORY_ORDER:
        n = int(cat_counts.get(cat, 0))
        summary_lines.append(f"  {cat}\t{n}\t{n/total*100:.1f}%")
    summary_lines.append("")
    summary_lines.append("Integron subtype among qualified sequences:")
    subtype_counts = merged["integron_subtype"].value_counts()
    for subtype in ["In0", "complete integron", "CALIN", "none"]:
        n = int(subtype_counts.get(subtype, 0))
        summary_lines.append(f"  {subtype}\t{n}")

    with open(out_prefix + ".summary.txt", "w") as fh:
        fh.write("\n".join(summary_lines) + "\n")

    make_plot(wide_df, out_prefix + ".png", out_prefix + ".pdf", len(merged))

if __name__ == "__main__":
    main()
