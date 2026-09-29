#!/bin/bash
#SBATCH -J codo
#SBATCH -o cfr_codon_profiles_%j.out
#SBATCH -e cfr_codon_profiles_%j.err
#SBATCH -N 1
#SBATCH -p cu
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=2-00:00:00

set -euo pipefail

# =====================================================================
# Step 4: Calculate three codon-usage profiles
#   1) cfr-family gene profile (one row per QC-passed cfr locus)
#   2) host chromosome profile (one row per isolate chromosome)
#   3) cfr-carrying replicon profile (one row per cfr-positive
#      chromosome or plasmid)
#
# Metrics:
#   - 61 sense-codon counts
#   - 61 RSCU values
#   - GC content and GC3
#
# To avoid circular similarity, Prodigal CDSs overlapping any cfr locus
# are excluded from host/carrier profiles.
# =====================================================================

BASE_DIR="/public/home/220322096/FYM/cfr-SRR3/2499"
PROJECT_DIR="${BASE_DIR}/cfr_codon_compatibility"

LOCI_DIR="${PROJECT_DIR}/03_cfr_locus_extraction"
PRODIGAL_DIR="${PROJECT_DIR}/02_mob_carrier_prodigal"
OUTDIR="${PROJECT_DIR}/04_codon_profiles"

LOCI_MANIFEST="${LOCI_DIR}/cfr_locus_manifest.tsv"
CFR_CODON_FASTA="${LOCI_DIR}/all_cfr_loci_codon_ready.fna"
PRODIGAL_SUMMARY="${PRODIGAL_DIR}/prodigal_summary_clean.tsv"

PYTHON_BIN="/public/home/220322096/miniconda3/envs/resfinder/bin/python"

mkdir -p "${OUTDIR}"

for f in "${LOCI_MANIFEST}" "${CFR_CODON_FASTA}" "${PRODIGAL_SUMMARY}"; do
    [[ -s "${f}" ]] || {
        echo "[ERROR] Required input missing or empty: ${f}" >&2
        exit 1
    }
done

[[ -x "${PYTHON_BIN}" ]] || {
    echo "[ERROR] Python not found: ${PYTHON_BIN}" >&2
    exit 1
}

echo "[INFO] Python            : $(${PYTHON_BIN} --version 2>&1)"
echo "[INFO] Locus manifest    : ${LOCI_MANIFEST}"
echo "[INFO] cfr codon FASTA   : ${CFR_CODON_FASTA}"
echo "[INFO] Prodigal summary  : ${PRODIGAL_SUMMARY}"
echo "[INFO] Output directory  : ${OUTDIR}"
echo "[INFO] Start time        : $(date)"

"${PYTHON_BIN}" - \
    "${LOCI_MANIFEST}" \
    "${CFR_CODON_FASTA}" \
    "${PRODIGAL_SUMMARY}" \
    "${OUTDIR}" <<'PY'
from __future__ import print_function

import csv
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

locus_manifest_path = Path(sys.argv[1])
cfr_fasta_path = Path(sys.argv[2])
prodigal_summary_path = Path(sys.argv[3])
outdir = Path(sys.argv[4])

# ---------------------------------------------------------------------
# Standard bacterial genetic code (NCBI translation table 11 uses the
# same amino-acid assignments as the standard code; alternative starts
# do not affect synonymous codon families).
# Stop codons are excluded from all codon-usage profiles.
# ---------------------------------------------------------------------
CODON_TO_AA = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

SENSE_CODONS = sorted(c for c, aa in CODON_TO_AA.items() if aa != "*")
INFORMATIVE_CODONS = [c for c in SENSE_CODONS if c not in {"ATG", "TGG"}]

AA_TO_CODONS = defaultdict(list)
for codon in SENSE_CODONS:
    AA_TO_CODONS[CODON_TO_AA[codon]].append(codon)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def normalize_contig_id(value):
    value = (value or "").strip().split()[0]
    if value.startswith("lcl|"):
        value = value[4:]
    return value

def parse_fasta(path):
    header = None
    seq_parts = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper().replace("U", "T")
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts).upper().replace("U", "T")

def parse_gff(path):
    """
    Map Prodigal CDS ID to:
      contig, start, end, strand, partial
    """
    result = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            contig = normalize_contig_id(fields[0])
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            strand = fields[6]
            attributes = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = value
            gene_id = attributes.get("ID", "")
            if not gene_id:
                continue
            result[gene_id] = {
                "contig": contig,
                "start": min(start, end),
                "end": max(start, end),
                "strand": strand,
                "partial": attributes.get("partial", ""),
            }
    return result

def extract_prodigal_id(header):
    match = re.search(r"(?:^|[;\s])ID=([^;\s]+)", header)
    if match:
        return match.group(1)
    return ""

def intervals_overlap(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

def prepare_coding_sequence(seq):
    """
    Return (clean_sequence_without_terminal_stop, reason).
    Only complete, unambiguous, in-frame CDSs are retained.
    """
    seq = seq.upper().replace("U", "T")
    if not seq:
        return None, "EMPTY"
    if any(base not in {"A", "C", "G", "T"} for base in seq):
        return None, "AMBIGUOUS_NT"
    if len(seq) % 3 != 0:
        return None, "NOT_DIVISIBLE_BY_3"

    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in {"TAA", "TAG", "TGA"}:
        codons = codons[:-1]

    if not codons:
        return None, "NO_SENSE_CODONS"
    if any(CODON_TO_AA.get(codon) == "*" for codon in codons):
        return None, "INTERNAL_STOP"
    if any(codon not in CODON_TO_AA for codon in codons):
        return None, "UNKNOWN_CODON"

    return "".join(codons), "PASS"

def count_codons(sequences):
    counts = Counter()
    total_nt = 0
    total_gc = 0
    total_third = 0
    total_gc3 = 0

    for seq in sequences:
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            if codon not in SENSE_CODONS:
                continue
            counts[codon] += 1
            total_nt += 3
            total_gc += codon.count("G") + codon.count("C")
            total_third += 1
            if codon[2] in {"G", "C"}:
                total_gc3 += 1

    total_codons = sum(counts.values())
    gc = float(total_gc) / total_nt if total_nt else ""
    gc3 = float(total_gc3) / total_third if total_third else ""
    return counts, total_codons, gc, gc3

def calculate_rscu(counts):
    rscu = {}
    for aa, codons in AA_TO_CODONS.items():
        total = sum(counts.get(codon, 0) for codon in codons)
        n_synonymous = len(codons)
        expected = float(total) / n_synonymous if total else 0.0
        for codon in codons:
            if expected > 0:
                rscu[codon] = counts.get(codon, 0) / expected
            else:
                rscu[codon] = 0.0
    return rscu

def profile_columns():
    return (
        ["total_sense_codons", "GC", "GC3"]
        + ["COUNT_" + codon for codon in SENSE_CODONS]
        + ["RSCU_" + codon for codon in SENSE_CODONS]
    )

def profile_values(sequences):
    counts, total_codons, gc, gc3 = count_codons(sequences)
    rscu = calculate_rscu(counts)
    values = {
        "total_sense_codons": total_codons,
        "GC": round(gc, 6) if gc != "" else "",
        "GC3": round(gc3, 6) if gc3 != "" else "",
    }
    for codon in SENSE_CODONS:
        values["COUNT_" + codon] = counts.get(codon, 0)
    for codon in SENSE_CODONS:
        values["RSCU_" + codon] = round(rscu[codon], 6)
    return values

# ---------------------------------------------------------------------
# Read cfr locus manifest.
# Build:
#   locus metadata
#   exclusion intervals by carrier and contig
# ---------------------------------------------------------------------
with locus_manifest_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
    locus_rows = list(csv.DictReader(handle, delimiter="\t"))

locus_by_id = {}
exclusions = defaultdict(lambda: defaultdict(list))
cfr_positive_carriers = set()

for row in locus_rows:
    locus_id = (row.get("locus_id") or "").strip()
    if locus_id:
        locus_by_id[locus_id] = row

    carrier_id = (row.get("sample_field") or "").strip()
    if carrier_id:
        cfr_positive_carriers.add(carrier_id)

    contig = normalize_contig_id(
        (row.get("contig_id_resolved") or "").strip()
        or (row.get("contig_id_reported") or "").strip()
    )
    try:
        start = int(float(row.get("start") or ""))
        stop = int(float(row.get("stop") or ""))
    except ValueError:
        continue

    if carrier_id and contig:
        exclusions[carrier_id][contig].append(
            (min(start, stop), max(start, stop), locus_id)
        )

# ---------------------------------------------------------------------
# 1. cfr gene profiles
# ---------------------------------------------------------------------
cfr_output = outdir / "cfr_gene_codon_profiles.tsv"
cfr_fields = [
    "locus_id", "isolate_id", "variant", "carrier_type", "carrier_name",
    "strict_QC_pass", "sequence_length_nt"
] + profile_columns()

cfr_rows = []
cfr_skipped = []

for header, raw_seq in parse_fasta(cfr_fasta_path):
    locus_id = header.split("|", 1)[0].split()[0]
    meta = locus_by_id.get(locus_id, {})

    clean_seq, reason = prepare_coding_sequence(raw_seq)
    if clean_seq is None:
        cfr_skipped.append((locus_id, reason))
        continue

    row = {
        "locus_id": locus_id,
        "isolate_id": meta.get("isolate_id", ""),
        "variant": meta.get("variant", ""),
        "carrier_type": meta.get("carrier_type", ""),
        "carrier_name": meta.get("carrier_name", ""),
        "strict_QC_pass": meta.get("strict_QC_pass", ""),
        "sequence_length_nt": len(clean_seq),
    }
    row.update(profile_values([clean_seq]))
    cfr_rows.append(row)

with cfr_output.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=cfr_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(cfr_rows)

# ---------------------------------------------------------------------
# Read clean Prodigal summary and verify there are no failures.
# ---------------------------------------------------------------------
with prodigal_summary_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
    prodigal_rows = list(csv.DictReader(handle, delimiter="\t"))

bad_status = [
    row for row in prodigal_rows
    if (row.get("status") or "").strip() not in {"OK", "SKIPPED_EXISTING"}
]
if bad_status:
    raise SystemExit(
        "[ERROR] prodigal_summary_clean.tsv still contains %d failed rows."
        % len(bad_status)
    )

host_rows = []
carrier_rows = []
profile_failures = []

host_fields = [
    "carrier_id", "isolate_id", "carrier_type", "carrier_name",
    "ffn_file", "gff_file",
    "CDS_total", "CDS_used", "CDS_partial_excluded",
    "CDS_cfr_overlap_excluded", "CDS_invalid_excluded",
    "basic_profile_pass", "strict_profile_pass"
] + profile_columns()

carrier_fields = list(host_fields)

# ---------------------------------------------------------------------
# 2 & 3. Host chromosome and cfr-carrier profiles
# ---------------------------------------------------------------------
for index, prow in enumerate(prodigal_rows, start=1):
    carrier_id = (prow.get("carrier_id") or "").strip()
    isolate_id = (prow.get("isolate_id") or "").strip()
    carrier_type = (prow.get("carrier_type") or "").strip()
    carrier_name = (prow.get("carrier_name") or "").strip()
    ffn_path = Path((prow.get("ffn_file") or "").strip())
    gff_path = Path((prow.get("gff_file") or "").strip())

    print(
        "[%d/%d] %s" % (index, len(prodigal_rows), carrier_id),
        flush=True
    )

    if not ffn_path.is_file() or not gff_path.is_file():
        profile_failures.append({
            "carrier_id": carrier_id,
            "reason": "FFN_OR_GFF_MISSING",
            "ffn_file": str(ffn_path),
            "gff_file": str(gff_path),
        })
        continue

    gff_map = parse_gff(gff_path)
    sequences_used = []

    cds_total = 0
    cds_used = 0
    cds_partial = 0
    cds_cfr_overlap = 0
    cds_invalid = 0

    for header, seq in parse_fasta(ffn_path):
        cds_total += 1
        gene_id = extract_prodigal_id(header)
        feature = gff_map.get(gene_id)

        if feature is None:
            cds_invalid += 1
            continue

        # Exclude incomplete Prodigal CDS predictions.
        partial = feature.get("partial", "")
        if partial and partial != "00":
            cds_partial += 1
            continue

        # Exclude any CDS overlapping a cfr-family locus in the same carrier.
        overlap = False
        contig = normalize_contig_id(feature["contig"])
        for cfr_start, cfr_end, _locus_id in exclusions[carrier_id].get(contig, []):
            if intervals_overlap(
                feature["start"], feature["end"], cfr_start, cfr_end
            ):
                overlap = True
                break
        if overlap:
            cds_cfr_overlap += 1
            continue

        clean_seq, reason = prepare_coding_sequence(seq)
        if clean_seq is None:
            cds_invalid += 1
            continue

        sequences_used.append(clean_seq)
        cds_used += 1

    metrics = profile_values(sequences_used)

    # Transparent QC flags. All profiles are retained; flags are used later.
    total_codons = metrics["total_sense_codons"]
    if carrier_type == "chromosome":
        basic_pass = cds_used >= 100 and total_codons >= 10000
        strict_pass = cds_used >= 500 and total_codons >= 100000
    else:
        basic_pass = cds_used >= 5 and total_codons >= 300
        strict_pass = cds_used >= 10 and total_codons >= 1000

    base_row = {
        "carrier_id": carrier_id,
        "isolate_id": isolate_id,
        "carrier_type": carrier_type,
        "carrier_name": carrier_name,
        "ffn_file": str(ffn_path),
        "gff_file": str(gff_path),
        "CDS_total": cds_total,
        "CDS_used": cds_used,
        "CDS_partial_excluded": cds_partial,
        "CDS_cfr_overlap_excluded": cds_cfr_overlap,
        "CDS_invalid_excluded": cds_invalid,
        "basic_profile_pass": "Yes" if basic_pass else "No",
        "strict_profile_pass": "Yes" if strict_pass else "No",
    }
    base_row.update(metrics)

    # Every chromosome is a host profile.
    if carrier_type == "chromosome":
        host_rows.append(dict(base_row))

    # Only replicons actually carrying a cfr-family locus are carrier profiles.
    if carrier_id in cfr_positive_carriers:
        carrier_rows.append(dict(base_row))

host_output = outdir / "host_chromosome_codon_profiles.tsv"
carrier_output = outdir / "cfr_carrier_codon_profiles.tsv"

with host_output.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=host_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(host_rows)

with carrier_output.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=carrier_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(carrier_rows)

# ---------------------------------------------------------------------
# QC summary
# ---------------------------------------------------------------------
summary_output = outdir / "codon_profile_QC_summary.tsv"
failure_output = outdir / "codon_profile_failures.tsv"
informative_output = outdir / "informative_RSCU_columns.txt"

def count_flag(rows, field, value="Yes"):
    return sum(1 for row in rows if row.get(field) == value)

summary_rows = [
    {
        "profile_type": "cfr_gene",
        "total_profiles": len(cfr_rows),
        "basic_profile_pass": len(cfr_rows),
        "strict_profile_pass": sum(
            1 for row in cfr_rows if row.get("strict_QC_pass") == "Yes"
        ),
        "failed_or_skipped": len(cfr_skipped),
    },
    {
        "profile_type": "host_chromosome",
        "total_profiles": len(host_rows),
        "basic_profile_pass": count_flag(host_rows, "basic_profile_pass"),
        "strict_profile_pass": count_flag(host_rows, "strict_profile_pass"),
        "failed_or_skipped": sum(
            1 for row in prodigal_rows
            if (row.get("carrier_type") or "").strip() == "chromosome"
        ) - len(host_rows),
    },
    {
        "profile_type": "cfr_carrier_all",
        "total_profiles": len(carrier_rows),
        "basic_profile_pass": count_flag(carrier_rows, "basic_profile_pass"),
        "strict_profile_pass": count_flag(carrier_rows, "strict_profile_pass"),
        "failed_or_skipped": len(cfr_positive_carriers) - len(carrier_rows),
    },
    {
        "profile_type": "cfr_carrier_chromosome",
        "total_profiles": sum(
            1 for row in carrier_rows if row["carrier_type"] == "chromosome"
        ),
        "basic_profile_pass": sum(
            1 for row in carrier_rows
            if row["carrier_type"] == "chromosome"
            and row["basic_profile_pass"] == "Yes"
        ),
        "strict_profile_pass": sum(
            1 for row in carrier_rows
            if row["carrier_type"] == "chromosome"
            and row["strict_profile_pass"] == "Yes"
        ),
        "failed_or_skipped": 0,
    },
    {
        "profile_type": "cfr_carrier_plasmid",
        "total_profiles": sum(
            1 for row in carrier_rows if row["carrier_type"] == "plasmid"
        ),
        "basic_profile_pass": sum(
            1 for row in carrier_rows
            if row["carrier_type"] == "plasmid"
            and row["basic_profile_pass"] == "Yes"
        ),
        "strict_profile_pass": sum(
            1 for row in carrier_rows
            if row["carrier_type"] == "plasmid"
            and row["strict_profile_pass"] == "Yes"
        ),
        "failed_or_skipped": 0,
    },
]

with summary_output.open("w", encoding="utf-8", newline="") as handle:
    fields = [
        "profile_type", "total_profiles", "basic_profile_pass",
        "strict_profile_pass", "failed_or_skipped"
    ]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

with failure_output.open("w", encoding="utf-8", newline="") as handle:
    fields = ["profile_type", "profile_id", "reason", "ffn_file", "gff_file"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for locus_id, reason in cfr_skipped:
        writer.writerow({
            "profile_type": "cfr_gene",
            "profile_id": locus_id,
            "reason": reason,
            "ffn_file": "",
            "gff_file": "",
        })
    for row in profile_failures:
        writer.writerow({
            "profile_type": "host_or_carrier",
            "profile_id": row["carrier_id"],
            "reason": row["reason"],
            "ffn_file": row["ffn_file"],
            "gff_file": row["gff_file"],
        })

with informative_output.open("w", encoding="utf-8") as handle:
    for codon in INFORMATIVE_CODONS:
        handle.write("RSCU_" + codon + "\n")

print("=" * 78)
print("[RESULT] cfr gene profiles             :", len(cfr_rows))
print("[RESULT] host chromosome profiles       :", len(host_rows))
print("[RESULT] cfr-positive carrier profiles  :", len(carrier_rows))
print("[RESULT] profile failures               :", len(profile_failures))
print("[OUTPUT] cfr gene profiles              :", cfr_output)
print("[OUTPUT] host chromosome profiles        :", host_output)
print("[OUTPUT] cfr carrier profiles            :", carrier_output)
print("[OUTPUT] QC summary                      :", summary_output)
print("[OUTPUT] informative RSCU columns        :", informative_output)
print("=" * 78)
PY

echo "[INFO] End time: $(date)"
