#!/bin/bash
#SBATCH -J cfr_locus_QC
#SBATCH -o cfr_locus_QC_%j.out
#SBATCH -e cfr_locus_QC_%j.err
#SBATCH -N 1
#SBATCH -p cu
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00

set -euo pipefail

# =====================================================================
# Step 3: Extract and QC all cfr-family loci from MOB-recon carrier FASTAs
#
# Input:
#   /public/home/220322096/FYM/cfr-SRR3/2499/moball_amrfinder_results.tsv
#   /public/home/220322096/FYM/cfr-SRR3/2499/mobrecon/mob_recon_*/...
#
# Output:
#   one locus per row, strand-corrected nucleotide sequences,
#   translated proteins, codon-ready and strict-QC FASTA sets.
#
# This step is independent of Prodigal and can run while Step 2b is running.
# =====================================================================

BASE_DIR="/public/home/220322096/FYM/cfr-SRR3/2499"
AMR_TSV="${BASE_DIR}/moball_amrfinder_results.tsv"
MOB_ROOT="${BASE_DIR}/mobrecon"
OUTDIR="${BASE_DIR}/cfr_codon_compatibility/03_cfr_locus_extraction"

# Use Python directly from the existing environment; no conda activate/deactivate.
PYTHON_BIN="/public/home/220322096/miniconda3/envs/resfinder/bin/python"

mkdir -p "${OUTDIR}"

[[ -x "${PYTHON_BIN}" ]] || {
    echo "[ERROR] Python not found: ${PYTHON_BIN}" >&2
    exit 1
}
[[ -s "${AMR_TSV}" ]] || {
    echo "[ERROR] AMRFinder table missing or empty: ${AMR_TSV}" >&2
    exit 1
}
[[ -d "${MOB_ROOT}" ]] || {
    echo "[ERROR] MOB-recon root directory missing: ${MOB_ROOT}" >&2
    exit 1
}

echo "[INFO] Python      : $(${PYTHON_BIN} --version 2>&1)"
echo "[INFO] AMR table   : ${AMR_TSV}"
echo "[INFO] MOB root    : ${MOB_ROOT}"
echo "[INFO] Output      : ${OUTDIR}"
echo "[INFO] Start time  : $(date)"

"${PYTHON_BIN}" - "${AMR_TSV}" "${MOB_ROOT}" "${OUTDIR}" <<'PY'
from __future__ import print_function

import csv
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

amr_tsv = Path(sys.argv[1])
mob_root = Path(sys.argv[2])
outdir = Path(sys.argv[3])

EXPECTED_LOCI = 4085
ALLOWED_VARIANTS = {"cfr", "cfr(B)", "cfr(C)", "cfr(D)", "cfr(E)"}
VALID_STARTS_CODE11 = {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

manifest_path = outdir / "cfr_locus_manifest.tsv"
failure_path = outdir / "cfr_locus_failures.tsv"
counts_path = outdir / "cfr_locus_QC_counts.tsv"
carrier_summary_path = outdir / "carrier_processing_summary.tsv"

raw_fna_path = outdir / "all_cfr_loci_extracted_raw.fna"
codon_fna_path = outdir / "all_cfr_loci_codon_ready.fna"
strict_fna_path = outdir / "all_cfr_loci_strict_QC.fna"
codon_faa_path = outdir / "all_cfr_loci_codon_ready.faa"

CODON_TABLE = {
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

def safe_float(value):
    try:
        return float(str(value).strip())
    except Exception:
        return None

def safe_int(value):
    try:
        return int(float(str(value).strip()))
    except Exception:
        return None

def reverse_complement(seq):
    table = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn",
                          "TGCAYRKMVHDBNtgcayrkmvhdbn")
    return seq.translate(table)[::-1]

def translate_nt(seq):
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)

def parse_fasta(path):
    """
    Return:
      sequence_by_alias: aliases -> sequence
      canonical_headers: first-token IDs in original order
      contig_lengths: first-token ID -> length
    """
    sequence_by_alias = {}
    canonical_headers = []
    contig_lengths = {}

    header = None
    chunks = []

    def save_record(h, parts):
        if h is None:
            return
        seq = "".join(parts).replace(" ", "").replace("\r", "").replace("\n", "").upper()
        full = h.strip()
        token = full.split()[0]
        aliases = {full, token}
        if token.startswith("lcl|"):
            aliases.add(token[4:])
        if full.startswith("lcl|"):
            aliases.add(full[4:])
        for alias in aliases:
            if alias and alias not in sequence_by_alias:
                sequence_by_alias[alias] = seq
        canonical_headers.append(token)
        contig_lengths[token] = len(seq)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                save_record(header, chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
        save_record(header, chunks)

    return sequence_by_alias, canonical_headers, contig_lengths

def resolve_contig(contig_id, seq_by_alias):
    if contig_id in seq_by_alias:
        return seq_by_alias[contig_id], contig_id

    stripped = contig_id
    if stripped.startswith("lcl|"):
        stripped = stripped[4:]
    if stripped in seq_by_alias:
        return seq_by_alias[stripped], stripped

    # Conservative fallback: use a unique alias whose first token matches.
    candidates = []
    for alias, seq in seq_by_alias.items():
        token = alias.split()[0]
        token2 = token[4:] if token.startswith("lcl|") else token
        if token == contig_id or token2 == stripped:
            candidates.append((seq, alias))
    unique = {}
    for seq, alias in candidates:
        unique[alias] = seq
    if len(unique) == 1:
        alias = next(iter(unique))
        return unique[alias], alias

    return None, ""

def fasta_header(locus_id, isolate_id, variant, carrier_type, carrier_name,
                 contig_id, start, stop, strand, qc_label):
    return (
        locus_id
        + "|isolate=" + isolate_id
        + "|variant=" + variant
        + "|carrier_type=" + carrier_type
        + "|carrier=" + carrier_name
        + "|contig=" + contig_id
        + "|coords=" + str(start) + "-" + str(stop)
        + "|strand=" + strand
        + "|qc=" + qc_label
    )

def write_fasta_record(handle, header, seq, width=80):
    handle.write(">" + header + "\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i+width] + "\n")

# ------------------------------------------------------------------
# Read and filter cfr-family hits
# ------------------------------------------------------------------
with amr_tsv.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {
        "Sample", "Contig id", "Start", "Stop", "Strand",
        "Element symbol", "Target length", "Reference sequence length",
        "% Coverage of reference", "% Identity to reference"
    }
    missing_columns = required.difference(reader.fieldnames or [])
    if missing_columns:
        raise SystemExit(
            "[ERROR] Required columns missing from AMRFinder table: "
            + ", ".join(sorted(missing_columns))
        )

    hits = []
    for source_row_number, row in enumerate(reader, start=2):
        symbol = (row.get("Element symbol") or "").strip()
        if symbol not in ALLOWED_VARIANTS:
            continue
        row["_source_row_number"] = source_row_number
        hits.append(row)

print("[INFO] cfr-family AMRFinder rows:", len(hits))
if len(hits) != EXPECTED_LOCI:
    print("[WARNING] Expected %d loci but found %d." % (EXPECTED_LOCI, len(hits)))

# Group by carrier so each FASTA is read only once.
hits_by_carrier = defaultdict(list)
for row in hits:
    hits_by_carrier[(row.get("Sample") or "").strip()].append(row)

manifest_fields = [
    "locus_id",
    "source_AMRFinder_row",
    "sample_field",
    "isolate_id",
    "variant",
    "carrier_type",
    "carrier_name",
    "carrier_fasta",
    "contig_id_reported",
    "contig_id_resolved",
    "contig_length_nt",
    "start",
    "stop",
    "strand",
    "coordinate_length_nt",
    "extracted_length_nt",
    "length_mod_3",
    "ambiguous_nt_count",
    "start_codon",
    "valid_bacterial_start_codon",
    "terminal_stop_codon",
    "internal_stop_count",
    "translated_aa_length",
    "target_length_aa",
    "reference_length_aa",
    "target_length_match",
    "coverage_reference_pct",
    "identity_reference_pct",
    "at_contig_edge",
    "frame_QC_pass",
    "strict_QC_pass",
    "QC_status",
    "QC_reasons"
]

failure_fields = [
    "locus_id", "sample_field", "isolate_id", "variant",
    "carrier_type", "carrier_name", "carrier_fasta",
    "contig_id_reported", "start", "stop", "strand",
    "QC_status", "QC_reasons"
]

manifest_rows = []
failure_rows = []
carrier_rows = []

overall_counts = Counter()
variant_counts = defaultdict(Counter)
carrier_type_counts = defaultdict(Counter)

with raw_fna_path.open("w", encoding="utf-8") as raw_fna, \
     codon_fna_path.open("w", encoding="utf-8") as codon_fna, \
     strict_fna_path.open("w", encoding="utf-8") as strict_fna, \
     codon_faa_path.open("w", encoding="utf-8") as codon_faa:

    locus_number = 0

    for carrier_index, sample_field in enumerate(sorted(hits_by_carrier), start=1):
        carrier_hits = hits_by_carrier[sample_field]

        if "__" not in sample_field:
            folder_name = sample_field
            carrier_name = ""
            isolate_id = sample_field.replace("mob_recon_", "", 1)
            carrier_type = "unknown"
            carrier_fasta = Path("")
        else:
            folder_name, carrier_name = sample_field.rsplit("__", 1)
            isolate_id = folder_name
            if isolate_id.startswith("mob_recon_"):
                isolate_id = isolate_id[len("mob_recon_"):]
            carrier_type = "chromosome" if carrier_name == "chromosome" else (
                "plasmid" if carrier_name.startswith("plasmid_") else "unknown"
            )
            carrier_fasta = mob_root / folder_name / (carrier_name + ".fasta")

        carrier_status = "OK"
        seq_by_alias = {}
        canonical_headers = []
        contig_lengths = {}

        if not carrier_fasta.is_file() or carrier_fasta.stat().st_size == 0:
            carrier_status = "CARRIER_FASTA_MISSING"
        else:
            try:
                seq_by_alias, canonical_headers, contig_lengths = parse_fasta(carrier_fasta)
                if not seq_by_alias:
                    carrier_status = "CARRIER_FASTA_EMPTY_OR_INVALID"
            except Exception as exc:
                carrier_status = "CARRIER_FASTA_PARSE_ERROR:" + str(exc)

        print(
            "[%d/%d] %s | loci=%d | %s"
            % (carrier_index, len(hits_by_carrier), sample_field, len(carrier_hits), carrier_status),
            flush=True
        )

        carrier_rows.append({
            "sample_field": sample_field,
            "isolate_id": isolate_id,
            "carrier_type": carrier_type,
            "carrier_name": carrier_name,
            "carrier_fasta": str(carrier_fasta),
            "number_of_cfr_loci": len(carrier_hits),
            "number_of_fasta_records": len(canonical_headers),
            "carrier_status": carrier_status,
        })

        for row in carrier_hits:
            locus_number += 1
            locus_id = "CFRL%06d" % locus_number
            variant = (row.get("Element symbol") or "").strip()
            contig_id = (row.get("Contig id") or "").strip()
            start = safe_int(row.get("Start"))
            stop = safe_int(row.get("Stop"))
            strand = (row.get("Strand") or "").strip()
            target_len = safe_int(row.get("Target length"))
            reference_len = safe_int(row.get("Reference sequence length"))
            coverage = safe_float(row.get("% Coverage of reference"))
            identity = safe_float(row.get("% Identity to reference"))

            reasons = []
            qc_status = "EXTRACTED"
            resolved_id = ""
            contig_len = None
            coordinate_len = None
            extracted_len = None
            length_mod3 = None
            ambiguous_count = None
            start_codon = ""
            valid_start = None
            terminal_stop = None
            internal_stop_count = None
            aa_len = None
            target_match = None
            at_edge = None
            frame_pass = False
            strict_pass = False
            oriented_seq = ""
            aa_core = ""

            if carrier_status != "OK":
                qc_status = carrier_status
                reasons.append(carrier_status)
            elif start is None or stop is None or strand not in {"+", "-"}:
                qc_status = "INVALID_COORDINATES_OR_STRAND"
                reasons.append(qc_status)
            else:
                contig_seq, resolved_id = resolve_contig(contig_id, seq_by_alias)
                if contig_seq is None:
                    qc_status = "CONTIG_NOT_FOUND"
                    reasons.append(qc_status)
                else:
                    contig_len = len(contig_seq)
                    lo = min(start, stop)
                    hi = max(start, stop)
                    coordinate_len = hi - lo + 1

                    if lo < 1 or hi > contig_len or lo > hi:
                        qc_status = "COORDINATES_OUT_OF_RANGE"
                        reasons.append(qc_status)
                    else:
                        oriented_seq = contig_seq[lo-1:hi].upper().replace("U", "T")
                        if strand == "-":
                            oriented_seq = reverse_complement(oriented_seq)

                        extracted_len = len(oriented_seq)
                        length_mod3 = extracted_len % 3
                        ambiguous_count = sum(base not in {"A", "C", "G", "T"} for base in oriented_seq)
                        start_codon = oriented_seq[:3] if extracted_len >= 3 else ""
                        valid_start = start_codon in VALID_STARTS_CODE11

                        translated = translate_nt(oriented_seq) if length_mod3 == 0 else ""
                        terminal_stop = bool(translated.endswith("*")) if translated else False
                        aa_core = translated[:-1] if terminal_stop else translated
                        internal_stop_count = aa_core.count("*") if translated else None
                        aa_len = len(aa_core) if translated else None
                        target_match = (
                            aa_len is not None
                            and target_len is not None
                            and abs(aa_len - target_len) <= 1
                        )
                        at_edge = (lo <= 3 or hi >= contig_len - 2)

                        if coordinate_len != extracted_len:
                            reasons.append("COORDINATE_LENGTH_MISMATCH")
                        if length_mod3 != 0:
                            reasons.append("LENGTH_NOT_DIVISIBLE_BY_3")
                        if ambiguous_count != 0:
                            reasons.append("AMBIGUOUS_NUCLEOTIDES")
                        if internal_stop_count not in (None, 0):
                            reasons.append("INTERNAL_STOP_CODON")
                        if aa_len is not None and aa_len < 300:
                            reasons.append("PROTEIN_SHORTER_THAN_300_AA")
                        if target_match is False:
                            reasons.append("TARGET_LENGTH_MISMATCH")
                        if coverage is None or coverage < 95.0:
                            reasons.append("REFERENCE_COVERAGE_BELOW_95")
                        if not valid_start:
                            # Informative flag only; it does not fail frame QC because
                            # AMRFinder BLASTX coordinates may begin at an aligned codon.
                            reasons.append("NONCANONICAL_OR_UNRESOLVED_START")
                        if at_edge:
                            # Informative flag only.
                            reasons.append("LOCUS_AT_CONTIG_EDGE")

                        frame_pass = (
                            coordinate_len == extracted_len
                            and length_mod3 == 0
                            and ambiguous_count == 0
                            and internal_stop_count == 0
                            and aa_len is not None
                            and aa_len >= 300
                        )

                        strict_pass = (
                            frame_pass
                            and target_match is True
                            and coverage is not None
                            and coverage >= 95.0
                        )

                        if frame_pass:
                            qc_status = "FRAME_QC_PASS"
                        else:
                            qc_status = "FRAME_QC_FAIL"

                        header_raw = fasta_header(
                            locus_id, isolate_id, variant, carrier_type, carrier_name,
                            contig_id, start, stop, strand, qc_status
                        )
                        write_fasta_record(raw_fna, header_raw, oriented_seq)

                        if frame_pass:
                            # Exclude a terminal stop codon from codon-usage input.
                            codon_seq = oriented_seq[:-3] if terminal_stop else oriented_seq
                            header_codon = fasta_header(
                                locus_id, isolate_id, variant, carrier_type, carrier_name,
                                contig_id, start, stop, strand,
                                "STRICT_PASS" if strict_pass else "FRAME_PASS"
                            )
                            write_fasta_record(codon_fna, header_codon, codon_seq)
                            write_fasta_record(codon_faa, header_codon, aa_core)
                            if strict_pass:
                                write_fasta_record(strict_fna, header_codon, codon_seq)

            reason_text = ";".join(reasons) if reasons else "PASS"

            manifest_row = {
                "locus_id": locus_id,
                "source_AMRFinder_row": row.get("_source_row_number"),
                "sample_field": sample_field,
                "isolate_id": isolate_id,
                "variant": variant,
                "carrier_type": carrier_type,
                "carrier_name": carrier_name,
                "carrier_fasta": str(carrier_fasta),
                "contig_id_reported": contig_id,
                "contig_id_resolved": resolved_id,
                "contig_length_nt": contig_len if contig_len is not None else "",
                "start": start if start is not None else "",
                "stop": stop if stop is not None else "",
                "strand": strand,
                "coordinate_length_nt": coordinate_len if coordinate_len is not None else "",
                "extracted_length_nt": extracted_len if extracted_len is not None else "",
                "length_mod_3": length_mod3 if length_mod3 is not None else "",
                "ambiguous_nt_count": ambiguous_count if ambiguous_count is not None else "",
                "start_codon": start_codon,
                "valid_bacterial_start_codon": (
                    "Yes" if valid_start is True else ("No" if valid_start is False else "")
                ),
                "terminal_stop_codon": (
                    "Yes" if terminal_stop is True else ("No" if terminal_stop is False else "")
                ),
                "internal_stop_count": internal_stop_count if internal_stop_count is not None else "",
                "translated_aa_length": aa_len if aa_len is not None else "",
                "target_length_aa": target_len if target_len is not None else "",
                "reference_length_aa": reference_len if reference_len is not None else "",
                "target_length_match": (
                    "Yes" if target_match is True else ("No" if target_match is False else "")
                ),
                "coverage_reference_pct": coverage if coverage is not None else "",
                "identity_reference_pct": identity if identity is not None else "",
                "at_contig_edge": "Yes" if at_edge is True else ("No" if at_edge is False else ""),
                "frame_QC_pass": "Yes" if frame_pass else "No",
                "strict_QC_pass": "Yes" if strict_pass else "No",
                "QC_status": qc_status,
                "QC_reasons": reason_text,
            }
            manifest_rows.append(manifest_row)

            overall_counts["total"] += 1
            variant_counts[variant]["total"] += 1
            carrier_type_counts[carrier_type]["total"] += 1

            if oriented_seq:
                overall_counts["extracted"] += 1
                variant_counts[variant]["extracted"] += 1
                carrier_type_counts[carrier_type]["extracted"] += 1
            if frame_pass:
                overall_counts["frame_pass"] += 1
                variant_counts[variant]["frame_pass"] += 1
                carrier_type_counts[carrier_type]["frame_pass"] += 1
            if strict_pass:
                overall_counts["strict_pass"] += 1
                variant_counts[variant]["strict_pass"] += 1
                carrier_type_counts[carrier_type]["strict_pass"] += 1
            if not frame_pass:
                failure_rows.append({
                    "locus_id": locus_id,
                    "sample_field": sample_field,
                    "isolate_id": isolate_id,
                    "variant": variant,
                    "carrier_type": carrier_type,
                    "carrier_name": carrier_name,
                    "carrier_fasta": str(carrier_fasta),
                    "contig_id_reported": contig_id,
                    "start": start if start is not None else "",
                    "stop": stop if stop is not None else "",
                    "strand": strand,
                    "QC_status": qc_status,
                    "QC_reasons": reason_text,
                })

# ------------------------------------------------------------------
# Write tables
# ------------------------------------------------------------------
with manifest_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

with failure_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=failure_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(failure_rows)

carrier_fields = [
    "sample_field", "isolate_id", "carrier_type", "carrier_name",
    "carrier_fasta", "number_of_cfr_loci", "number_of_fasta_records",
    "carrier_status"
]
with carrier_summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=carrier_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(carrier_rows)

count_rows = []
count_rows.append({
    "group_type": "overall",
    "group_name": "all",
    "total_loci": overall_counts["total"],
    "successfully_extracted": overall_counts["extracted"],
    "frame_QC_pass": overall_counts["frame_pass"],
    "strict_QC_pass": overall_counts["strict_pass"],
    "frame_QC_failed": overall_counts["total"] - overall_counts["frame_pass"],
})
for variant in ["cfr", "cfr(B)", "cfr(C)", "cfr(D)", "cfr(E)"]:
    c = variant_counts[variant]
    count_rows.append({
        "group_type": "variant",
        "group_name": variant,
        "total_loci": c["total"],
        "successfully_extracted": c["extracted"],
        "frame_QC_pass": c["frame_pass"],
        "strict_QC_pass": c["strict_pass"],
        "frame_QC_failed": c["total"] - c["frame_pass"],
    })
for carrier_type in ["chromosome", "plasmid", "unknown"]:
    c = carrier_type_counts[carrier_type]
    if c["total"] == 0:
        continue
    count_rows.append({
        "group_type": "carrier_type",
        "group_name": carrier_type,
        "total_loci": c["total"],
        "successfully_extracted": c["extracted"],
        "frame_QC_pass": c["frame_pass"],
        "strict_QC_pass": c["strict_pass"],
        "frame_QC_failed": c["total"] - c["frame_pass"],
    })

count_fields = [
    "group_type", "group_name", "total_loci",
    "successfully_extracted", "frame_QC_pass",
    "strict_QC_pass", "frame_QC_failed"
]
with counts_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=count_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(count_rows)

print("=" * 76)
print("[RESULT] Total cfr-family loci      :", overall_counts["total"])
print("[RESULT] Successfully extracted    :", overall_counts["extracted"])
print("[RESULT] Frame-QC/codon-ready pass  :", overall_counts["frame_pass"])
print("[RESULT] Strict-QC pass             :", overall_counts["strict_pass"])
print("[RESULT] Frame-QC failures          :", overall_counts["total"] - overall_counts["frame_pass"])
print("[OUTPUT] Manifest                   :", manifest_path)
print("[OUTPUT] Failure table              :", failure_path)
print("[OUTPUT] QC counts                  :", counts_path)
print("[OUTPUT] Raw nucleotide FASTA       :", raw_fna_path)
print("[OUTPUT] Codon-ready nucleotide FASTA:", codon_fna_path)
print("[OUTPUT] Strict-QC nucleotide FASTA :", strict_fna_path)
print("[OUTPUT] Codon-ready protein FASTA  :", codon_faa_path)
print("=" * 76)
PY

echo "[INFO] End time: $(date)"
