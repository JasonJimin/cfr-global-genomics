#!/usr/bin/env python3
import os
import glob
import pandas as pd
import networkx as nx

ROOT = "/public/home/220322096/FYM/cfr-SRR3/2499"
OUTROOT = f"{ROOT}/snv_cluster/all_groups"

all_rows = []

for group_dir in sorted(glob.glob(f"{OUTROOT}/*")):
    if not os.path.isdir(group_dir):
        continue

    manifest_file = os.path.join(group_dir, "group_manifest.csv")
    dist_file = os.path.join(group_dir, "distances.distances.tsv")
    thr_file = os.path.join(group_dir, "threshold_info.txt")

    if not (os.path.exists(manifest_file) and os.path.exists(dist_file) and os.path.exists(thr_file)):
        continue

    manifest = pd.read_csv(manifest_file, dtype=str).fillna("")
    d = pd.read_csv(dist_file, sep="\t", dtype=str).fillna("")

    info = {}
    with open(thr_file) as f:
        for line in f:
            k, v = line.rstrip("\n").split("\t", 1)
            info[k] = v

    group_id = info["group_id"]
    cutoff = float(info["approx_snp_cutoff"])

    sample_ids = set()
    for x in manifest["matched_fasta"]:
        bn = os.path.basename(x)
        sid = bn.replace(".fasta", "").replace(".fa", "").replace(".fna", "")
        sample_ids.add(sid)

    G = nx.Graph()
    G.add_nodes_from(sample_ids)

    d["SNPs_num"] = pd.to_numeric(d["SNPs"], errors="coerce")

    keep = d[d["SNPs_num"] <= cutoff].copy()

    for _, r in keep.iterrows():
        s1 = r["Sample 1"]
        s2 = r["Sample 2"]
        if s1 in sample_ids and s2 in sample_ids:
            G.add_edge(s1, s2)

    components = list(nx.connected_components(G))
    comp_rows = []

    for i, comp in enumerate(sorted(components, key=lambda x: (-len(x), sorted(x)[0])), start=1):
        cid = f"{group_id}__cfr_cluster_{i}"
        for sid in sorted(comp):
            comp_rows.append({"sample_id_for_ska": sid, "cfr_cluster_id": cid, "group_id": group_id})

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(os.path.join(group_dir, "cfr_clusters_5snv_perMb.tsv"), sep="\t", index=False)

    merged = manifest.copy()
    merged["sample_id_for_ska"] = merged["matched_fasta"].map(lambda x: os.path.basename(x).replace(".fasta","").replace(".fa","").replace(".fna",""))
    merged = merged.merge(comp_df, on=["sample_id_for_ska", "group_id"], how="left")

    merged.to_csv(os.path.join(group_dir, "group_manifest_with_cfr_clusters.tsv"), sep="\t", index=False)

    all_rows.append(merged)

if all_rows:
    all_df = pd.concat(all_rows, ignore_index=True)
    all_df.to_csv(f"{ROOT}/snv_cluster/cfr_clustered_manifest_all.tsv", sep="\t", index=False)

    summary = (
        all_df.groupby(["group_id", "cfr_cluster_id"], as_index=False)
        .agg(
            n=("strain_id", "size"),
            n_human=("Isolation_source_cat", lambda s: (s == "human").sum()),
            n_animal=("Isolation_source_cat", lambda s: (s == "animal").sum()),
            n_environment=("Isolation_source_cat", lambda s: (s == "environment").sum())
        )
        .sort_values(["group_id", "n"], ascending=[True, False])
    )
    summary.to_csv(f"{ROOT}/snv_cluster/cfr_cluster_summary.tsv", sep="\t", index=False)

    print("Wrote:")
    print(f"{ROOT}/snv_cluster/cfr_clustered_manifest_all.tsv")
    print(f"{ROOT}/snv_cluster/cfr_cluster_summary.tsv")
else:
    print("No group outputs found.")
