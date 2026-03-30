#!/usr/bin/env python3
import sys
import re
import pandas as pd

"""
Map both all_pass.stats and all_failed.stats onto CD-HIT clusters.
Clusters with >=1 passed entry go to clusters_with_pass.tsv (contains ALL stats).
Clusters with 0 passed entries go to clusters_without_pass.tsv (contains ALL stats).

Usage:
  python passing_clusters.py all_pass.stats all_failed.stats prefixed_30_clustered_Universal_80_80.clstr

Outputs:
  clusters_with_pass.tsv
  clusters_without_pass.tsv
  (both with new 'status' column: 'passed', 'failed', or 'no_stats')
"""

def parse_clstr(clstr_path):
    cluster_id = None
    rows = []

    with open(clstr_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                cluster_id = line.split()[1]
            else:
                m = re.search(r">(.+?)\.\.\.", line)
                if m and cluster_id is not None:
                    header = m.group(1).strip()
                    name = header.lstrip('>')
                    rows.append((cluster_id, name))

    return pd.DataFrame(rows, columns=["cluster_id", "name"])

def main():
    if len(sys.argv) != 4:
        print("Usage: python passing_clusters.py all_pass.stats all_failed.stats in.clstr")
        sys.exit(1)

    pass_path = sys.argv[1]
    failed_path = sys.argv[2]
    clstr_path = sys.argv[3]

    # 1 Read stats tables separately with source tracking
    passed = pd.read_csv(pass_path, sep="\t")
    passed["name_clean"] = passed["name"].astype(str)
    passed["source"] = "passed"
    
    failed = pd.read_csv(failed_path, sep="\t")
    failed["name_clean"] = failed["name"].astype(str)
    failed["source"] = "failed"

    # 2 Parse clusters
    clusters = parse_clstr(clstr_path)
    clusters["name_clean"] = clusters["name"].astype(str)

    # 3 Create comprehensive stats lookup
    all_stats = pd.concat([passed, failed], ignore_index=True, sort=False)
    stats_lookup = all_stats.set_index("name_clean")["source"].to_dict()

    # 4 Merge ALL stats onto clusters
    merged = clusters.merge(all_stats.drop(columns=["name", "source"]), 
                           on="name_clean", how="left")
    
    # 5 Add status column
    merged["status"] = merged["name_clean"].map(stats_lookup).fillna("no_stats")

    # 6 Identify passed entries for cluster-level decision
    merged["is_pass"] = merged["status"] == "passed"
    
    # 7 Cluster-level flag
    cluster_has_pass = merged.groupby("cluster_id")["is_pass"].any().rename("cluster_has_pass")
    merged = merged.merge(cluster_has_pass, on="cluster_id", how="left")

    # 8 Split by cluster decision (KEEP status column)
    with_pass = merged[merged["cluster_has_pass"]].copy()
    without_pass = merged[~merged["cluster_has_pass"]].copy()

    # Reorder columns (status near name for visibility)
    cols = ['cluster_id', 'name', 'name_clean', 'status', 'is_pass', 'cluster_has_pass'] + \
           [col for col in merged.columns if col not in ['cluster_id', 'name', 'name_clean', 'status', 'is_pass', 'cluster_has_pass']]
    with_pass = with_pass[cols]
    without_pass = without_pass[cols]

    # Save
    with_pass.to_csv("clusters_with_pass.tsv", sep="\t", index=False)
    without_pass.to_csv("clusters_without_pass.tsv", sep="\t", index=False)

    print("Status breakdown:")
    print(merged["status"].value_counts())
    print("\nTotal clusters:", merged["cluster_id"].nunique())
    print("Passed entries:", (merged["status"] == "passed").sum())
    print("Failed entries:", (merged["status"] == "failed").sum())
    print("No stats entries:", (merged["status"] == "no_stats").sum())
    print("Clusters with ≥1 passed:", with_pass["cluster_id"].nunique())
    print("Clusters w/o passed:", without_pass["cluster_id"].nunique())
    print("Wrote clusters_with_pass.tsv and clusters_without_pass.tsv")

if __name__ == "__main__":
    main()


