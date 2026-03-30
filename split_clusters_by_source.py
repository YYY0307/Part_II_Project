#!/usr/bin/env python3
"""
Split clusters_with_pass.tsv into:
  1) clusters_with_pass_mch_pantera_only.tsv
  2) clusters_with_pass_pepi_berg_only.tsv
  3) clusters_with_pass_mixed.tsv  (both groups present)

Within each cluster, rows are ordered as:
  1) Flybase(Berg) passed
  2) Gonzalez(Pepi) passed
  3) Pantera passed
  4) MCH passed
  5) All failed (any source, original order among them)
"""

import pandas as pd

# Read clusters_with_pass from previous step
df = pd.read_csv("clusters_with_pass.tsv", sep="\t")

# Restrict to passed entries to determine which sources are represented
passed = df[df["is_pass"]].copy()

name_lower_passed = passed["name"].astype(str).str.lower()

# Group flags per TE (for presence classification)
passed["is_groupA"] = name_lower_passed.str.startswith("pepi") | name_lower_passed.str.startswith("berg")
passed["is_groupB"] = name_lower_passed.str.startswith("mch")  | name_lower_passed.str.startswith("pantera")

# Cluster-level flags: does this cluster have any passed TE from each group?
groupA_by_cluster = passed.groupby("cluster_id")["is_groupA"].any().rename("cluster_has_groupA")
groupB_by_cluster = passed.groupby("cluster_id")["is_groupB"].any().rename("cluster_has_groupB")

cluster_flags = pd.concat([groupA_by_cluster, groupB_by_cluster], axis=1).fillna(False)

# Merge flags back to all rows of those clusters
df = df.merge(cluster_flags, on="cluster_id", how="left").fillna(
    {"cluster_has_groupA": False, "cluster_has_groupB": False}
)


# Define within-cluster sort key: Flybase pass > Gonzalez pass > mch pass > pantera pass > all fails
name_lower_all = df["name"].astype(str).str.lower()
df["orig_idx"] = df.index  # to preserve order among fails

def sort_key(row) -> int:
    if row["is_pass"]:
        n = str(row["name"]).lower()
        if n.startswith("berg"):
            return 0
        if n.startswith("pepi"):
            return 1
        if n.startswith("pantera"):
            return 2
        if n.startswith("mch"):
            return 3
    # all fails go last
    return 5

df["sort_group"] = df.apply(sort_key, axis=1).astype(int)

# Sort globally by cluster_id, then sort_group, then original index
df_sorted = df.sort_values(
    by=["cluster_id", "sort_group", "orig_idx"],
    ascending=[True, True, True]
).reset_index(drop=True)


# Split into categories (clusters; keep all members of each cluster)

only_mch_pantera = df_sorted[(df_sorted["cluster_has_groupB"]) & (~df_sorted["cluster_has_groupA"])]
only_pepi_berg   = df_sorted[(df_sorted["cluster_has_groupA"]) & (~df_sorted["cluster_has_groupB"])]
mixed            = df_sorted[(df_sorted["cluster_has_groupA"]) & (df_sorted["cluster_has_groupB"])]

# Write outputs
only_mch_pantera.to_csv("clusters_with_pass_mch_pantera_only.tsv", sep="\t", index=False)
only_pepi_berg.to_csv("clusters_with_pass_pepi_berg_only.tsv", sep="\t", index=False)
mixed.to_csv("clusters_with_pass_mixed.tsv", sep="\t", index=False)