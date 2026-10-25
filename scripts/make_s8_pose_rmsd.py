#!/usr/bin/env python3
"""Table S8: 本文 2 段落目のポーズ RMSD（1.02 Å と 0.64 Å）を再導出する。

  python3 scripts/make_s8_pose_rmsd.py      # numpy が要る

出力: data/tables/S8_pose_rmsd.tsv

HADDOCK は run ごとに座標系が違うので、座標を直に比べてはいけない。
共通の受容体（鎖 A の Cα）で Kabsch 重ね合わせしてから、リガンド重原子の RMSD を測る
（`scripts/pose_rmsd.py` と同じ手続き。重ね合わせの実装はそこから読み込む）。

この表が示すこと
----------------
* 無拘束 run の最良クラスタは、拘束あり run の最良クラスタと **1.02 Å** で一致する
* 電荷モデルを AM1-BCC から Gasteiger に替えると、
  **最良クラスタ同士は 1.64 Å 離れる**が、**対応するクラスタは 0.64 Å しか動かない**。
  動いたのはポーズではなく**順位**で、そのクラスタは 1 位から 3 位に落ちる
  （MMG_chainB の capri_clt: cluster_rank 1→cluster_id 2、cluster_rank 3→cluster_id 1）
"""
import glob
import importlib.util
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "data", "runs")
OUT = os.path.join(ROOT, "data", "tables", "S8_pose_rmsd.tsv")

_spec = importlib.util.spec_from_file_location(
    "pr", os.path.join(ROOT, "scripts", "pose_rmsd.py"))
pr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pr)


def clusters(run):
    """cluster_rank -> (cluster_id, score) を capri_clt.tsv から読む。"""
    p = os.path.join(RUNS, run, "capri_clt.tsv")
    L = [l for l in open(p) if not l.startswith("#")]
    h = L[0].rstrip("\n").split("\t")
    out = {}
    for l in L[1:]:
        d = dict(zip(h, l.rstrip("\n").split("\t")))
        out[d["cluster_rank"]] = (d["cluster_id"], d["score"])
    return out


def compare(a, b):
    A, B = pr.read(a), pr.read(b)
    ca_a = {x["resnum"]: x["xyz"] for x in A if x["chain"] == "A" and x["name"] == "CA"}
    ca_b = {x["resnum"]: x["xyz"] for x in B if x["chain"] == "A" and x["name"] == "CA"}
    com = sorted(set(ca_a) & set(ca_b))
    P = np.array([ca_a[r] for r in com])
    Q = np.array([ca_b[r] for r in com])
    tf = pr.kabsch(P, Q)
    fit = float(np.sqrt((((np.array([tf(p) for p in P])) - Q) ** 2).sum(1).mean()))
    la = [x for x in A if x["chain"] == "B" and x["el"] != "H"]
    lb = {x["name"]: x["xyz"] for x in B if x["chain"] == "B" and x["el"] != "H"}
    pairs = [(x["xyz"], lb[x["name"]]) for x in la if x["name"] in lb]
    LA = np.array([tf(p) for p, _ in pairs])
    LB = np.array([q for _, q in pairs])
    rms = float(np.sqrt(((LA - LB) ** 2).sum(1).mean()))
    cen = float(np.linalg.norm(LA.mean(0) - LB.mean(0)))
    return len(com), fit, rms, cen, len(pairs)


def path(run, rank, model=1):
    return os.path.join(RUNS, run, f"cluster_{rank}_model_{model}.pdb.gz")


# (説明, A run, A rank, B run, B rank, 本文のどこ)
JOBS = [
    ("blind vs restrained, both top-ranked",
     "MMT_blind", 1, "MMT_chainB", 1,
     "Letter, paragraph 2: 'agreed to 1.02 A with a restrained docking'"),
    # 電荷モデルを替えたとき、AM1-BCC の最良ポーズが Gasteiger のどのクラスタに
    # 対応するかを全クラスタと比べて示す。対応先は rank 3（= cluster_id 1）である。
    ("AM1-BCC top vs Gasteiger rank 1", "MMT_chainB", 1, "MMG_chainB", 1,
     "the two best-scoring poses are NOT the same pose"),
    ("AM1-BCC top vs Gasteiger rank 2", "MMT_chainB", 1, "MMG_chainB", 2, ""),
    ("AM1-BCC top vs Gasteiger rank 3", "MMT_chainB", 1, "MMG_chainB", 3,
     "Letter, paragraph 2: 'an equivalent pose 0.64 A away, though it fell from first to third'"),
    ("AM1-BCC top vs Gasteiger rank 4", "MMT_chainB", 1, "MMG_chainB", 4, ""),
]


def main():
    cl = {r: clusters(r) for r in {j[1] for j in JOBS} | {j[3] for j in JOBS}}
    rows = []
    for lab, ra, ka, rb, kb, note in JOBS:
        n, fit, rms, cen, na = compare(path(ra, ka), path(rb, kb))
        ida, sca = cl[ra][str(ka)]
        idb, scb = cl[rb][str(kb)]
        rows.append([lab, ra, ka, ida, sca, rb, kb, idb, scb,
                     n, f"{fit:.3f}", na, f"{rms:.3f}", f"{cen:.3f}", note])

    head = ["comparison", "run_A", "A_cluster_rank", "A_cluster_id", "A_cluster_score",
            "run_B", "B_cluster_rank", "B_cluster_id", "B_cluster_score",
            "n_CA_superposed", "receptor_fit_rmsd_A", "n_ligand_heavy_atoms",
            "ligand_rmsd_A", "ligand_centroid_shift_A", "note"]
    with open(OUT, "w") as f:
        f.write("# Table S8. Ligand RMSD after Kabsch superposition on the common receptor "
                "(chain A CA atoms).\n")
        f.write("# Model 1 of each cluster. RMSD is symmetry-unaware.\n")
        f.write("# Cluster files are named by rank; the cluster_id column shows that the "
                "charge model\n")
        f.write("# reorders the clusters without moving the pose.\n")
        for r in ("MMT_chainB", "MMG_chainB"):
            order = "; ".join(f"rank {k} = cluster_id {v[0]} (score {v[1]})"
                              for k, v in sorted(cl[r].items(), key=lambda x: int(x[0])))
            f.write(f"# {r}: {order}\n")
        f.write("\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(str(c) for c in r) + "\n")
    print(f"{OUT} に {len(rows)} 行を書いた\n")
    for r in rows:
        print(f"  {r[0]:52s} ligand RMSD {r[12]:>6s} Å   (receptor fit {r[10]} Å, {r[9]} CA)")
    print()
    for r in ("MMT_chainB", "MMG_chainB"):
        print(f"  {r}: " + ", ".join(f"rank{k}=id{v[0]}"
                                     for k, v in sorted(cl[r].items(), key=lambda x: int(x[0]))))


if __name__ == "__main__":
    main()
