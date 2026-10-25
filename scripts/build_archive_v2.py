#!/usr/bin/env python3
"""レター v2 の Supplementary Data アーカイブを組む。

  python3 scripts/build_archive_v2.py

出力: supplementary/Supplementary_Data_v2/ と同名の .tar.gz

v1 の `Supplementary_Data/`（2026-09-11）は欠陥リガンド時代のもので、v2 の数値と合わない。
こちらが v2 用。本文が依存する 97 run だけを入れる（探索的な run は入れない）。
"""
import glob
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "data", "runs")
OUT = os.path.join(ROOT, "supplementary", "Supplementary_Data_v2")
SEEDS = ["", "_i20000", "_i40000", "_i60000", "_i80000"]

LADDER = ["Y3B", "Y3BD184N", "Y3BD184A", "Y3BD184E", "Y3BD184E2",
          "Y3BD184HP", "Y3BD184HN", "Y3BD184HP2",
          "MMAREP_MMT_blind", "MMAREP_MMT_D184N_blind", "MMAREP_MMT_D184A_blind",
          "MMAREP_MMT_E164A_blind",
          "D184H_rot1P_MMT_blind", "D184H_rot1N_MMT_blind", "D184H_rot2P_MMT_blind",
          "Y3MWT", "Y3ME", "Y3ME2"]
METHODS = ["MMT_blind", "MMT_blind_seed1000", "MMT_blind_seed2000", "MMT_blind_seed3000",
           "MMT_chainB", "MMG_chainB", "PSB-172656_blind"]
TABLES = ["S1_run_inventory.tsv", "S2_score_terms_wt.tsv", "S3_D184_ladder.tsv",
          "S4_contact_residues_wt.tsv", "S5_distance_to_184.tsv",
          "S6_hbond_geometry.tsv", "S7_jaccard_per_pose.tsv",
          "S8_pose_rmsd.tsv", "S9_cryoem_contacts.tsv",
          # verify.py は v1 の表も読む（レター v1 からの連続性を保つため）
          "T1_energies.tsv", "T4_contacts_per_ligand_atom.tsv",
          "T6_seed_replicates.tsv", "T7_D184_dissection.tsv",
          "T8_binding_mode_comparison.tsv", "T9_D184_ladder.tsv"]
SCRIPTS = ["make_tables.py", "make_ladder_table.py", "make_si_tables.py",
           "make_s6_hbond_geometry.py", "make_s8_pose_rmsd.py", "pose_rmsd.py",
           "build_c4880_trimer.py", "build_charge_mutants.py",
           "build_d184e_rotamers.py", "build_d184h_rotamers.py",
           "gen_d184_ladder.py", "n_to_sidechain.py",
           "verify.py", "build_archive_v2.py"]
PROTOCOL = ["PREREG_TRIMER.md", "PREREG_replication.md", "PREREG_D184H.md"]

README = """Supplementary Data for

  "Two routes into one pocket: compound 48/80 engages MRGPRX2 through charge and
   monomethyl methacrylate through van der Waals contact"

Contents
--------
data/runs/       every run cited in the Letter: its HADDOCK3 configuration, the
                 per-cluster scores, the cluster assignment of all refined
                 models, and the four models of the rank-1 cluster
data/raw/        the receptor, the variant receptors and the ligands as docked
data/params/     ACPYPE parameter and topology files for every ligand
data/restraints/ the restraint files of the restrained reference run
data/tables/     Supplementary Tables S1-S9 as tab-separated values
cache/           7vv6.cif and 7vdm.cif, read by Table S9
scripts/         the analysis code
protocol/        the protocol written before the runs were submitted.
                 NOTE: section 16 of PREREG_TRIMER.md is mislabelled as written
                 before submission; the correction is at the top of that section.

The layout mirrors the project repository, so the scripts below run unchanged
from the top of this archive.

Runs
----
%(nrun)d runs: 18 conditions x 5 replicate seeds for the D184 substitution series,
plus %(nmeth)d runs used only in the methods (the unrestrained replicates, the
restrained reference run, its Gasteiger-charge counterpart, and the antagonist
control).

Each run directory holds
  <run>.cfg                  the HADDOCK3 configuration
  capri_clt.tsv              per-cluster scores from caprieval
  cluster_assignment.tsv     which cluster each refined model fell into
  cluster_1_model_[1-4].pdb.gz  the four models of the rank-1 cluster

The cluster_<N>_model_<M> files are named by cluster RANK, not by cluster id.
Where the two differ, Table S1 gives both.

Regenerating the numbers
------------------------
With Python 3 and NumPy, from the directory that contains data/ and scripts/:

  python3 scripts/make_ladder_table.py      # T9 = S3
  python3 scripts/make_si_tables.py         # S1, S2, S4, S5, S7, S9
  python3 scripts/make_s6_hbond_geometry.py # S6
  python3 scripts/make_s8_pose_rmsd.py      # S8  (needs NumPy)
  python3 scripts/verify.py                 # checks every quoted number against
                                            # the tables (46 items)

scripts/make_tables.py regenerates the earlier tables T1-T8. Those draw on runs
from the first campaign that are not part of this deposit, so T1-T8 are included
as generated; verify.py reads them for the checks that the two versions share.

HADDOCK3 (2025.11.0) is required only to repeat the docking itself.

Checksums for every file are in MANIFEST.sha256.
"""


def copy(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


def main():
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    runs = [p + s for p in LADDER for s in SEEDS] + METHODS
    missing = [r for r in runs if not os.path.isdir(os.path.join(RUNS, r))]
    if missing:
        sys.exit(f"run が無い: {missing}")

    n_model = 0
    for r in runs:
        for f in sorted(glob.glob(os.path.join(RUNS, r, "*"))):
            b = os.path.basename(f)
            # Table S8 は MMG_chainB の全クラスタ代表を使うので、
            # cluster_1 だけでなく harvest 済みのクラスタをすべて入れる
            if b.startswith("cluster_") and b.endswith(".pdb.gz"):
                n_model += 1
            elif not (b.endswith(".cfg") or b in ("capri_clt.tsv", "cluster_assignment.tsv")):
                continue
            # リポジトリと同じ配置にする。そうしないとスクリプトがそのまま走らない
            copy(f, os.path.join(OUT, "data", "runs", r, b))

    for sub in ("raw", "params", "restraints"):
        for f in sorted(glob.glob(os.path.join(ROOT, "data", sub, "*"))):
            copy(f, os.path.join(OUT, "data", sub, os.path.basename(f)))
    for f in sorted(glob.glob(os.path.join(ROOT, "cache", "7v*.cif"))):   # S9 が読む
        copy(f, os.path.join(OUT, "cache", os.path.basename(f)))
    for t in TABLES:
        copy(os.path.join(ROOT, "data", "tables", t), os.path.join(OUT, "data", "tables", t))
    for s in SCRIPTS:
        copy(os.path.join(ROOT, "scripts", s), os.path.join(OUT, "scripts", s))
    for p in PROTOCOL:
        copy(os.path.join(ROOT, p), os.path.join(OUT, "protocol", p))

    with open(os.path.join(OUT, "README.txt"), "w") as f:
        f.write(README % {"nrun": len(runs), "nmeth": len(METHODS)})

    files = sorted(f for f in glob.glob(os.path.join(OUT, "**", "*"), recursive=True)
                   if os.path.isfile(f))
    with open(os.path.join(OUT, "MANIFEST.sha256"), "w") as f:
        for p in files:
            h = hashlib.sha256(open(p, "rb").read()).hexdigest()
            f.write(f"{h}  {os.path.relpath(p, OUT)}\n")

    tgz = OUT + ".tar.gz"
    if os.path.exists(tgz):
        os.remove(tgz)
    with tarfile.open(tgz, "w:gz") as t:
        t.add(OUT, arcname=os.path.basename(OUT))

    total = sum(os.path.getsize(p) for p in files)
    print(f"{OUT}")
    print(f"  run {len(runs)} 本（ラダー {len(runs)-len(METHODS)} + 方法 {len(METHODS)}）、"
          f"構造 {n_model} 本、ファイル {len(files)+1} 個、{total/1e6:.1f} MB")
    print(f"  -> {tgz}  {os.path.getsize(tgz)/1e6:.1f} MB")


if __name__ == "__main__":
    main()
