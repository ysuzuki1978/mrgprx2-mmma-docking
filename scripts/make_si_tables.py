#!/usr/bin/env python3
"""レター v2 の Supplementary Tables S1・S2・S4・S5・S7・S9 を生成する。

  python3 scripts/make_si_tables.py      # letter_MMA/ 直下で実行

残りの表:
  S3 = data/tables/T9_D184_ladder.tsv      （`make_ladder_table.py`）
  S6 = data/tables/S6_hbond_geometry.tsv   （`make_s6_hbond_geometry.py`）
  S8 = data/tables/S8_pose_rmsd.tsv        （`make_s8_pose_rmsd.py`、numpy が要る）

接触・同一部位・電荷寄与の定義は `make_ladder_table.py` から読み込む。
定義を 2 か所に書かないため（本文・Figure 2・SI がすべて同じ閾値で動く）。

このスクリプトは標準ライブラリだけで動く。
"""
import csv
import glob
import importlib.util
import math
import os
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "data", "runs")
TAB = os.path.join(ROOT, "data", "tables")
CACHE = os.path.join(ROOT, "cache")

_spec = importlib.util.spec_from_file_location(
    "ml", os.path.join(ROOT, "scripts", "make_ladder_table.py"))
ml = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ml)          # main() は __main__ ガードの中なので走らない

SEEDS = ml.SEEDS
CONTACT = ml.CONTACT
SAME_SITE = ml.SAME_SITE
W = ml.W_EMREF
BYRESN = ml.BYRESN

# (リガンド表示名, 受容体表示名, run 接頭辞)。T9 と同じ 18 条件。
LADDER = [(lig, rec, pre) for lig, rec, pre, _ in ml.SETS]
# 方法の記述だけに使う 7 run（反復の接尾辞を持たない）
METHODS_RUNS = [
    ("MMMA", "WT", "MMT_blind", "無拘束。本文の収束率 379/400 の出所"),
    ("MMMA", "WT", "MMT_blind_seed1000", "同上、独立 seed"),
    ("MMMA", "WT", "MMT_blind_seed2000", "同上、独立 seed"),
    ("MMMA", "WT", "MMT_blind_seed3000", "同上、独立 seed"),
    ("MMMA", "WT", "MMT_chainB", "拘束あり。1.02 Å の比較相手"),
    ("MMMA", "WT", "MMG_chainB", "拘束あり・Gasteiger 電荷。0.64 Å の比較相手"),
    ("PSB-172656", "WT", "PSB-172656_blind", "無拘束。拮抗薬の対照"),
]
HOSTS = [("/home/yakuri/trimer_check", "yakuri3 (trimer_check)"),
         ("/home/yakuri/haddock3_transfer", "yakuri3 (haddock3_transfer)"),
         ("/home/yasu1978/haddock3_nig", "NIG (epyc)")]


def write(name, header, rows, notes=()):
    p = os.path.join(TAB, name)
    with open(p, "w") as f:
        for n in notes:
            f.write(f"# {n}\n")
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if c is None else str(c) for c in r) + "\n")
    print(f"  {name:28s} {len(rows):4d} 行")
    return p


def cfg_host(run):
    for p in glob.glob(os.path.join(RUNS, run, "*.cfg")):
        t = open(p).read()
        for pre, lab in HOSTS:
            if pre in t:
                return lab
    return "?"


def cfg_iniseed(run):
    for p in glob.glob(os.path.join(RUNS, run, "*.cfg")):
        for l in open(p):
            if l.startswith("iniseed"):
                return l.split("=")[1].strip()
        return "917"          # cfg に行が無い＝既定値
    return "?"


def assignment(run):
    """cluster_assignment.tsv から (総模型数, 未クラスタ数, {cluster_id: 個数})。"""
    p = os.path.join(RUNS, run, "cluster_assignment.tsv")
    if not os.path.exists(p):
        return None
    n, unc, cnt = 0, 0, {}
    for i, l in enumerate(open(p)):
        if i == 0:
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 4 or not f[3]:
            continue
        n += 1
        if f[3] == "-":
            unc += 1
        else:
            cnt[f[3]] = cnt.get(f[3], 0) + 1
    return n, unc, cnt


def poses(pre):
    """(seed ラベル, model 番号, パス) を 20 個返す。"""
    out = []
    for s in SEEDS:
        for m in sorted(glob.glob(os.path.join(RUNS, pre + s, "cluster_1_model_*.pdb.gz"))):
            out.append((s.lstrip("_i") or "917", os.path.basename(m).split("_model_")[1][0], m))
    return out


def contact_set(lig, prot):
    heavy = [a[4] for a in lig if a[3] != "H"]
    return {(a[0], a[1]) for a in prot
            if min(math.dist(a[4], h) for h in heavy) <= CONTACT}


# ---------------------------------------------------------------- S1
def s1():
    rows = []
    for lig, rec, pre in LADDER:
        for s in SEEDS:
            run = pre + s
            a = assignment(run)
            r = ml.rank1(run)
            if not a or not r:
                continue
            n, unc, cnt = a
            cid = r.get("cluster_id")
            pop = cnt.get(cid, 0)
            rows.append([lig, rec, run, cfg_host(run), cfg_iniseed(run),
                         n, unc, cid, pop, f"{100*pop/n:.1f}",
                         r["score"], r["vdw"], r["elec"], r["desolv"], "D184 ladder"])
    for lig, rec, run, note in METHODS_RUNS:
        a = assignment(run)
        r = ml.rank1(run)
        if not a or not r:
            continue
        n, unc, cnt = a
        cid = r.get("cluster_id")
        pop = cnt.get(cid, 0)
        rows.append([lig, rec, run, cfg_host(run), cfg_iniseed(run),
                     n, unc, cid, pop, f"{100*pop/n:.1f}",
                     r["score"], r["vdw"], r["elec"], r["desolv"], note])
    return write("S1_run_inventory.tsv",
                 ["ligand", "receptor", "run", "host", "iniseed",
                  "models_refined", "unclustered", "rank1_cluster_id",
                  "rank1_population", "rank1_pct_of_models",
                  "score", "vdw", "elec", "desolv", "role"],
                 rows,
                 ["Table S1. Every HADDOCK3 run cited in the Letter.",
                  "rank1_population counts models of the top-ranked cluster, identified by",
                  "cluster_rank==1 in capri_clt.tsv; seletopclusts names the files by rank, not by id.",
                  "Energies are the rank-1 cluster means reported by caprieval."])


# ---------------------------------------------------------------- S2
def s2():
    rows = []
    for lig, pre in (("MMMA", "MMAREP_MMT_blind"), ("compound 48/80", "Y3B")):
        per = []
        for s in SEEDS:
            r = ml.rank1(pre + s)
            if not r:
                continue
            sc, vd = float(r["score"]), float(r["vdw"])
            el, ds = W * float(r["elec"]), float(r["desolv"])
            per.append((sc, vd, el, ds))
            rows.append([lig, pre + s, s.lstrip("_i") or "917",
                         f"{sc:.3f}", f"{vd:.3f}", f"{el:.3f}", f"{ds:.3f}",
                         f"{100*abs(el)/abs(sc):.1f}", f"{100*abs(vd)/abs(sc):.1f}",
                         f"{100*abs(ds)/abs(sc):.1f}",
                         f"{100*(abs(el)+abs(vd)+abs(ds))/abs(sc):.1f}"])
        med = [statistics.median(x) for x in zip(*per)]
        sc, vd, el, ds = med
        rows.append([lig, "(median of 5)", "", f"{sc:.3f}", f"{vd:.3f}", f"{el:.3f}", f"{ds:.3f}",
                     f"{100*abs(el)/abs(sc):.1f}", f"{100*abs(vd)/abs(sc):.1f}",
                     f"{100*abs(ds)/abs(sc):.1f}",
                     f"{100*(abs(el)+abs(vd)+abs(ds))/abs(sc):.1f}"])
    return write("S2_score_terms_wt.tsv",
                 ["ligand", "run", "iniseed", "score", "vdw", "elec_contribution_0.2x",
                  "desolv", "elec_pct", "vdw_pct", "desolv_pct", "sum_pct"],
                 rows,
                 ["Table S2. HADDOCK score partition against the wild-type receptor.",
                  "elec_contribution is 0.2 x elec, the emref weight. Percentages are of |score|.",
                  "The median row takes the median of each term separately, so sum_pct is not",
                  "exactly 100 (Figure 2A carries the same note)."])


# ---------------------------------------------------------------- S4
def s4():
    freq = {}
    for lab, pre in (("MMMA", "MMAREP_MMT_blind"), ("C4880", "Y3B")):
        for _, _, m in poses(pre):
            lig, prot = ml.load(m)
            for key in contact_set(lig, prot):
                freq.setdefault(key, {"MMMA": 0, "C4880": 0})[lab] += 1
    rows = []
    for (ri, rn), d in sorted(freq.items()):
        if d["MMMA"] and d["C4880"]:
            cat = "both"
        elif d["C4880"]:
            cat = "compound 48/80 only"
        else:
            cat = "MMMA only"
        rows.append([ri, rn, d["MMMA"], d["C4880"], cat])
    only = [r for r in rows if r[4] == "MMMA only"]
    print(f"    MMMA のみが接触する残基: {len(only)} 個 {only if only else '（本文の主張どおり 0）'}")
    return write("S4_contact_residues_wt.tsv",
                 ["residue_number", "residue_name", "poses_contacting_MMMA_of_20",
                  "poses_contacting_C4880_of_20", "category"],
                 rows,
                 [f"Table S4. Receptor residues within {CONTACT} A of a ligand heavy atom,",
                  "counted over the 20 poses of each ligand (rank-1 cluster, 4 models x 5 seeds).",
                  "A residue counts as contacted if any receptor atom is within the cutoff."])


# ---------------------------------------------------------------- S5
def s5():
    """184 への距離。T9 と同じ定義を主列に置き、補助の測り方を併記する。

    T9（本文・Figure 2B）の定義:
      陽イオン性窒素を持つリガンド（48/80）→ **窒素**から 184 の側鎖極性原子まで
      窒素を持たないリガンド（MMMA）      → **重原子**から同じ原子まで
    48/80 では炭素のほうが窒素より近いポーズがあるので、両方を出さないと表が読めない。
    """
    rows, by = [], {}
    for lig, rec, pre in LADDER:
        for seed, model, m in poses(pre):
            l_, prot = ml.load(m)
            rr = {a[1] for a in prot if a[0] == 184}
            resn = rr.pop() if rr else "?"
            tgt = [a for a in prot if a[0] == 184 and a[2] in BYRESN.get(resn, ())]
            heavy = [a for a in l_ if a[3] != "H"]
            run = pre + ("" if seed == "917" else "_i" + seed)
            if not tgt:
                rows.append([lig, rec, run, seed, model, resn,
                             "none (no polar side-chain atom)", "", "", "", "", "", ""])
                continue
            has_n = any(x[3] == "N" for x in heavy)
            measure = "ligand N" if has_n else "ligand heavy atom"
            prim = [a for a in heavy if a[3] == "N"] if has_n else heavy
            d1, a1 = min(((math.dist(a[4], t[4]), a) for a in prim for t in tgt),
                         key=lambda x: x[0])
            dh, ah = min(((math.dist(a[4], t[4]), a) for a in heavy for t in tgt),
                         key=lambda x: x[0])
            ox = [a for a in heavy if a[3] == "O"]
            do, ao = min(((math.dist(a[4], t[4]), a) for a in ox for t in tgt),
                         key=lambda x: x[0]) if ox else (None, None)
            rows.append([lig, rec, run, seed, model, resn, measure,
                         a1[2], f"{d1:.2f}", ah[2], ah[3], f"{dh:.2f}",
                         "" if do is None else f"{do:.2f}"])
            by.setdefault((lig, rec), []).append(d1)

    # T9 の dist_184_range_A と一致することを確かめる（定義がずれたら気づけるように）
    t9 = {(r["ligand"], r["receptor"]): r["dist_184_range_A"]
          for r in csv.DictReader(open(os.path.join(TAB, "T9_D184_ladder.tsv")), delimiter="\t")}
    bad = [k for k, v in by.items()
           if t9.get(k, "-") != f"{min(v):.2f}-{max(v):.2f}"]
    print(f"    T9 の dist_184_range_A と照合: {'全条件一致' if not bad else '不一致 ' + str(bad)}")

    return write("S5_distance_to_184.tsv",
                 ["ligand", "receptor", "run", "iniseed", "model", "residue184_name",
                  "primary_measure", "primary_ligand_atom", "primary_distance_A",
                  "nearest_heavy_atom", "nearest_heavy_element", "nearest_heavy_distance_A",
                  "nearest_oxygen_distance_A"],
                 rows,
                 ["Table S5. Distances from the ligand to the side-chain polar atoms of residue 184",
                  "(Asp OD1,OD2; Glu OE1,OE2; Asn OD1,ND2; Gln OE1,NE2; His ND1,NE2;",
                  "Ala and Gly have none).",
                  "primary_distance_A uses the measure quoted in the Letter and in Figure 2B:",
                  "the cationic nitrogen for compound 48/80, any heavy atom for MMMA (which has",
                  "no nitrogen). The nearest_heavy and nearest_oxygen columns give the same pose",
                  "measured the other two ways; for compound 48/80 a carbon is sometimes closer",
                  "than the nitrogen, so the two differ."])


# ---------------------------------------------------------------- S7
def s7():
    ref = {}
    for lig, rec, pre in LADDER:
        if rec == "WT":
            l_, prot = ml.load(poses(pre)[0][2])
            ref[lig] = contact_set(l_, prot)
    rows = []
    for lig, rec, pre in LADDER:
        for seed, model, m in poses(pre):
            l_, prot = ml.load(m)
            c = contact_set(l_, prot)
            j = len(c & ref[lig]) / len(c | ref[lig])
            rows.append([lig, rec, pre + ("" if seed == "917" else "_i" + seed), seed, model,
                         len(c), f"{j:.2f}", "yes" if j >= SAME_SITE else "no"])
    return write("S7_jaccard_per_pose.tsv",
                 ["ligand", "receptor", "run", "iniseed", "model",
                  "n_contact_residues", "jaccard_vs_wt", "same_site"],
                 rows,
                 ["Table S7. Per-pose overlap of the contact-residue set with the wild type.",
                  f"Reference = model 1 of the wild-type run at iniseed 917, per ligand.",
                  f"same_site is Jaccard >= {SAME_SITE}, the threshold fixed before the runs",
                  "(PREREG_TRIMER.md section 5, T4)."])


# ---------------------------------------------------------------- S9
def load_cif(p):
    rows, hdr, inloop = [], [], False
    for l in open(p):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
            inloop = True
            continue
        if inloop:
            if l.startswith("#") or not l.strip():
                break
            f = l.split()
            if len(f) != len(hdr):
                break
            rows.append(dict(zip(hdr, f)))
    return rows


def s9():
    def xyz(r):
        return (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
    rows = []
    for pdb, binder, pick in (
            ("7vv6", "compound 48/80 (het 6IB)",
             lambda r: r["label_comp_id"] == "6IB"),
            ("7vdm", "substance P (chain L, Arg1-Pro2-Lys3 modelled)",
             lambda r: r["auth_asym_id"] == "L")):
        cif = load_cif(os.path.join(CACHE, pdb + ".cif"))
        lig = [r for r in cif if pick(r) and r["type_symbol"] != "H"]
        rec = [r for r in cif if r["auth_asym_id"] == "R" and r["type_symbol"] != "H"
               and not pick(r)]
        L = [xyz(r) for r in lig]
        best = {}
        for r in rec:
            d = min(math.dist(xyz(r), p) for p in L)
            k = (int(r["auth_seq_id"]), r["label_comp_id"])
            if k not in best or d < best[k][0]:
                best[k] = (d, r["label_atom_id"])
        for res in (164, 168, 170, 184, 243, 248, 254):
            hit = [(k, v) for k, v in best.items() if k[0] == res]
            if not hit:
                rows.append([pdb.upper(), binder, res, "-", "-", "-"])
                continue
            (ri, rn), (d, an) = hit[0]
            rows.append([pdb.upper(), binder, f"{rn}{ri}", an, f"{d:.2f}",
                         "yes" if d <= CONTACT else "no"])
    return write("S9_cryoem_contacts.tsv",
                 ["pdb", "binder", "receptor_residue", "receptor_atom",
                  "min_distance_A", f"within_{CONTACT}A"],
                 rows,
                 ["Table S9. Distances in the deposited cryo-EM structures, measured the same",
                  "way as the docking poses (shortest heavy-atom to heavy-atom distance).",
                  "In 7VDM only residues 1-3 of substance P are modelled."])


def main():
    print("Supplementary Tables を書く")
    for fn in (s1, s2, s4, s5, s7, s9):
        fn()
    # S3 は T9 そのもの。SI の番号でも引けるようにコピーを置く
    import shutil
    src = os.path.join(TAB, "T9_D184_ladder.tsv")
    dst = os.path.join(TAB, "S3_D184_ladder.tsv")
    shutil.copyfile(src, dst)
    print(f"  {'S3_D184_ladder.tsv':28s} （T9_D184_ladder.tsv のコピー）")
    print("\nS6 は scripts/make_s6_hbond_geometry.py、S8 は scripts/make_s8_pose_rmsd.py で作る")


if __name__ == "__main__":
    main()
