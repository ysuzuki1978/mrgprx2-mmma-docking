#!/usr/bin/env python3
"""レター v2 が引用する D184 ラダーの数値台帳を作る。

  python3 scripts/make_ladder_table.py      # letter_MMA/ 直下で実行

出力: data/tables/T9_D184_ladder.tsv

T7（v1 用）との違い:
  * T7 は seed 917 の 1 本。T9 は **5 反復の中央値**
  * T7 の compound 48/80 は欠陥リガンド（二量体 + 遊離単量体）。
    T9 は**正しい共有結合三量体**（PDB chem comp 6IB の結合配座）
  * T9 は D184E と D184H も含む（ラダー全体）

T7 は v1 が何を主張したかの記録として残す。消さない。

数値の出所はすべて `data/runs/<run>/capri_clt.tsv` の rank-1 クラスタ
（`cluster_rank == 1` の行）と、同ディレクトリの `cluster_1_model_*.pdb.gz`。
"""
import glob
import gzip
import math
import os
import statistics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "data", "runs")
OUT = os.path.join(ROOT, "data", "tables", "T9_D184_ladder.tsv")
SEEDS = ["", "_i20000", "_i40000", "_i60000", "_i80000"]
W_EMREF = 0.2                      # emref の静電の重み
CONTACT = 4.0                      # 接触の閾値 (A)
SAME_SITE = 0.56                   # 同一部位とみなす Jaccard の下限（PREREG_TRIMER.md section 5 の T4）

# (リガンド, 受容体ラベル, run 接頭辞, 計算機)
SETS = [
    ("Compound4880", "WT",          "Y3B",          "yakuri3"),
    ("Compound4880", "D184N",       "Y3BD184N",     "yakuri3"),
    ("Compound4880", "D184A",       "Y3BD184A",     "yakuri3"),
    ("Compound4880", "D184E_rot1",  "Y3BD184E",     "yakuri3"),
    ("Compound4880", "D184E_rot2",  "Y3BD184E2",    "yakuri3"),
    ("Compound4880", "D184H_rot1P", "Y3BD184HP",    "yakuri3"),
    ("Compound4880", "D184H_rot1N", "Y3BD184HN",    "yakuri3"),
    ("Compound4880", "D184H_rot2P", "Y3BD184HP2",   "yakuri3"),
    ("MMT", "WT",          "MMAREP_MMT_blind",       "EPYC"),
    ("MMT", "D184N",       "MMAREP_MMT_D184N_blind", "EPYC"),
    ("MMT", "D184A",       "MMAREP_MMT_D184A_blind", "EPYC"),
    ("MMT", "E164A",       "MMAREP_MMT_E164A_blind", "EPYC"),
    ("MMT", "D184H_rot1P", "D184H_rot1P_MMT_blind",  "EPYC"),
    ("MMT", "D184H_rot1N", "D184H_rot1N_MMT_blind",  "EPYC"),
    ("MMT", "D184H_rot2P", "D184H_rot2P_MMT_blind",  "EPYC"),
    # 2026-10-23 追加（事前登録 section 17）。yakuri3 で取った MMMA。
    # EPYC 系列と分けておく（WT 対照も yakuri3 で取ったので、Δ は同一機械で算出される）。
    ("MMT_y3", "WT",         "Y3MWT",  "yakuri3"),
    ("MMT_y3", "D184E_rot1", "Y3ME",   "yakuri3"),
    ("MMT_y3", "D184E_rot2", "Y3ME2",  "yakuri3"),
]
# 残基種ごとの側鎖末端の極性原子（変異で残基種が変わるため、PDB の残基名から決める）
BYRESN = {"GLU": ("OE1", "OE2"), "ASP": ("OD1", "OD2"), "GLN": ("OE1", "NE2"),
          "ASN": ("OD1", "ND2"), "HIS": ("ND1", "NE2"), "ALA": (), "GLY": ()}


def rank1(run):
    p = os.path.join(RUNS, run, "capri_clt.tsv")
    if not os.path.exists(p):
        return None
    lines = [l for l in open(p) if not l.startswith("#")]
    head = lines[0].rstrip("\n").split("\t")
    for l in lines[1:]:
        d = dict(zip(head, l.rstrip("\n").split("\t")))
        if d.get("cluster_rank") == "1":
            return d
    return None


def models(run):
    return sorted(glob.glob(os.path.join(RUNS, run, "cluster_1_model_*.pdb.gz")))


def load(path):
    lig, prot = [], []
    with gzip.open(path, "rt") as f:
        for l in f:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            an, rn, ch = l[12:16].strip(), l[17:20].strip(), l[21]
            ri = int(l[22:26])
            el = (l[76:78].strip() or an[0]).upper()
            x = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            (lig if (l.startswith("HETATM") or ch == "B") else prot).append((ri, rn, an, el, x))
    return lig, prot


def geom(run_prefix, ref_contacts=None):
    """20 モデル分の「184 への最短距離」「WT 接触集合との Jaccard」を集める。"""
    d184, jac, csets = [], [], []   # d184 は 48/80 では N、MMMA では重原子からの距離
    for s in SEEDS:
        for m in models(run_prefix + s):
            lig, prot = load(m)
            ns = [a[4] for a in lig if a[3] == "N"]
            heavy = [a[4] for a in lig if a[3] != "H"]
            # 184 への最短距離（リガンド窒素 -> 側鎖末端の極性原子）
            rr = {a[1] for a in prot if a[0] == 184}
            resn = rr.pop() if rr else "?"
            names = BYRESN.get(resn, ())
            tgt = [a[4] for a in prot if a[0] == 184 and a[2] in names]
            if tgt:
                if ns:          # 陽イオン性窒素を持つリガンド（48/80）
                    d184.append(min(math.dist(n, t) for n in ns for t in tgt))
                else:           # 窒素を持たないリガンド（MMMA）は重原子で測る
                    d184.append(min(math.dist(h, t) for h in heavy for t in tgt))
            # 接触残基集合
            c = {(a[0], a[1]) for a in prot
                 if min(math.dist(a[4], h) for h in heavy) <= CONTACT}
            csets.append(c)
            if ref_contacts is not None:
                jac.append(len(c & ref_contacts) / len(c | ref_contacts))
    return d184, jac, csets


def f3(x):
    return "" if x is None else f"{x:.3f}"


def main():
    med = statistics.median
    rows, base = [], {}
    # 各リガンドの WT 接触集合（seed 917 の model 1）を Jaccard の基準にする
    ref = {}
    for lig, lab, pre, _ in SETS:
        if lab == "WT":
            _, _, cs = geom(pre)
            ref[lig] = cs[0]

    for lig, lab, pre, cpu in SETS:
        rs = [r for r in (rank1(pre + s) for s in SEEDS) if r]
        if not rs:
            print(f"  ! {pre} が無い"); continue
        sc = [float(r["score"]) for r in rs]
        vd = [float(r["vdw"]) for r in rs]
        el = [float(r["elec"]) for r in rs]
        ds = [float(r["desolv"]) for r in rs]
        ec = [W_EMREF * abs(v) for v in el]          # 電荷の寄与（正の値）
        m_sc, m_ec, m_vd = med(sc), med(ec), med(vd)
        d184, jac, _ = geom(pre, ref[lig])
        n_contact = sum(1 for d in d184 if d < CONTACT)
        if lab == "WT":
            base[lig] = (m_sc, m_ec)
            d_sc = d_ec = pct_sc = pct_ec = None
        else:
            b_sc, b_ec = base[lig]
            d_sc, d_ec = m_sc - b_sc, m_ec - b_ec
            pct_sc = 100 * (abs(m_sc) - abs(b_sc)) / abs(b_sc)
            pct_ec = 100 * (m_ec - b_ec) / b_ec
        rows.append([
            lig, lab, cpu, str(len(rs)),
            f3(m_sc), f3(m_vd), f3(med(el)), f3(med(ds)), f3(m_ec),
            f"{100*m_ec/abs(m_sc):.1f}", f"{100*abs(m_vd)/abs(m_sc):.1f}",
            f3(d_sc), f3(d_ec),
            "" if pct_sc is None else f"{pct_sc:.1f}",
            "" if pct_ec is None else f"{pct_ec:.1f}",
            f"{n_contact}/{len(d184)}" if d184 else "-",
            f"{min(d184):.2f}-{max(d184):.2f}" if d184 else "-",
            f"{min(jac):.2f}-{max(jac):.2f}" if jac else "-",
            f"{sum(1 for v in jac if v >= SAME_SITE)}/{len(jac)}" if jac else "-",
            f"{max(sc)-min(sc):.3f}",
        ])

    head = ["ligand", "receptor", "cpu", "n_replicates",
            "median_score", "median_vdw", "median_elec", "median_desolv",
            "elec_contribution", "elec_pct_of_score", "vdw_pct_of_score",
            "d_score", "d_elec_contribution",
            "d_score_pct_of_WT", "d_elec_contribution_pct_of_WT",
            "contacts_184_within_4A", "dist_184_range_A", "jaccard_vs_WT_range",
            "site_kept_of_20",
            "score_range"]
    with open(OUT, "w") as f:
        f.write("\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"{OUT} に {len(rows)} 行を書いた")
    for r in rows:
        print(f"  {r[0]:13s}{r[1]:12s}n={r[3]}  score={r[4]:>9s}  "
              f"elec寄与={r[8]:>7s} ({r[9]:>5s}%)  Δ%={r[13]:>6s}  184接触={r[15]:>6s}")


if __name__ == "__main__":
    main()
