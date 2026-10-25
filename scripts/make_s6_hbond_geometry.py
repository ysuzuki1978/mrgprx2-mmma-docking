#!/usr/bin/env python3
"""Table S6: Asp184 のカルボキシラートに向く水素結合の幾何を、ポーズごとに測る。

  python3 scripts/make_s6_hbond_geometry.py      # letter_MMA/ 直下で実行

出力: data/tables/S6_hbond_geometry.tsv（40 行 = 両リガンド × 20 ポーズ）

なぜ要るか
----------
レターの訂正の核心は「MMMA と D184 の接触は水素結合ではなく**弱い C–H⋯O** であり、
向きが通常の水素結合と逆」という点である。距離だけでは C–H⋯O と呼べるかを示せないので、
**供与体–H⋯受容体の角度**まで出す。比較のため compound 48/80 の N–H⋯O（塩橋）も同じ基準で測る。

判定の目安（Desiraju & Steiner, *The Weak Hydrogen Bond*, IUCr Monographs 9, 1999）:
  弱い C–H⋯O … H⋯O 2.0–3.0 Å、C–H⋯O 角 110–180 度

測り方
------
* ポーズ = rank-1 クラスタの 4 モデル × 5 反復 seed = 20（本文の n の定義と同じ）
* 供与体の割り当て: 各水素について、同じリガンド内で最も近い重原子をその水素の結合先とする
  （リガンドに CONECT が無いため。C–H も N–H も 1.0–1.1 Å なので一意に決まる）
* 受容体 = Asp184 の OD1 / OD2
* 各ポーズについて H⋯O が最短になる組を採る。48/80 は N–H⋯O に限った最短も併記する
  （最短が C–H になるポーズがあっても、塩橋そのものを追えるようにするため）
"""
import glob
import gzip
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "data", "runs")
OUT = os.path.join(ROOT, "data", "tables", "S6_hbond_geometry.tsv")
SEEDS = ["", "_i20000", "_i40000", "_i60000", "_i80000"]
ACCEPTOR = ("OD1", "OD2")          # Asp184 のカルボキシラート酸素
RESI = 184

SETS = [("MMMA", "MMAREP_MMT_blind"),
        ("compound 48/80", "Y3B")]


def load(path):
    lig, acc = [], []
    with gzip.open(path, "rt") as f:
        for l in f:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            an, ch, ri = l[12:16].strip(), l[21], int(l[22:26])
            el = (l[76:78].strip() or an[0]).upper()
            x = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            if ch == "B":
                lig.append((an, el, x))
            elif ri == RESI and an in ACCEPTOR:
                acc.append((an, el, x))
    return lig, acc


def angle(a, b, c):
    """a-b-c の b における角度（度）。"""
    u = [a[i] - b[i] for i in range(3)]
    v = [c[i] - b[i] for i in range(3)]
    nu = math.sqrt(sum(t * t for t in u))
    nv = math.sqrt(sum(t * t for t in v))
    cos = sum(u[i] * v[i] for i in range(3)) / (nu * nv)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def geometry(lig, acc, donor_elements=None):
    """H⋯O が最短になる (供与体, H, 受容体) を返す。"""
    heavy = [a for a in lig if a[1] != "H"]
    best = None
    for h in (a for a in lig if a[1] == "H"):
        d = min(heavy, key=lambda t: math.dist(h[2], t[2]))      # H の結合先
        if donor_elements and d[1] not in donor_elements:
            continue
        for o in acc:
            dho = math.dist(h[2], o[2])
            if best is None or dho < best[0]:
                best = (dho, d, h, o, angle(d[2], h[2], o[2]), math.dist(d[2], o[2]))
    return best


def main():
    rows, summary = [], {}
    for lab, pre in SETS:
        got = []
        for s in SEEDS:
            for m in sorted(glob.glob(os.path.join(RUNS, pre + s, "cluster_1_model_*.pdb.gz"))):
                lig, acc = load(m)
                if not acc:
                    continue
                model = os.path.basename(m).split("_model_")[1][0]
                b = geometry(lig, acc)
                # 48/80 は N–H⋯O（塩橋そのもの）も別に出す
                n = geometry(lig, acc, donor_elements={"N"}) if lab != "MMMA" else None
                rows.append([lab, pre + s, s.lstrip("_i") or "917", model,
                             f"{b[1][0]}({b[1][1]})", b[2][0], b[3][0],
                             f"{b[0]:.2f}", f"{b[4]:.0f}", f"{b[5]:.2f}",
                             "" if n is None else f"{n[1][0]}({n[1][1]})",
                             "" if n is None else f"{n[0]:.2f}",
                             "" if n is None else f"{n[4]:.0f}",
                             "" if n is None else f"{n[5]:.2f}"])
                got.append((b[0], b[4], b[5], b[1][1], n))
        summary[lab] = got
        # リガンドに O–H / N–H（＝水素結合供与体）があるかを確かめる
        lig, _ = load(sorted(glob.glob(os.path.join(RUNS, pre, "cluster_1_model_1.pdb.gz")))[0])
        heavy = [a for a in lig if a[1] != "H"]
        donors = sorted({min(heavy, key=lambda t: math.dist(h[2], t[2]))[1]
                         for h in lig if h[1] == "H"})
        summary[lab + "/donor_elements"] = donors

    head = ["ligand", "run", "iniseed", "model",
            "donor", "H", "acceptor", "d_H_O_A", "angle_D_H_O_deg", "d_D_O_A",
            "NH_donor", "NH_d_H_O_A", "NH_angle_deg", "NH_d_N_O_A"]
    with open(OUT, "w") as f:
        f.write("# Table S6. Hydrogen-bond geometry at Asp184 (20 poses per ligand).\n")
        f.write("# Poses = the 4 models of the rank-1 cluster x 5 replicate seeds.\n")
        f.write("# Donor assignment: nearest ligand heavy atom to each hydrogen.\n")
        f.write("# Weak C-H...O criteria (Desiraju & Steiner 1999): H...O 2.0-3.0 A, angle 110-180 deg.\n")
        f.write("\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"{OUT} に {len(rows)} 行を書いた\n")

    for lab, _ in SETS:
        g = summary[lab]
        dho = [x[0] for x in g]
        ang = [x[1] for x in g]
        ddo = [x[2] for x in g]
        els = sorted({x[3] for x in g})
        print(f"{lab}（n={len(g)} ポーズ）")
        print(f"  水素の結合先に現れる元素 : {', '.join(summary[lab + '/donor_elements'])}"
              f"  -> 水素結合供与体 (O–H / N–H) は"
              f"{'ない' if set(summary[lab + '/donor_elements']) <= {'C'} else 'ある'}")
        print(f"  最短の H⋯O の供与体      : {', '.join(els)}")
        print(f"  H⋯O                      : {min(dho):.2f}–{max(dho):.2f} Å"
              f"  （中央値 {sorted(dho)[len(dho)//2]:.2f}）")
        print(f"  供与体–H⋯O の角度        : {min(ang):.0f}–{max(ang):.0f}°"
              f"  （中央値 {sorted(ang)[len(ang)//2]:.0f}）")
        print(f"  供与体⋯O                 : {min(ddo):.2f}–{max(ddo):.2f} Å")
        if g[0][4] is not None:
            nd = [x[4][0] for x in g]
            na = [x[4][4] for x in g]
            nn = [x[4][5] for x in g]
            print(f"  N–H⋯O に限ると           : H⋯O {min(nd):.2f}–{max(nd):.2f} Å、"
                  f"角 {min(na):.0f}–{max(na):.0f}°、N⋯O {min(nn):.2f}–{max(nn):.2f} Å")
        print()


if __name__ == "__main__":
    main()
