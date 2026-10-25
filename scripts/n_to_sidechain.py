#!/usr/bin/env python3
"""リガンドの窒素から、指定残基の側鎖末端の極性原子までの最短距離を測る。

  python3 scripts/n_to_sidechain.py Y3B Y3BD184N ...

rank-1 クラスタの 4 モデル（`cluster_1_model_*.pdb.gz`）それぞれについて
164 / 184 / 254 番までの距離を出す。変異体では残基種が変わるので、
PDB に実際に入っている残基名から対象原子を決める（固定表を引かない）。

numpy を使わない（素の python3 で動く）。
"""
import glob
import gzip
import math
import os
import sys

# 残基種ごとの「側鎖末端の極性原子」。Ala/Gly は側鎖に極性原子が無い。
BYRESN = {"GLU": ("OE1", "OE2"), "ASP": ("OD1", "OD2"),
          "GLN": ("OE1", "NE2"), "ASN": ("OD1", "ND2"),
          "HIS": ("ND1", "NE2"),
          "ALA": (), "GLY": ()}
TARGETS = (164, 184, 254)


def load(path):
    """(リガンド原子, 受容体原子) を返す。リガンドは HETATM か chain B。"""
    lig, prot = [], []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for l in f:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            an = l[12:16].strip()
            rn = l[17:20].strip()
            ch = l[21]
            ri = int(l[22:26])
            el = (l[76:78].strip() or an[0]).upper()
            x = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            if l.startswith("HETATM") or ch == "B":
                lig.append((an, el, x))
            else:
                prot.append((ri, rn, an, x))
    return lig, prot


def mindist(lig, prot, resi):
    ns = [x for _, el, x in lig if el == "N"]
    rr = {r for ri, r, _, _ in prot if ri == resi}
    resn = rr.pop() if rr else "?"
    if not ns or resn not in BYRESN:
        return None, resn
    names = BYRESN[resn]
    tgt = [x for ri, _, a, x in prot if ri == resi and a in names]
    if not tgt:
        return None, resn                 # Ala/Gly など
    return min(math.dist(n, t) for n in ns for t in tgt), resn


def main():
    runs = sys.argv[1:]
    print(f"{'run':22s}{'model':7s}" + "".join(f"{'N..'+str(r):>9s}" for r in TARGETS)
          + "   (残基名)")
    for run in runs:
        for m in sorted(glob.glob(f"data/runs/{run}/cluster_1_model_*.pdb.gz")):
            lig, prot = load(m)
            cols, names = [], []
            for resi in TARGETS:
                dd, rn = mindist(lig, prot, resi)
                cols.append("  -  " if dd is None else f"{dd:5.2f}")
                names.append(rn)
            tag = os.path.basename(m).split("_")[-1][:-7]
            print(f"{run:22s}{tag:7s}" + "".join(f"{c:>9s}" for c in cols)
                  + "   " + "/".join(names))


if __name__ == "__main__":
    main()
