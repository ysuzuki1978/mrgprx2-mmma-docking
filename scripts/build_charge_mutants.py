#!/usr/bin/env python3
"""カルボキシラート 3 部位の電荷を消した受容体モデルを作る（2^3-1 = 7 通り）。

  python3 build_charge_mutants.py <WT.pdb> <出力ディレクトリ>

対象は MRGPRX2 のオルソステリック部位を囲む 3 つのカルボキシラート:

  Glu164 → Gln  (OE2 を NE2 に改名)
  Asp184 → Asn  (OD2 を ND2 に改名)
  Asp254 → Asn  (OD2 を ND2 に改名)

**原子は 1 つも動かさない。** 酸素を窒素に改名してアミドに変えるだけで、負電荷だけが
消える。レターの既存変異体 D184N がこの方法で作られており（2026-09-28 に現物で確認、
座標移動 0 原子）、本スクリプトはそれを 3 部位に一般化したもの。

側鎖を作り直さないので配座の自由度がなく、D184H で問題になった配座依存
（PREREG_D184H.md §2）が原理的に起きない。

出力は重原子のみ。水素は haddock3 の topoaa が付ける。
"""
import itertools
import os
import sys

# 部位 -> (元の残基名, 新しい残基名, 改名する原子, 改名後の原子名)
SITES = {
    "E164Q": (164, "GLU", "GLN", "OE2", "NE2"),
    "D184N": (184, "ASP", "ASN", "OD2", "ND2"),
    "D254N": (254, "ASP", "ASN", "OD2", "ND2"),
}
ORDER = ["E164Q", "D184N", "D254N"]


def mutate(lines, muts):
    """指定した部位だけ残基名と原子名を書き換える。座標は触らない。"""
    out, hits = [], {m: 0 for m in muts}
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            resi = int(l[22:26])
            name = l[12:16].strip()
            for m in muts:
                site, old_res, new_res, old_at, new_at = SITES[m]
                if resi != site:
                    continue
                assert l[17:20].strip() == old_res, f"{site} 番が {old_res} でない"
                l = l[:17] + f"{new_res:<3s}" + l[20:]
                if name == old_at:
                    l = l[:12] + f" {new_at:<3s}" + l[16:]
                    # 元素欄（77-78 桁）と形式電荷欄（79-80 桁）も直す。
                    # WT の OD2/OE2 は "O1-"（カルボキシラートの負電荷）なので、
                    # 窒素に変えるときに電荷の表記を残してはいけない。
                    l = l.rstrip("\n").ljust(80)
                    l = l[:76] + f"{new_at[0]:>2s}" + "  " + "\n"
                    hits[m] += 1
        out.append(l)
    for m, n in hits.items():
        assert n == 1, f"{m}: 改名した原子が {n} 個（1 個であるべき）"
    return out


def main():
    wt_path = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else "."
    os.makedirs(outdir, exist_ok=True)
    lines = open(wt_path).read().splitlines(keepends=True)

    made = []
    for n in (1, 2, 3):
        for combo in itertools.combinations(ORDER, n):
            tag = "_".join(combo)
            out = os.path.join(outdir, f"MRGPRX2_{tag}_chainA.pdb")
            open(out, "w").writelines(mutate(lines, combo))
            made.append((tag, out))

    # 検証: 座標が 1 つも動いていないこと、原子数が変わっていないこと
    def coords(path):
        d = {}
        for l in open(path):
            if l.startswith(("ATOM", "HETATM")):
                d[(int(l[22:26]), l[30:54])] = l[12:16].strip()
        return d

    wt = coords(wt_path)
    print(f"{'モデル':24s}{'原子数':>7s}{'座標が動いた原子':>18s}")
    for tag, out in made:
        mu = coords(out)
        moved = len(set(wt) - set(mu))        # 座標キーが一致しない = 動いた
        print(f"{tag:24s}{len(mu):7d}{moved:18d}")
        assert len(mu) == len(wt), f"{tag}: 原子数が変わった"
        assert moved == 0, f"{tag}: 座標が動いた"
    print(f"\n{len(made)} モデルを {outdir} に作成した（すべて座標移動 0）")


if __name__ == "__main__":
    main()
