#!/usr/bin/env python3
"""D184 ラダー（A / E rot1 / E rot2 / H protonated / H neutral）の cfg を生成する。

  ssh yakuri3 'python3 ~/trimer_check/gen_d184_ladder.py'

事前登録: letter_MMA/PREREG_TRIMER.md §14（2026-10-08、投入前に固定）

Y3BE164Q.cfg（正しい三量体 × 変異受容体）をそのまま土台にする。変えるのは 3 つだけ:
  1. run_dir
  2. 受容体 PDB
  3. 184 番のプロトン化（H の中性アームのみ。autohis を切って HISE を明示）

HISE の指定は gen_d184h.py（2026-09-28、NIG）と同じ:
  P -> 144 HISE / 184 HIS(+1) / 259,261 HIS(+1)
  N -> 144 HISE / 184 HISE    / 259,261 HIS(+1)
天然の His は WT run と同一で、違いは 184 番だけ。
"""
import os
import re
import sys

B = "/home/yakuri/trimer_check"
CFG = f"{B}/cfg"
# seed ごとに、対応する E164Q の cfg をそのまま土台にする。
# iniseed は既定の 917 のときだけ cfg に書かれないので、土台を seed 別に取るのが安全。
def base_cfg(seed):
    return f"{CFG}/Y3BE164Q{seed}.cfg"
E164Q_RECEPTOR = f"{B}/inputs/MRGPRX2_E164Q_chainA.pdb"
SEEDS = ["", "_i20000", "_i40000", "_i60000", "_i80000"]

# (run 接頭辞, 受容体, アーム)   アーム N = 184 を中性 HISE にする
CONDITIONS = [
    ("Y3BD184A",  f"{B}/inputs/MRGPRX2_D184A_chainA.pdb",      "P"),
    ("Y3BD184E",  f"{B}/inputs/MRGPRX2_D184E_chainA.pdb",      "P"),
    ("Y3BD184E2", f"{B}/inputs/MRGPRX2_D184E_chainA_rot2.pdb", "P"),
    ("Y3BD184HP", f"{B}/inputs/MRGPRX2_D184H_chainA.pdb",      "P"),
    ("Y3BD184HN", f"{B}/inputs/MRGPRX2_D184H_chainA.pdb",      "N"),
    # 2026-10-21 追加（事前登録 PREREG_TRIMER.md section 16）。rot2 の配座依存を確かめる。
    ("Y3BD184HP2", f"{B}/inputs/MRGPRX2_D184H_chainA_rot2.pdb", "P"),
]

HISE_BLOCK = """[topoaa.mol1]
# 184 番だけ中性（NE2 にプロトン）にする。144 は WT run の autohis と同じ HISE。
# 259 と 261 は指定しない＝HIS(+1) のままで、これも WT run と同じ。
nhise = 2
hise_1 = 144
hise_2 = 184

"""

HEADER = """# D184 ラダー x 正しい共有結合三量体。
# 事前登録: letter_MMA/PREREG_TRIMER.md section 14（2026-10-08、投入前）。
# 受容体と 184 のプロトン化以外は B 群（Y3B）と完全同一。
"""


def build(name, receptor, arm, seed):
    text = open(base_cfg(seed)).read()

    run = name + seed
    text, n = re.subn(r'run_dir\s*=\s*"[^"]*"', f'run_dir = "{B}/runs/{run}/run1"', text)
    assert n == 1, f"{run}: run_dir が {n} 箇所"

    assert text.count(E164Q_RECEPTOR) == 1, f"{run}: 受容体のパスが 1 箇所でない"
    text = text.replace(E164Q_RECEPTOR, receptor)

    # iniseed は土台の cfg が持っている（[rigidbody] の中。既定 917 のときは行が無い）。
    if seed:
        want = int(seed.lstrip("_i"))
        assert re.search(rf"^iniseed\s*=\s*{want}\b", text, flags=re.M), \
            f"{run}: iniseed {want} が土台に無い"
    else:
        assert not re.search(r"^iniseed\s*=", text, flags=re.M), \
            f"{run}: 既定 seed の土台に iniseed 行がある"

    if arm == "N":
        text, n = re.subn(r"^autohis\s*=\s*true", "autohis = false", text, flags=re.M)
        assert n == 1, f"{run}: autohis 行が {n} 箇所"
        lines = text.splitlines(keepends=True)
        i = next(k for k, l in enumerate(lines) if l.startswith("[topoaa]"))
        j = next(k for k in range(i + 1, len(lines)) if lines[k].startswith("["))
        text = "".join(lines[:j]) + HISE_BLOCK + "".join(lines[j:])
    else:
        assert re.search(r"^autohis\s*=\s*true", text, flags=re.M), f"{run}: autohis=true が無い"

    return run, HEADER + text


def main():
    # 引数を与えると、その接頭辞の条件だけを生成する（既存 cfg を触らない）。
    only = sys.argv[1] if len(sys.argv) > 1 else None
    runs = []
    for name, receptor, arm in CONDITIONS:
        if only and name != only:
            continue
        for seed in SEEDS:
            run, text = build(name, receptor, arm, seed)
            with open(f"{CFG}/{run}.cfg", "w") as f:
                f.write(text)
            runs.append(run)
    if only:
        rl = f"{B}/runlist_{only}.txt"
        with open(rl, "w") as f:
            f.write("\n".join(runs) + "\n")
        print(f"cfg {len(runs)} 本を書いた -> {rl}")
        for r in runs:
            print(f"  {r}")
        return
    # 2 本並走させるので runlist を 2 つに割る。条件が偏らないよう交互に配る。
    for k in (0, 1):
        with open(f"{B}/runlist_ladder{k + 1}.txt", "w") as f:
            f.write("\n".join(runs[k::2]) + "\n")
    print(f"cfg {len(runs)} 本を書いた")
    for k in (0, 1):
        got = runs[k::2]
        print(f"  runlist_ladder{k + 1}.txt : {len(got)} 本  {got[0]} ... {got[-1]}")


if __name__ == "__main__":
    main()
