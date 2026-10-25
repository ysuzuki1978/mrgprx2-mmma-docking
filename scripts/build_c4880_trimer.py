#!/usr/bin/env python3
"""compound 48/80 の正しい共有結合三量体を SMILES から構築する（遺伝研で実行）。

  <acpype env python> build_c4880_trimer.py <出力ディレクトリ>

事前登録: letter_MMA/PREREG_TRIMER.md §2（投入前に固定）

既存の `Compound4880_chainB.pdb` は共有結合した三量体ではなく、メチレン橋 1 本を欠く
**二量体 + 単量体の 2 分子**だった（2026-10-06 判明、CRYOEM_CHECK.md §1.1）。
本スクリプトはそれを PDB 化学成分 6IB の SMILES から作り直す。

**実験構造 7VV6 の結合配座は入力に使わない。** 結合配座から始めて実験の姿勢を再現しても
循環論法になるため、SMILES から独立に 3D 配座を起こす。

配座生成は RDKit の ETKDGv3 を**固定シード**で使う。Open Babel の `--gen3d` は
(1) 速度引数を付けると座標を生成せず、(2) 引数なしだと多スレッドで 3.3 CPU 時間でも
終わらなかったため使わない（2026-10-06 に両方を実測）。

PREREG_TRIMER.md §2.2 の検査をすべてここで実施する。1 つでも落ちたら異常終了する。
"""
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# PDB 化学成分 6IB（RCSB 化学成分 API、2026-10-06 取得）
SMILES_NEUTRAL = "CNCCc1ccc(OC)c(Cc2cc(CCNC)cc(Cc3cc(CCNC)ccc3OC)c2OC)c1"
# 上の 3 つの第二級アミンを pH 7.4 で帯電させた形（正味 +3）。
# 置換箇所は CNCC が 1 つ、CCNC が 2 つ。目視でなく下の検査で確認する。
SMILES_CATION = ("C[NH2+]CCc1ccc(OC)c(Cc2cc(CC[NH2+]C)cc"
                 "(Cc3cc(CC[NH2+]C)ccc3OC)c2OC)c1")
INCHIKEY_SKELETON = "DXEVTBAWIGXNAW"      # 中性形の InChIKey 前半 14 文字
SEED = 917                                 # 既存 run の iniseed と同じ値を使う
EXPECT_HEAVY = {"C": 32, "N": 3, "O": 3}


def fail(msg):
    sys.exit(f"検査不合格: {msg} → 投入しない")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)

    # 1. 中性形の同一性を InChIKey で確認する
    neu = Chem.MolFromSmiles(SMILES_NEUTRAL)
    assert neu is not None, "中性 SMILES を解釈できない"
    key_neu = Chem.MolToInchiKey(neu)
    print(f"[1] 中性形 InChIKey = {key_neu}")
    if not key_neu.startswith(INCHIKEY_SKELETON):
        fail(f"中性形の InChIKey が {key_neu}（{INCHIKEY_SKELETON} で始まるべき）")

    # 2. カチオン形が、中性形と同じ骨格であることを確認する
    cat = Chem.MolFromSmiles(SMILES_CATION)
    assert cat is not None, "カチオン SMILES を解釈できない"
    charge = Chem.GetFormalCharge(cat)
    print(f"[2] カチオン形の正味電荷 = {charge:+d}、"
          f"分子式 = {rdMolDescriptors.CalcMolFormula(cat)}")
    if charge != 3:
        fail(f"正味電荷が {charge:+d}（+3 であるべき）")
    # 電荷と付加プロトンを外して中性形に戻し、InChIKey が一致するかを見る
    back = Chem.RWMol(cat)
    for a in back.GetAtoms():
        if a.GetFormalCharge() != 0:
            a.SetFormalCharge(0)
            a.SetNumExplicitHs(0)
            a.SetNoImplicit(False)
    back = back.GetMol()
    Chem.SanitizeMol(back)
    key_back = Chem.MolToInchiKey(back)
    if key_back != key_neu:
        fail(f"カチオン形を中性に戻すと InChIKey が {key_back}（{key_neu} であるべき）")
    print(f"[3] カチオン形を中性に戻すと InChIKey が一致した")

    # 3. 3D 配座を生成（ETKDGv3、固定シード）して MMFF94 で最適化する
    m = Chem.AddHs(cat)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = SEED
    if AllChem.EmbedMolecule(m, ps) != 0:
        fail("ETKDGv3 で 3D 配座を生成できなかった")
    conv = AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    print(f"[4] ETKDGv3 (seed={SEED}) で配座を生成、MMFF94 で最適化"
          f"（収束コード {conv}: 0 なら収束）")
    if conv not in (0, 1):
        fail(f"MMFF94 最適化の戻り値が {conv}")

    # 4. 検査（PREREG_TRIMER.md §2.2）
    heavy = [a for a in m.GetAtoms() if a.GetSymbol() != "H"]
    comp = {}
    for a in heavy:
        comp[a.GetSymbol()] = comp.get(a.GetSymbol(), 0) + 1
    frags = Chem.GetMolFrags(m)
    rings = rdMolDescriptors.CalcNumRings(m)
    arom = rdMolDescriptors.CalcNumAromaticRings(m)
    # メトキシ（炭素 2 個に挟まれた O）と第二級アンモニウム（炭素 2 個に挟まれた N）
    o2c = sum(1 for a in heavy if a.GetSymbol() == "O"
              and sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "C") == 2)
    n2c = sum(1 for a in heavy if a.GetSymbol() == "N"
              and sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "C") == 2)
    nh = m.GetNumAtoms() - len(heavy)

    print(f"[5] 重原子 {len(heavy)}、元素 {comp}、水素 {nh}")
    print(f"[6] 連結成分 {len(frags)}、環 {rings}（うち芳香環 {arom}）")
    print(f"[7] メトキシ型の O {o2c} 個、第二級アミン型の N {n2c} 個")

    checks = [("1 重原子数 38", len(heavy) == 38),
              ("2 元素構成 C32/N3/O3", comp == EXPECT_HEAVY),
              ("3 連結成分 1", len(frags) == 1),
              ("4 環 3（すべて芳香環）", rings == 3 and arom == 3),
              ("5 メトキシ 3", o2c == 3),
              ("6 第二級アミン 3", n2c == 3),
              ("7 正味電荷 +3", Chem.GetFormalCharge(m) == 3),
              ("8 InChIKey 一致", True)]
    bad = [n for n, ok in checks if not ok]
    for n, ok in checks:
        print(f"    {'OK ' if ok else 'NG '} 検査 {n}")
    if bad:
        fail(" / ".join(bad))

    # 5. 結合長の妥当性（配座生成が壊れていないことの確認）
    confm = m.GetConformer()
    import math
    lens = []
    for b in m.GetBonds():
        i, j = b.GetBeginAtom(), b.GetEndAtom()
        if i.GetSymbol() == "H" or j.GetSymbol() == "H":
            continue
        pi, pj = confm.GetAtomPosition(i.GetIdx()), confm.GetAtomPosition(j.GetIdx())
        lens.append(math.dist((pi.x, pi.y, pi.z), (pj.x, pj.y, pj.z)))
    print(f"[8] 重原子間の結合長 {min(lens):.3f}–{max(lens):.3f} A "
          f"(平均 {sum(lens)/len(lens):.3f})")
    if not (1.2 <= min(lens) and max(lens) <= 1.7):
        fail(f"結合長が 1.2–1.7 A の範囲外（{min(lens):.3f}–{max(lens):.3f}）")

    # 6. PDB で書き出す（残基名 6IB、chain/segid B）
    path = os.path.join(out, "C4880_trimer_chainB.pdb")
    info_done = set()
    for idx, a in enumerate(m.GetAtoms(), start=1):
        nm = f"{a.GetSymbol()}{idx}"
        ri = Chem.AtomPDBResidueInfo()
        ri.SetResidueName("6IB")
        ri.SetResidueNumber(1)
        ri.SetChainId("B")
        ri.SetName(f"{nm:<4s}"[:4])
        ri.SetIsHeteroAtom(True)
        a.SetMonomerInfo(ri)
        info_done.add(nm)
    assert len(info_done) == m.GetNumAtoms(), "原子名が重複した"
    Chem.MolToPDBFile(m, path, flavor=4)

    # segid 欄 (73-76) を B にし、CONECT を落とす（haddock3 は使わない）
    lines = [f"REMARK   1 compound 48/80 trimer. PDB chem comp 6IB SMILES, RDKit ETKDGv3\n",
             f"REMARK   1 seed {SEED}, MMFF94 optimised, 3 amines protonated (net +3).\n",
             f"REMARK   1 Pre-registered: letter_MMA/PREREG_TRIMER.md\n"]
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            l = l.rstrip("\n").ljust(80)
            lines.append("HETATM" + l[6:72] + "B   " + l[76:] + "\n")
    lines.append("END\n")
    with open(path, "w") as f:
        f.writelines(lines)
    n_at = sum(1 for l in open(path) if l.startswith("HETATM"))
    print(f"\n[出力] {path}  HETATM {n_at} 行")
    assert n_at == m.GetNumAtoms(), "書き出した原子数が合わない"

    # acpype に渡す正味電荷を記録する
    with open(os.path.join(out, "NET_CHARGE"), "w") as f:
        f.write("3\n")
    print(f"[出力] {os.path.join(out, 'NET_CHARGE')}  (acpype の -n に渡す)")


if __name__ == "__main__":
    main()
