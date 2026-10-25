#!/usr/bin/env python3
"""D184H 受容体モデルを「他の原子を一切動かさずに」作る。

  python3 scripts/build_d184h_rotamers.py <WT.pdb> <出力の接頭辞>

なぜ専用スクリプトかというと、レターの既存変異体が**原子の書き換えだけ**で作られて
いるため（2026-09-28 に実機確認）:

  D184N : OD2 を ND2 に改名しただけ。座標移動ゼロ
  D184A : CG/OD1/OD2 を削除しただけ。座標移動ゼロ
  E164A : CD/CG/OE1/OE2 を削除しただけ。座標移動ゼロ

His は Asp より大きいのでこの方法が使えない。側鎖パッカー（EvoEF2 BuildMutant）を
試したところ **TRP248 の側鎖を最大 4.44 A 動かした**。W248 はレターが名指しする
芳香族サブポケットの構成残基なので、それを動かしたモデルでは比較にならない。

そこで本スクリプトは:
  1. 184 番の N/CA/C/O/CB を WT のまま残す（Asp と His で共通の原子）
  2. His の側鎖（CG/ND1/CD2/CE1/NE2）を、構造内の実在 His（既定 144 番）の
     側鎖内部座標をテンプレートにして chi1/chi2 の格子で探索する
  3. 他の全原子は WT の座標のまま。動かすのは 184 番の側鎖だけ
  4. 立体障害が最小の配座と、それと十分に異なる第 2 配座を出力する
     （W243R で rot0/rot3 を比べたのと同じ頑健性確認をするため）

出力は重原子のみ（他の入力 PDB と同じ）。水素は haddock3 の topoaa が付ける。
"""
import math
import sys

import numpy as np

SIDE = ["CG", "ND1", "CD2", "CE1", "NE2"]   # 置き換える側鎖原子
KEEP = ["N", "CA", "C", "O", "CB"]          # WT から引き継ぐ原子
# 許容する最近接距離は元素の組み合わせで変える。
# N/O 同士は水素結合になりうるので 2.45 A まで許し、炭素・硫黄が絡む対は 3.00 A を下限にする。
# （WT でも ASP184 側鎖と GLU164 OE2 は 3.03 A しか離れていない＝この部位は元々詰まっている）
HB_MIN = 2.45
VDW_MIN = 3.00
CLASH_SOFT = 3.2      # この距離未満の接触数を数える（参考値）
STEP = 10             # chi の刻み（度）


def read_pdb(path):
    """PDB の ATOM/HETATM 行を (行, 残基番号, 原子名, 残基名, 座標) で返す。"""
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            out.append(dict(
                line=line.rstrip("\n"),
                resi=int(line[22:26]),
                name=line[12:16].strip(),
                resn=line[17:20].strip(),
                xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            ))
    return out


def kabsch(P, Q):
    """P を Q に重ねる回転と並進を返す。"""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, pc, qc


def rotate_about(axis_a, axis_b, angle_deg, points):
    """axis_a→axis_b の軸まわりに points を回す（Rodrigues）。"""
    k = axis_b - axis_a
    k = k / np.linalg.norm(k)
    t = math.radians(angle_deg)
    c, s = math.cos(t), math.sin(t)
    v = points - axis_a
    rot = v * c + np.cross(k, v) * s + np.outer(v @ k, k) * (1 - c)
    return rot + axis_a


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return math.degrees(math.atan2(np.cross(b1, v) @ w, v @ w))


def main():
    wt_path = sys.argv[1] if len(sys.argv) > 1 else "MRGPRX2_chainA.pdb"
    prefix = sys.argv[2] if len(sys.argv) > 2 else "MRGPRX2_D184H_chainA"
    target = 184
    template_resi = int(sys.argv[3]) if len(sys.argv) > 3 else 144

    atoms = read_pdb(wt_path)
    tgt = {a["name"]: a for a in atoms if a["resi"] == target}
    tpl = {a["name"]: a for a in atoms if a["resi"] == template_resi}
    assert tgt and tgt["N"]["resn"] == "ASP", f"{target} 番が ASP でない"
    assert tpl and tpl["N"]["resn"] == "HIS", f"{template_resi} 番が HIS でない（テンプレート）"
    for n in SIDE:
        assert n in tpl, f"テンプレート His に {n} が無い"

    # テンプレート His の N/CA/CB を 184 番の N/CA/CB に重ねる
    P = np.array([tpl[n]["xyz"] for n in ("N", "CA", "CB")])
    Q = np.array([tgt[n]["xyz"] for n in ("N", "CA", "CB")])
    R, pc, qc = kabsch(P, Q)
    base = {n: (tpl[n]["xyz"] - pc) @ R.T + qc for n in SIDE}

    # 184 番の側鎖以外すべて（判定の相手）
    others = np.array([a["xyz"] for a in atoms if a["resi"] != target])
    other_meta = [(a["resi"], a["resn"], a["name"]) for a in atoms if a["resi"] != target]
    # 同じ残基の主鎖も相手に含める（CB は結合相手なので除く）
    # 184 番自身の原子は相手に含めない。CG と CA のような結合を隔てた対は
    # 距離が幾何学的に決まっており（約 2.5 A）、衝突判定に入れると全配座が失格になる。
    # 側鎖内部の形はテンプレートにした実在の His の座標をそのまま使っているので担保される。
    partners = others
    # 元素ごとの許容下限（側鎖 5 原子 × 相手）
    p_el = np.array([m[2][0] for m in other_meta])
    s_el = np.array([n[0] for n in SIDE])
    allowed = np.where((np.isin(s_el, ["N", "O"])[:, None]) & (np.isin(p_el, ["N", "O"])[None, :]),
                       HB_MIN, VDW_MIN)

    ca, cb = tgt["CA"]["xyz"], tgt["CB"]["xyz"]
    results = []
    for chi1 in range(-180, 180, STEP):
        p1 = {n: rotate_about(ca, cb, chi1, base[n][None, :])[0] for n in SIDE}
        cg = p1["CG"]
        ring = [n for n in SIDE if n != "CG"]
        for chi2 in range(-180, 180, STEP):
            arr = rotate_about(cb, cg, chi2, np.array([p1[n] for n in ring]))
            pos = {"CG": cg}
            pos.update({n: arr[i] for i, n in enumerate(ring)})
            xyz = np.array([pos[n] for n in SIDE])
            d = np.linalg.norm(xyz[:, None, :] - partners[None, :, :], axis=2)
            margin = d - allowed          # 負なら許容下限を割っている
            dmin = float(d.min())
            nsoft = int((d < CLASH_SOFT).sum())
            # 実際の chi 角（テンプレート由来の初期値を含む値）を測り直す
            c1 = dihedral(tgt["N"]["xyz"], ca, cb, pos["CG"])
            c2 = dihedral(ca, cb, pos["CG"], pos["ND1"])
            i_w, j_w = np.unravel_index(margin.argmin(), margin.shape)
            results.append(dict(chi1=c1, chi2=c2, dmin=dmin, nsoft=nsoft, pos=pos,
                                near=other_meta[j_w], near_atom=SIDE[int(i_w)],
                                near_d=float(d[i_w, j_w]), worst=float(margin.min())))

    # 許容下限からの余裕が大きい順
    results.sort(key=lambda r: (-r["worst"], r["nsoft"]))
    print(f"探索した配座 {len(results)} 個。上位 8 個:")
    for r in results[:8]:
        ok = "可" if r["worst"] >= 0 else "衝突"
        print(f"   chi1={r['chi1']:7.1f} chi2={r['chi2']:7.1f}  最近接={r['dmin']:5.2f} A  "
              f"基準を最も割る対: {r['near_atom']}–{r['near'][1]}{r['near'][0]} {r['near'][2]} "
              f"{r['near_d']:.2f} A (余裕{r['worst']:+5.2f})  <{CLASH_SOFT}A 接触={r['nsoft']:3d}  {ok}")
    feasible = [r for r in results if r["worst"] >= 0]
    print(f"許容下限（N/O 同士 {HB_MIN} A・それ以外 {VDW_MIN} A）を満たす配座: {len(feasible)} 個")
    if feasible:
        feasible.sort(key=lambda r: (r["nsoft"], -r["worst"]))
        results = feasible
    else:
        # 基準を満たす配座が無い＝His184 はこのポケットで本質的に窮屈。
        # 事実として記録したうえで、違反が最小の配座を採用する（下流の flexref で緩和される）。
        print("※ 基準を満たす配座は無い。違反が最小の配座を採用する。"
              " この窮屈さ自体が結果なので、事前登録と報告に明記すること。")
    best = results[0]

    # 第 2 配座: chi1 か chi2 が 60 度以上離れているもののうち最良
    def far(a, b):
        return min(abs(a - b), 360 - abs(a - b))
    alt = None
    for r in results[1:]:
        if far(r["chi1"], best["chi1"]) >= 60 or far(r["chi2"], best["chi2"]) >= 60:
            alt = r
            break

    for tag, r in (("rot1", best), ("rot2", alt)):
        if r is None:
            print(f"{tag}: 十分に異なる配座が見つからない")
            continue
        out = f"{prefix}.pdb" if tag == "rot1" else f"{prefix}_rot2.pdb"
        write(atoms, target, r["pos"], out)
        print(f"{tag} -> {out}")
        print(f"   chi1={r['chi1']:7.1f}  chi2={r['chi2']:7.1f}  "
              f"最近接={r['dmin']:.2f} A ({r['near'][1]}{r['near'][0]} {r['near'][2]})  "
              f"{CLASH_SOFT} A 未満の接触={r['nsoft']}")


def write(atoms, target, pos, out):
    """184 番の側鎖だけ差し替えた PDB を書く。他の行は元のまま（原子番号のみ振り直す）。"""
    tpl_line = next(a["line"] for a in atoms if a["resi"] == target and a["name"] == "CB")
    lines = []
    for a in atoms:
        if a["resi"] != target:
            lines.append(a["line"])
            continue
        if a["name"] not in KEEP:
            continue                      # Asp の CG/OD1/OD2 は落とす
        lines.append(a["line"][:17] + "HIS" + a["line"][20:])
        if a["name"] == "CB":             # CB の直後に His 側鎖を挿入
            for n in SIDE:
                x, y, z = pos[n]
                el = n[0]                 # C か N
                line = (tpl_line[:12] + f" {n:<3s}" + tpl_line[16:17] + "HIS" + tpl_line[20:30]
                        + f"{x:8.3f}{y:8.3f}{z:8.3f}" + tpl_line[54:76] + f"{el:>2s}")
                lines.append(line)
    with open(out, "w") as f:
        for i, line in enumerate(lines, start=1):   # 原子番号を振り直す
            f.write(f"{line[:6]}{i:5d}{line[11:]}\n")
        f.write("END\n")
    # 原子名の重複や欠落が無いか確認する
    chk = read_pdb(out)
    names = [a["name"] for a in chk if a["resi"] == target]
    assert sorted(names) == sorted(KEEP + SIDE), f"184 番の原子が想定と違う: {names}"


if __name__ == "__main__":
    main()
