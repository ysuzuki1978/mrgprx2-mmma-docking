#!/usr/bin/env python3
"""レター本文に書く数値が、data/tables/ の内容と一致するかを検査する回帰テスト。

使い方: python3 scripts/verify.py     (letter_MMA/ 直下)
本文の数値を変えたら、ここの EXPECTED も更新すること。
不一致が出たら「本文が古い」か「データが変わった」のどちらか。必ず原因を特定する。
"""
import csv, os, sys
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T=os.path.join(ROOT,"data","tables")

def load(name, key=None):
    with open(os.path.join(T,name)) as f:
        # 先頭の `#` 行は表の注記なので飛ばす（S6 など）
        return list(csv.DictReader((l for l in f if not l.startswith("#")),delimiter="\t"))

# (説明, 実測値を取る関数, 期待値, 許容誤差)
def g7(lig,rec,col):
    for r in load("T7_D184_dissection.tsv"):
        if r["ligand"]==lig and r["receptor"]==rec: return float(r[col])
    raise KeyError((lig,rec,col))
def g8(m):
    for r in load("T8_binding_mode_comparison.tsv"):
        if r["measure"]==m: return r["value"]
    raise KeyError(m)
def g9(lig,rec,col):
    """T9_D184_ladder.tsv（レター v2 用、5 反復中央値・正しい三量体）から引く。"""
    for r in load("T9_D184_ladder.tsv"):
        if r["ligand"]==lig and r["receptor"]==rec: return r[col]
    raise KeyError((lig,rec,col))
def g9f(lig,rec,col):
    return float(g9(lig,rec,col))

def g1(run,col):
    for r in load("T1_energies.tsv"):
        if r["run"]==run: return float(r[col])
    raise KeyError((run,col))

def s6_ok():
    """Table S6 の MMMA 20 ポーズのうち、弱い C-H...O の基準を満たす数。"""
    n = 0
    for r in load("S6_hbond_geometry.tsv"):
        if r["ligand"] != "MMMA" or not r["angle_D_H_O_deg"]:
            continue
        if 110 <= float(r["angle_D_H_O_deg"]) <= 180 and 2.0 <= float(r["d_H_O_A"]) <= 3.0:
            n += 1
    return n


def offset_pair():
    """静電損失を vdW がどれだけ埋めるか。(中央値の差での %, seed ごとの範囲)。

    項ごとの中央値は加法的でないので、この 2 つは一致しない。図 1 の注記はその両方を出す。
    """
    import statistics as st
    import sys
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import make_ladder_table as ml

    def terms(pre):
        out = []
        for s in ml.SEEDS:
            r = ml.rank1(pre + s)
            if r:
                out.append((float(r["vdw"]), ml.W_EMREF * abs(float(r["elec"]))))
        return out

    wt, dn = terms("MMAREP_MMT_blind"), terms("MMAREP_MMT_D184N_blind")
    med = 100 * (st.median(w[0] for w in wt) - st.median(d[0] for d in dn)) \
              / (st.median(w[1] for w in wt) - st.median(d[1] for d in dn))
    per = sorted(round(100 * (w[0] - d[0]) / (w[1] - d[1])) for w, d in zip(wt, dn))
    return med, f"{per[0]}-{per[-1]}"


def pvalue_exponent():
    """本文の p 値の指数。F(1,32) = 24.55 なら p = 2.266e-05 なので -5。

    2026-11-02: 文献番号を繰り下げるスクリプトの正規表現が `10⁻⁵` の `⁵` を引用番号と
    誤認し、黙って `10⁻⁶` に書き換えた。上付き文字は引用番号だけとは限らないので、
    本文から直接読み取って見張る。
    """
    import re
    body = open(os.path.join(ROOT, "DRAFT_letter_v2.md"), encoding="utf-8").read()
    m = re.search(r"p = 2\.26 × 10([⁻⁻])([⁰¹²³⁴⁵⁶⁷⁸⁹]+)", body)
    if not m:
        raise KeyError("本文に p = 2.26 × 10^-n が見つからない")
    digits = {c: str(i) for i, c in enumerate("⁰¹²³⁴⁵⁶⁷⁸⁹")}
    return "-" + "".join(digits[c] for c in m.group(2))


def max48():
    """compound 48/80 が全条件で示した重心移動の最大値（T9 の範囲列の上端）。"""
    return max(float(r["centroid_shift_range_A"].split("-")[1])
               for r in load("T9_D184_ladder.tsv") if r["ligand"] == "Compound4880")


def s7max(lig, rec):
    """Table S7（ポーズごとの Jaccard）の最大値を、表に書いてある文字列のまま返す。"""
    v = [r["jaccard_vs_wt"] for r in load("S7_jaccard_per_pose.tsv")
         if r["ligand"] == lig and r["receptor"] == rec]
    if not v:
        raise KeyError((lig, rec))
    return max(v, key=float)

def s6(lig):
    """Table S6（Asp184 への水素結合の幾何、20 ポーズ）の行を返す。"""
    rows=[r for r in load("S6_hbond_geometry.tsv") if r["ligand"]==lig]
    if not rows: raise KeyError(lig)
    return rows
def s6rng(lig,col):
    v=[float(r[col]) for r in s6(lig)]
    d=1 if col.endswith("_A") else 0
    w=2 if d else 0
    return f"{min(v):.{d*2}f}-{max(v):.{d*2}f}" if d else f"{min(v):.0f}-{max(v):.0f}"
def s6donors(lig):
    """最短の H⋯O を作る供与体の元素。MMMA は C だけ（＝水素結合供与体を持たない）。"""
    return ",".join(sorted({r["donor"].split("(")[1].rstrip(")") for r in s6(lig)}))

def s5rng(lig,rec,col,resi="184"):
    """Table S5（酸性ペアへの距離、ポーズ単位）のレンジ。"""
    v=[float(r[col]) for r in load("S5_distance_to_acidic_pair.tsv")
       if r["ligand"]==lig and r["receptor"]==rec and r["residue_number"]==resi and r[col]]
    return f"{min(v):.2f}-{max(v):.2f}"
def s8(comparison):
    for r in load("S8_pose_rmsd.tsv"):
        if r["comparison"]==comparison: return f'{float(r["ligand_rmsd_A"]):.2f}'
    raise KeyError(comparison)
def s9(pdb,res):
    for r in load("S9_cryoem_contacts.tsv"):
        if r["pdb"]==pdb and r["receptor_residue"]==res: return r["min_distance_A"]
    raise KeyError((pdb,res))

CHECKS=[
 ("MMA WT score",                     lambda: g1("MMT_blind","score"),            -19.235, 0.001),
 ("MMA WT vdW % of score",            lambda: g1("MMT_blind","vdw_pct_of_score"),  81.4,   0.15),
 ("MMA WT elec % of score",           lambda: g1("MMT_blind","elec_pct_of_score"), 12.4,   0.15),
 ("C48/80 WT elec % of score",        lambda: g1("Compound4880_blind","elec_pct_of_score"), 53.7, 0.15),
 ("C48/80 D184N d_score",             lambda: g7("Compound4880","D184N","d_score"), 13.697, 0.01),
 ("C48/80 D184N d_score % of WT",     lambda: g7("Compound4880","D184N","d_score_pct_of_WT"), 12.3, 0.1),
 ("C48/80 D184A d_score",             lambda: g7("Compound4880","D184A","d_score"), 12.288, 0.01),
 ("MMA D184N d_score",                lambda: g7("MMT","D184N","d_score"),          -0.238, 0.01),
 ("MMA D184A d_score",                lambda: g7("MMT","D184A","d_score"),           0.818, 0.01),
 ("MMA D184A d_score % of WT",        lambda: g7("MMT","D184A","d_score_pct_of_WT"),  4.3,  0.1),
 ("MMA E164A d_score",                lambda: g7("MMT","E164A","d_score"),          -1.365, 0.01),
 ("MMA aromatic contact %",           lambda: float(g8("MMA_aromatic_pct")),         68.2,  0.15),
 ("C48/80 aromatic contact %",        lambda: float(g8("C4880_aromatic_pct")),       35.6,  0.15),
 ("MMA acidic contact %",             lambda: float(g8("MMA_acidic_pct")),           14.9,  0.15),
 ("C48/80 acidic contact %",          lambda: float(g8("C4880_acidic_pct")),         15.3,  0.15),
 ("MMA atoms within 2A of C48/80 (%)",lambda: float(g8("MMA_atoms_within_2A_of_C4880_volume_pct")), 100, 0.1),
]
# ---- レター v2 (DRAFT_letter_v2.md) が引用する数値 ----
# 48/80 は正しい三量体（6IB 結合配座）、各条件 5 反復の中央値。出所は T9_D184_ladder.tsv。
CHECKS_V2=[
 ("v2 C48/80 WT elec % of score", lambda: g9f("Compound4880","WT","elec_pct_of_score"), 54.2, 0.05),
 ("v2 C48/80 WT vdW % of score",  lambda: g9f("Compound4880","WT","vdw_pct_of_score"),  44.3, 0.05),
 ("v2 MMA WT elec % of score",    lambda: g9f("MMT","WT","elec_pct_of_score"),          12.6, 0.05),
 ("v2 MMA WT vdW % of score",     lambda: g9f("MMT","WT","vdw_pct_of_score"),           79.6, 0.05),
 ("v2 C48/80 D184N Δelec %",      lambda: g9f("Compound4880","D184N","d_elec_contribution_pct_of_WT"), -24.0, 0.05),
 ("v2 C48/80 D184N Δscore %",     lambda: g9f("Compound4880","D184N","d_score_pct_of_WT"),             -11.9, 0.05),
 ("v2 C48/80 D184A Δelec %",      lambda: g9f("Compound4880","D184A","d_elec_contribution_pct_of_WT"), -23.0, 0.05),
 ("v2 C48/80 D184A Δscore %",     lambda: g9f("Compound4880","D184A","d_score_pct_of_WT"),             -13.3, 0.05),
 ("v2 C48/80 D184E rot1 Δelec %", lambda: g9f("Compound4880","D184E_rot1","d_elec_contribution_pct_of_WT"), -2.5, 0.05),
 ("v2 C48/80 D184E rot2 Δelec %", lambda: g9f("Compound4880","D184E_rot2","d_elec_contribution_pct_of_WT"), -0.6, 0.05),
 ("v2 C48/80 D184H rot1P Δscore %", lambda: g9f("Compound4880","D184H_rot1P","d_score_pct_of_WT"), -22.5, 0.05),
 ("v2 C48/80 D184H rot2P Δscore %", lambda: g9f("Compound4880","D184H_rot2P","d_score_pct_of_WT"), -24.2, 0.05),
 ("v2 MMA WT score",              lambda: g9f("MMT","WT","median_score"),              -19.367, 0.001),
 ("v2 MMA D184N score",           lambda: g9f("MMT","D184N","median_score"),           -19.457, 0.001),
 ("v2 MMA D184A Δscore %",        lambda: g9f("MMT","D184A","d_score_pct_of_WT"),        -5.1,  0.05),
 ("v2 MMA E164A Δscore (units)",  lambda: -g9f("MMT","E164A","d_score"),                  1.179, 0.005),
 # 本文 ¶4 は 2026-11-01 に score 単位 (1.18) から % 表記 (+6.1%) に統一した
 ("v2 MMA E164A Δscore %",        lambda: g9f("MMT","E164A","d_score_pct_of_WT"),          6.1,  0.05),
 # 2026-10-24 追加（事前登録 section 17）。yakuri3 の MMMA。
 ("v2 MMA WT score (yakuri3)",    lambda: g9f("MMT_y3","WT","median_score"),            -19.516, 0.001),
 ("v2 MMA D184E rot1 Δscore %",   lambda: g9f("MMT_y3","D184E_rot1","d_score_pct_of_WT"), -11.6, 0.05),
 ("v2 MMA D184E rot2 Δscore %",   lambda: g9f("MMT_y3","D184E_rot2","d_score_pct_of_WT"), -10.3, 0.05),
]
# 部位の保持（Figure 2B の "stays at the site" 列と本文の根拠）
STR_CHECKS_V2=[
 ("v2 C48/80 D184N site kept",  lambda: g9("Compound4880","D184N","site_kept_of_20"), "20/20"),
 ("v2 MMA D184N site kept",     lambda: g9("MMT","D184N","site_kept_of_20"),          "20/20"),
 ("v2 MMA D184A site kept",     lambda: g9("MMT","D184A","site_kept_of_20"),          "16/20"),
 ("v2 MMA D184E rot1 site kept",lambda: g9("MMT_y3","D184E_rot1","site_kept_of_20"),   "4/20"),
 ("v2 MMA D184E rot2 site kept",lambda: g9("MMT_y3","D184E_rot2","site_kept_of_20"),  "16/20"),
 ("v2 MMA D184H rot1P site kept",lambda: g9("MMT","D184H_rot1P","site_kept_of_20"),    "0/20"),
 ("v2 MMA D184H rot1N site kept",lambda: g9("MMT","D184H_rot1N","site_kept_of_20"),    "0/20"),
 ("v2 MMA D184H rot2P site kept",lambda: g9("MMT","D184H_rot2P","site_kept_of_20"),   "20/20"),
 # 本文 ¶5「compound 48/80 は 20 ポーズ中 9〜13 で部位を保つ」の根拠（3 モデル）
 ("v2 C48/80 D184H rot1P site kept",lambda: g9("Compound4880","D184H_rot1P","site_kept_of_20"),"13/20"),
 ("v2 C48/80 D184H rot1N site kept",lambda: g9("Compound4880","D184H_rot1N","site_kept_of_20"), "9/20"),
 ("v2 C48/80 D184H rot2P site kept",lambda: g9("Compound4880","D184H_rot2P","site_kept_of_20"),"12/20"),
 # 本文 ¶5「3 モデル中 2 つは共有残基ゼロ、残り 1 つは 0.44-0.56」
 ("v2 MMA D184H rot1P jaccard",  lambda: g9("MMT","D184H_rot1P","jaccard_vs_WT_range"), "0.000-0.000"),
 ("v2 MMA D184H rot1N jaccard",  lambda: g9("MMT","D184H_rot1N","jaccard_vs_WT_range"), "0.000-0.000"),
 ("v2 MMA D184H rot2P jaccard",  lambda: g9("MMT","D184H_rot2P","jaccard_vs_WT_range"), "0.625-0.750"),
 # 本文 ¶4「E164A では 20 ポーズ中 13 でしか部位を保たない」
 ("v2 MMA E164A site kept",      lambda: g9("MMT","E164A","site_kept_of_20"),           "20/20"),
 # 本文 ¶5 の 0.556。2 桁だと 0.56 と出て閾値 0.56 と矛盾して見えるので S7 は 3 桁
 ("v2 MMA D184H rot2P Jaccard max", lambda: s7max("MMT","D184H_rot2P"),                 "0.750"),
 # 2026-11-02 追加。接触集合は変異で残基の形が変わると下がりうるので、
 # 本文は「動いたかどうか」を受容体 Cα 重ね合わせ後の重心移動で述べる。
 ("v2 MMA D184N 重心移動",        lambda: g9("MMT","D184N","centroid_shift_median_A"),      "0.98"),
 ("v2 MMA D184N 重心移動の範囲",   lambda: g9("MMT","D184N","centroid_shift_range_A"),  "0.07-1.29"),
 ("v2 MMA D184N 2A 以内",         lambda: g9("MMT","D184N","poses_within_2A_of_20"),       "20/20"),
 ("v2 C48/80 D184N 重心移動",     lambda: g9("Compound4880","D184N","centroid_shift_median_A"), "0.68"),
 ("v2 MMA D184H rot1P 重心移動",  lambda: g9("MMT","D184H_rot1P","centroid_shift_median_A"), "23.58"),
 ("v2 MMA D184H rot1N 重心移動",  lambda: g9("MMT","D184H_rot1N","centroid_shift_median_A"), "23.72"),
 ("v2 MMA D184H rot2P 重心移動",  lambda: g9("MMT","D184H_rot2P","centroid_shift_median_A"),  "3.10"),
 ("v2 MMA D184E rot1 重心移動",   lambda: g9("MMT_y3","D184E_rot1","centroid_shift_median_A"),"13.04"),
 ("v2 MMA D184E rot2 重心移動",   lambda: g9("MMT_y3","D184E_rot2","centroid_shift_median_A"), "2.19"),
 ("v2 MMA D184A 重心移動",        lambda: g9("MMT","D184A","centroid_shift_median_A"),        "1.10"),
 ("v2 MMA E164A 重心移動",        lambda: g9("MMT","E164A","centroid_shift_median_A"),        "1.83"),
 # 本文 ¶5「48/80 はどの置換でも 2.6 A 以上は動かない」
 ("v2 C48/80 重心移動の最大",     lambda: f"{max48():.2f}",                                   "2.52"),
 # SI の「18/20 が弱い C-H...O の基準を満たす」。残り 2 つは角度が 110 度をわずかに下回る
 ("S6 MMA 基準を満たすポーズ",    lambda: f"{s6_ok()}/20",                                   "18/20"),
 # 図 1 注記の補償率。中央値の差で 87%、seed ごとだと 34-121% に散る
 ("Fig1 vdW 補償（中央値の差）",  lambda: f"{offset_pair()[0]:.0f}",                            "87"),
 ("Fig1 vdW 補償（seed ごと）",   lambda: offset_pair()[1],                                 "34-121"),
]
# 原稿そのものを見る検査は、原稿が隣にあるときだけ走らせる。
# 寄託アーカイブには scripts/ と data/ しか入らないので、そこでは黙って飛ばす。
# 2025 年の論文の交互作用の p 値。F(1,32) = 24.55 -> p = 2.266e-05
if os.path.exists(os.path.join(ROOT, "DRAFT_letter_v2.md")):
    STR_CHECKS_V2.append(("¶1 の p 値の指数", pvalue_exponent, "-5"))
STR_CHECKS=[
 ("MMA-only contact residues", lambda: g8("contacts_MMA_only"), "(none)"),
 ("shared contact residues",   lambda: g8("contacts_shared"),   "CYS168;PHE170;ASP184;TRP243;TRP248"),
 ("v2 C48/80 WT 184 contacts",   lambda: g9("Compound4880","WT","contacts_184_within_4A"),   "20/20"),
 ("v2 C48/80 WT 184 dist (A)",   lambda: g9("Compound4880","WT","dist_184_range_A"),         "2.62-2.76"),
 ("v2 C48/80 D184H rot1P 184",   lambda: g9("Compound4880","D184H_rot1P","contacts_184_within_4A"), "0/20"),
 ("v2 C48/80 D184H rot1N 184",   lambda: g9("Compound4880","D184H_rot1N","contacts_184_within_4A"), "3/20"),
 ("v2 C48/80 D184H rot2P 184",   lambda: g9("Compound4880","D184H_rot2P","contacts_184_within_4A"), "0/20"),
 ("v2 MMA WT 184 dist (A)",      lambda: g9("MMT","WT","dist_184_range_A"),                  "3.08-3.84"),
 ("v2 C48/80 D184E rot1 184",    lambda: g9("Compound4880","D184E_rot1","contacts_184_within_4A"), "18/20"),
 ("v2 C48/80 D184E rot2 184",    lambda: g9("Compound4880","D184E_rot2","contacts_184_within_4A"), "19/20"),
 # 3 段落目の水素結合の幾何（Table S6。make_s6_hbond_geometry.py が生成する）
 ("v2 MMA C-H...O  H...O (A)",   lambda: s6rng("MMMA","d_H_O_A"),                "2.42-3.00"),
 ("v2 MMA C-H...O  角 (deg)",    lambda: s6rng("MMMA","angle_D_H_O_deg"),        "106-150"),
 ("v2 C48/80 N-H...O  H...O (A)",lambda: s6rng("compound 48/80","d_H_O_A"),      "1.62-1.86"),
 ("v2 C48/80 N-H...O  角 (deg)", lambda: s6rng("compound 48/80","angle_D_H_O_deg"), "139-165"),
 ("v2 MMA に水素結合供与体なし", lambda: s6donors("MMMA"),                        "C"),
 # 2 段落目のポーズ RMSD（Table S8）と cryo-EM の距離（Table S9）
 ("v2 無拘束 vs 拘束 RMSD (A)",  lambda: s8("blind vs restrained, both top-ranked"),      "1.02"),
 ("v2 電荷モデル差替 RMSD (A)",  lambda: s8("AM1-BCC top vs Gasteiger rank 3"),           "0.64"),
 ("v2 cryo-EM 7VV6 D184 (A)",    lambda: s9("7VV6","ASP184"),                            "2.47"),
 # 3 段落目の MMMA 酸素の距離（Table S5）
 ("v2 MMA WT 最近接酸素 D184",   lambda: s5rng("MMT","WT","nearest_oxygen_distance_A","184"), "5.20-5.92"),
 ("v2 MMA WT 最近接酸素 E164",   lambda: s5rng("MMT","WT","nearest_oxygen_distance_A","164"), "6.55-9.20"),
]
def main():
    bad=0
    print(f"{'項目':38s}{'実測':>12s}{'期待':>12s}  判定")
    print("-"*80)
    for lab,fn,exp,tol in CHECKS + CHECKS_V2:
        try: v=fn()
        except Exception as e: print(f"{lab:38s}{'ERR':>12s}{exp:>12}  !! {e}"); bad+=1; continue
        ok=abs(v-exp)<=tol
        bad += 0 if ok else 1
        print(f"{lab:38s}{v:12.3f}{exp:12.3f}  {'OK' if ok else '*** 不一致 ***'}")
    for lab,fn,exp in STR_CHECKS + STR_CHECKS_V2:
        v=fn(); ok=(v==exp); bad += 0 if ok else 1
        print(f"{lab:38s}{v[:28]:>28s}  {'OK' if ok else '*** 不一致: 期待 '+exp+' ***'}")
    print()
    # seed ノイズ床が本当に小さいか
    for r in load("T6_seed_replicates.tsv"):
        if r["run"]=="SUMMARY":
            print(f"seed ノイズ ({r['ligand']}): {r['score']} {r['vdw']} {r['elec']} {r['top_cluster_population']}")
    print()
    print(f"{'すべて一致' if bad==0 else str(bad)+' 件が不一致'}")
    sys.exit(1 if bad else 0)
if __name__=="__main__": main()
