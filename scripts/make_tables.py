#!/usr/bin/env python3
"""レターに載せる全数値を data/runs/ から再生成し、data/tables/*.tsv に書く。

使い方:  python3 scripts/make_tables.py     (letter_MMA/ 直下で実行)
このスクリプトが唯一の数値の出所。本文・図の数字はここから引くこと。
"""
import csv, glob, gzip, os, sys
import numpy as np
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "data", "runs")
RAW  = os.path.join(ROOT, "data", "raw")
OUT  = os.path.join(ROOT, "data", "tables")
os.makedirs(OUT, exist_ok=True)

AROM = {"PHE","TYR","TRP","HIS"}
ACID = {"ASP","GLU"}
BB   = {"N","CA","C","O","OXT"}
W_EMREF, W_FLEXREF = 0.2, 1.0

MMT_ATOM = {"C":"a-CH3","C1":"a-C(sp2)","C2":"=CH2","C3":"C=O carbon",
            "O":"C=O oxygen","O1":"ester O","C4":"O-CH3"}
MMT_CHG  = {"C":0.102,"C1":-0.151,"C2":0.119,"C3":0.632,
            "O":-0.539,"O1":-0.443,"C4":0.280}
DEEP = {164,168,170,180,184,243,244,247,248}
VEST = {19,20,21,22,25,103,169,251,252,253}

ALL_RUNS = ["MMT_blind","MMT_blind_seed1000","MMT_blind_seed2000","MMT_blind_seed3000",
            "MMT_D184A_blind","MMT_D184N_blind","MMT_E164A_blind","MMT_chainB","MMG_chainB",
            "Compound4880_blind","Compound4880_blind_seed1000","Compound4880_blind_seed2000",
            "Compound4880_blind_seed3000","Compound4880_D184A_blind","Compound4880_D184N_blind","PSB-172656_blind"]

def read_pdb(p):
    op = gzip.open if p.endswith(".gz") else open
    a=[]
    with op(p,"rt") as f:
        for l in f:
            if l.startswith(("ATOM","HETATM")):
                el=(l[76:78].strip() or l[12:16].strip()[0]).upper()
                a.append(dict(name=l[12:16].strip(),resn=l[17:20].strip(),ch=l[21],
                              resi=int(l[22:26]),el=el,
                              xyz=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
    return a

def models(run):
    return sorted(glob.glob(os.path.join(RUNS,run,"cluster_1_model_*.pdb.gz")))

def emref(run):
    c=os.path.join(RUNS,run,run+".cfg")
    return "[emref]" in open(c).read() if os.path.exists(c) else run.endswith("blind")

def rank1(run):
    p=os.path.join(RUNS,run,"capri_clt.tsv")
    if not os.path.exists(p): return None
    L=[l for l in open(p) if not l.startswith("#")]
    if len(L)<2: return None
    h=L[0].rstrip("\n").split("\t")
    for l in L[1:]:
        d=dict(zip(h,l.rstrip("\n").split("\t")))
        if d.get("cluster_rank")=="1": return d
    return None

def population(run):
    """rank1 クラスターの母集団。
    NOTE: seletopclusts の出力ファイル名は cluster_RANK であり cluster_id ではない。
    この取り違えが過去の解析で 47% のランに影響した既知のバグ（本体報告書 6.1）。"""
    p=os.path.join(RUNS,run,"cluster_assignment.tsv"); r=rank1(run)
    if not (os.path.exists(p) and r): return "",""
    pops=Counter()
    with open(p) as f:
        for row in csv.DictReader(f,delimiter="\t"):
            cid=(row.get("cluster_id") or "-").strip()
            if cid not in ("-","","None"): pops[cid]+=1
    return pops.get(r.get("cluster_id",""),""), sum(pops.values())

def kabsch(P,Q):
    pc,qc=P.mean(0),Q.mean(0)
    U,S,Vt=np.linalg.svd((P-pc).T@(Q-qc))
    d=np.sign(np.linalg.det(Vt.T@U.T))
    return Vt.T@np.diag([1,1,d])@U.T,pc,qc

REF=read_pdb(os.path.join(RAW,"MRGPRX2_chainA.pdb"))
REFCA={a["resi"]:a["xyz"] for a in REF if a["ch"]=="A" and a["name"]=="CA"}

def framed(run):
    out=[]
    for p in models(run):
        at=read_pdb(p)
        ca={a["resi"]:a["xyz"] for a in at if a["ch"]=="A" and a["name"]=="CA"}
        com=sorted(set(ca)&set(REFCA))
        R,pc,qc=kabsch(np.array([ca[i] for i in com]),np.array([REFCA[i] for i in com]))
        for a in at: a["fit"]=(a["xyz"]-pc)@R.T+qc
        out.append(at)
    return out

def w(name,header,rows):
    with open(os.path.join(OUT,name),"w",newline="") as f:
        wr=csv.writer(f,delimiter="\t"); wr.writerow(header); wr.writerows(rows)
    print(f"  {name}  ({len(rows)} 行)")

def t1():
    rows=[]
    for r in ALL_RUNS:
        d=rank1(r)
        if not d: continue
        we=W_EMREF if emref(r) else W_FLEXREF
        s,v,e,ds,b=(float(d[k]) for k in ("score","vdw","elec","desolv","bsa"))
        pop,tot=population(r)
        rows.append([r,"emref" if emref(r) else "flexref",f"{s:.3f}",d.get("score_std",""),
                     f"{v:.3f}",f"{e:.3f}",f"{ds:.3f}",f"{b:.1f}",f"{we:.1f}",f"{we*e:.3f}",
                     f"{100*we*e/s:.1f}",f"{100*v/s:.1f}",f"{100*ds/s:.1f}",pop,tot])
    w("T1_energies.tsv",["run","terminal_module","score","score_sd_within_cluster","vdw","elec",
       "desolv","bsa","w_elec","elec_contribution","elec_pct_of_score","vdw_pct_of_score",
       "desolv_pct_of_score","top_cluster_population","models_clustered"],rows)

def t2():
    rows=[]
    for r in ALL_RUNS:
        F=framed(r)
        if not F: continue
        rn=next((a["resn"] for a in F[0] if a["ch"]=="A" and a["resi"]==184),"")
        ds,deep,vest=[],Counter(),Counter()
        for at in F:
            L=np.array([a["fit"] for a in at if a["ch"]=="B" and a["el"]!="H"])
            side=[a["fit"] for a in at if a["ch"]=="A" and a["resi"]==184
                  and a["el"]!="H" and a["name"] not in BB]
            if not side:
                side=[a["fit"] for a in at if a["ch"]=="A" and a["resi"]==184 and a["el"]!="H"]
            ds.append(float(np.min(np.linalg.norm(L[:,None,:]-np.array(side)[None,:,:],axis=2))))
            for a in at:
                if a["ch"]!="A" or a["el"]=="H": continue
                if np.min(np.linalg.norm(L-a["fit"],axis=1))<=4.0:
                    if a["resi"] in DEEP: deep[a["resi"]]+=1
                    if a["resi"] in VEST: vest[a["resi"]]+=1
        nd,nv=len(deep),len(vest)
        rows.append([r,rn," ".join(f"{x:.2f}" for x in ds),
                     sum(1 for x in ds if x<=4.0),len(ds),nd,nv,
                     f"{nd/(nd+nv):.2f}" if nd+nv else ""])
    w("T2_pocket_access.tsv",["run","residue184_identity","dist_to_res184_per_model_A",
       "n_models_within_4A","n_models","n_deep_residues","n_vestibule_residues","deep_fraction"],rows)

def t3():
    rows=[]
    for r in ALL_RUNS:
        F=framed(r)
        if not F: continue
        c,nm=Counter(),{}
        for at in F:
            L=np.array([a["fit"] for a in at if a["ch"]=="B" and a["el"]!="H"])
            seen=set()
            for a in at:
                if a["ch"]!="A" or a["el"]=="H": continue
                if np.min(np.linalg.norm(L-a["fit"],axis=1))<=4.0:
                    seen.add(a["resi"]); nm[a["resi"]]=a["resn"]
            for x in seen: c[x]+=1
        for resi,n in sorted(c.items(),key=lambda x:(-x[1],x[0])):
            rows.append([r,nm[resi],resi,n,len(F),
                         "aromatic" if nm[resi] in AROM else
                         "acidic" if nm[resi] in ACID else "other"])
    w("T3_contacts_residue.tsv",["run","resname","resi","n_models_in_contact","n_models",
                                 "residue_class"],rows)

def t4():
    rows=[]
    for r in [x for x in ALL_RUNS if x.startswith(("MMT","Compound4880"))]:
        F=framed(r)
        if not F: continue
        tot,ar,mind=defaultdict(list),defaultdict(list),defaultdict(list)
        for at in F:
            rec=[a for a in at if a["ch"]=="A" and a["el"]!="H"]
            R=np.array([a["fit"] for a in rec])
            for a in at:
                if a["ch"]!="B" or a["el"]=="H": continue
                d=np.linalg.norm(R-a["fit"],axis=1)
                tot[a["name"]].append(int((d<=4.5).sum()))
                ar[a["name"]].append(int(sum(1 for i,x in enumerate(rec)
                                             if d[i]<=4.5 and x["resn"] in AROM)))
                mind[a["name"]].append(float(d.min()))
        T=sum(np.mean(v) for v in tot.values())
        for k in sorted(tot,key=lambda x:-np.mean(tot[x])):
            rows.append([r,k,MMT_ATOM.get(k,""),MMT_CHG.get(k,""),
                         f"{np.mean(tot[k]):.1f}",f"{100*np.mean(tot[k])/T:.1f}",
                         f"{np.mean(ar[k]):.1f}",f"{np.mean(mind[k]):.2f}"])
    w("T4_contacts_per_ligand_atom.tsv",["run","atom_name","atom_role","am1bcc_charge",
       "mean_contacts_within_4.5A","pct_of_total_contacts","mean_aromatic_contacts",
       "mean_min_dist_A"],rows)

def t5():
    rows=[]
    for r in ALL_RUNS:
        F=framed(r)
        if not F: continue
        sb=0; best=[]
        for at in F:
            lig=[a for a in at if a["ch"]=="B" and a["el"]!="H"]
            for a in lig:
                if a["el"]!="N": continue
                for b in at:
                    if (b["ch"]=="A" and b["resn"] in ACID
                            and b["name"] in ("OD1","OD2","OE1","OE2")
                            and np.linalg.norm(a["fit"]-b["fit"])<=4.0):
                        sb+=1; break
            o=[b for b in at if b["ch"]=="A" and b["resi"]==184 and b["name"] in ("OD1","OD2")]
            if o:
                best.append(min(((float(np.linalg.norm(x["fit"]-a["fit"])),a["name"],a["el"])
                                 for x in o for a in lig)))
        if best:
            md=float(np.mean([x[0] for x in best]))
            cls=("salt bridge" if best[0][2]=="N" and md<=4.0
                 else "weak C-H...O" if best[0][2]=="C" else "-")
            rows.append([r,"closest_ligand_atom_to_D184_carboxylate",best[0][1],best[0][2],
                         f"{md:.2f}",len(F),cls])
        rows.append([r,"salt_bridge_ligandN_to_AspGlu_count","","","",len(F),
                     f"{sb} model-instances"])
    w("T5_interactions.tsv",["run","measure","ligand_atom","element","mean_distance_A",
                             "n_models","classification"],rows)

def t6():
    rows=[]
    for lig,rs in (("MMA",["MMT_blind","MMT_blind_seed1000","MMT_blind_seed2000","MMT_blind_seed3000"]),
                   ("C48/80",["Compound4880_blind","Compound4880_blind_seed1000",
                              "Compound4880_blind_seed2000","Compound4880_blind_seed3000"])):
        sc=[]
        for r in rs:
            d=rank1(r)
            if not d: continue
            cfg=os.path.join(RUNS,r,r+".cfg"); seed=""
            if os.path.exists(cfg):
                for l in open(cfg):
                    if "iniseed" in l: seed=l.split("=")[1].strip()
            pop,tot=population(r); sc.append(float(d["score"]))
            rows.append([lig,r,seed or "default",f"{float(d['score']):.3f}",
                         f"{float(d['vdw']):.3f}",f"{float(d['elec']):.3f}",pop,tot])
        if len(sc)>1:
            a=np.array(sc)
            rows.append([lig,"SUMMARY","",f"mean={a.mean():.3f}",f"sd={a.std(ddof=1):.4f}",
                         f"range={a.max()-a.min():.4f}",f"n={len(a)}",""])
    w("T6_seed_replicates.tsv",["ligand","run","iniseed","score","vdw","elec",
                                "top_cluster_population","models_clustered"],rows)

def t7():
    rows=[]
    SETS=[("Compound4880",[("WT","Compound4880_blind"),("D184N","Compound4880_D184N_blind"),
                           ("D184A","Compound4880_D184A_blind")]),
          ("MMT",[("WT","MMT_blind"),("D184N","MMT_D184N_blind"),
                  ("D184A","MMT_D184A_blind"),("E164A","MMT_E164A_blind")])]
    for lig,rs in SETS:
        base=None
        for lab,run in rs:
            d=rank1(run)
            if not d: continue
            s,v,e,ds=(float(d[k]) for k in ("score","vdw","elec","desolv"))
            ec=W_EMREF*e
            if base is None:
                base=(s,ec,v,ds); dd=["","","",""]; pct=""
            else:
                dd=[f"{s-base[0]:.3f}",f"{ec-base[1]:.3f}",f"{v-base[2]:.3f}",f"{ds-base[3]:.3f}"]
                pct=f"{100*(s-base[0])/abs(base[0]):.1f}"
            rows.append([lig,lab,f"{s:.3f}",f"{v:.3f}",f"{e:.3f}",f"{ds:.3f}",f"{ec:.3f}",*dd,pct])
    w("T7_D184_dissection.tsv",["ligand","receptor","score","vdw","elec","desolv",
       "elec_contribution","d_score","d_elec_contribution","d_vdw","d_desolv",
       "d_score_pct_of_WT"],rows)

def t8():
    M=framed("MMT_blind"); C=framed("Compound4880_blind")
    Lm=np.vstack([[a["fit"] for a in at if a["ch"]=="B" and a["el"]!="H"] for at in M])
    Lc=np.vstack([[a["fit"] for a in at if a["ch"]=="B" and a["el"]!="H"] for at in C])
    d=np.linalg.norm(Lm[:,None,:]-Lc[None,:,:],axis=2).min(axis=1)
    def cset(F):
        c,nm=Counter(),{}
        for at in F:
            L=np.array([a["fit"] for a in at if a["ch"]=="B" and a["el"]!="H"])
            seen=set()
            for a in at:
                if a["ch"]!="A" or a["el"]=="H": continue
                if np.min(np.linalg.norm(L-a["fit"],axis=1))<=4.0:
                    seen.add(a["resi"]); nm[a["resi"]]=a["resn"]
            for x in seen: c[x]+=1
        return set(c),nm
    sm,nm1=cset(M); sc,nm2=cset(C)
    def chem(F):
        tot=ar=ac=0; n=0
        for at in F:
            n+=1
            L=np.array([a["fit"] for a in at if a["ch"]=="B" and a["el"]!="H"])
            for a in at:
                if a["ch"]!="A" or a["el"]=="H": continue
                k=int((np.linalg.norm(L-a["fit"],axis=1)<=4.5).sum())
                if k:
                    tot+=k
                    if a["resn"] in AROM: ar+=k
                    if a["resn"] in ACID and a["name"] not in BB: ac+=k
        return tot/n,100*ar/tot,100*ac/tot
    tm,am,acm=chem(M); tc,ac2,acc=chem(C)
    rows=[["MMA_atoms_within_2A_of_C4880_volume_pct",f"{100*(d<=2.0).mean():.0f}"],
          ["MMA_atom_to_C4880_cloud_mean_A",f"{d.mean():.2f}"],
          ["MMA_atom_to_C4880_cloud_max_A",f"{d.max():.2f}"],
          ["centroid_separation_A",f"{np.linalg.norm(Lm.mean(0)-Lc.mean(0)):.2f}"],
          ["contacts_MMA_only",";".join(f"{nm1[r]}{r}" for r in sorted(sm-sc)) or "(none)"],
          ["contacts_shared",";".join(f"{nm1[r]}{r}" for r in sorted(sm&sc))],
          ["contacts_C4880_only",";".join(f"{nm2[r]}{r}" for r in sorted(sc-sm))],
          ["MMA_total_contact_pairs",f"{tm:.1f}"],["MMA_aromatic_pct",f"{am:.1f}"],
          ["MMA_acidic_pct",f"{acm:.1f}"],["C4880_total_contact_pairs",f"{tc:.1f}"],
          ["C4880_aromatic_pct",f"{ac2:.1f}"],["C4880_acidic_pct",f"{acc:.1f}"]]
    w("T8_binding_mode_comparison.tsv",["measure","value"],rows)

if __name__=="__main__":
    print("data/tables/ に書き出し:")
    for fn in (t1,t2,t3,t4,t5,t6,t7,t8):
        try: fn()
        except Exception as e: print(f"  !! {fn.__name__} 失敗: {type(e).__name__}: {e}")
    print("完了")
