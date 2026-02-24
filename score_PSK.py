import csv, json, math, os, openpyxl

ROOT = r"""C:/Users/moshe/Dropbox/ISF 2025"""
SIG_DIR = os.path.join(ROOT, "RNA_seq_axes", "signatures")
SCORES_DIR = os.path.join(ROOT, "RNA_seq_axes", "scores")
CALIB_FILE = os.path.join(ROOT, "RNA_seq_axes", "calibration.csv")

AXES = ["P_axis","N_axis","ABA_axis","JA_axis","SA_axis",
        "Fe_axis","CK_axis","ET_axis","SA_axis_v2","Auxin_axis"]

def load_sig(an):
    fname = os.path.join(SIG_DIR, an+"_logFC.csv")
    sig = {}
    with open(fname) as fh:
        for r in csv.DictReader(fh):
            g = r.get("GeneID") or r.get("gene_id") or r.get("gene")
            v = r.get("logFC")  or r.get("log2FoldChange") or r.get("logfc")
            if g and v:
                try: sig[g.strip()] = float(v)
                except ValueError: pass
    return sig

def load_calib():
    c = {}
    with open(CALIB_FILE) as fh:
        for r in csv.DictReader(fh):
            c[r["axis"].strip()] = (float(r["ctrl_score"]), float(r["treat_score"]))
    return c

def score_samp(counts, sig):
    return sum(math.log2(counts.get(g,0)+1)*v for g,v in sig.items())

def delta(cs, ts, sig):
    cv = [score_samp(s,sig) for s in cs]
    tv = [score_samp(s,sig) for s in ts]
    return sum(tv)/len(tv)-sum(cv)/len(cv), cv, tv

def pct(d,c0,c1):
    r=c1-c0; return None if r==0 else d/r*100.0

def mk_samps(rows, idx):
    s=[{} for _ in idx]
    for row in rows[1:]:
        g=row[0]
        if not g or not str(g).startswith("AT"): continue
        g=str(g).strip()
        for k,ci in enumerate(idx):
            v=row[ci] if ci<len(row) else None
            s[k][g]=int(v) if v is not None else 0
    return s

print("Loading signatures...")
sigs={}
for ax in AXES:
    try:
        sigs[ax]=load_sig(ax)
        print("  "+ax+": "+str(len(sigs[ax]))+" genes")
    except FileNotFoundError as e:
        print("  "+ax+": MISSING "+str(e))

calib=load_calib()
print("Calibration:",list(calib.keys()))

print("Loading GSE254987...")
wb=openpyxl.load_workbook(os.path.join(ROOT,"RNA_seq_data","suppl_files","GSE254987","GSE254987_gene_count.xlsx"),read_only=True)
r987=list(wb.active.iter_rows(values_only=True)); wb.close()
print("Cols:",r987[0][:13])
c987=mk_samps(r987,[4,5,6]); p100=mk_samps(r987,[7,8,9]); p1uM=mk_samps(r987,[10,11,12])
print("Genes:",len(c987[0]))

print("Loading GSE254986...")
wb=openpyxl.load_workbook(os.path.join(ROOT,"RNA_seq_data","suppl_files","GSE254986","GSE254986_count_gene.xlsx"),read_only=True)
r986=list(wb.active.iter_rows(values_only=True)); wb.close()
print("All cols:",r986[0])
csh=mk_samps(r986,[4,5,6]); psh5=mk_samps(r986,[13,14,15])
crt=mk_samps(r986,[22,23,24]); prt5=mk_samps(r986,[31,32,33])
print("Shoot:",len(csh[0])," Root:",len(crt[0]))

comps={
    "GSE254987_PSK_1uM_5h_seedling":  {"c":c987,"t":p1uM, "desc":"GSE254987 tpst seedlings 1uM PSK 5h",  "tis":"seedling","dose":"1 uM", "tp":"5h"},
    "GSE254987_PSK_100nM_5h_seedling":{"c":c987,"t":p100, "desc":"GSE254987 tpst seedlings 100nM PSK 5h","tis":"seedling","dose":"100 nM","tp":"5h"},
    "GSE254986_PSK_10nM_5h_shoot":    {"c":csh, "t":psh5, "desc":"GSE254986 tpst shoot 10nM PSK 5h",    "tis":"shoot",   "dose":"10 nM", "tp":"5h"},
    "GSE254986_PSK_10nM_5h_root":     {"c":crt, "t":prt5, "desc":"GSE254986 tpst root 10nM PSK 5h",     "tis":"root",    "dose":"10 nM", "tp":"5h"},
}

print("Scoring...")
res={}
for cn,co in comps.items():
    print("  "+cn)
    asc={}
    for ax in AXES:
        if ax not in sigs: continue
        sig=sigs[ax]
        d_val,cv,tv=delta(co["c"],co["t"],sig)
        pc=pct(d_val,calib[ax][0],calib[ax][1]) if ax in calib else None
        ov=sum(1 for g in sig if g in co["c"][0])
        asc[ax]={"delta":round(d_val,2),"pct_calibration":round(pc,2) if pc is not None else None,
                 "ctrl_scores":[round(s,2) for s in cv],"treat_scores":[round(s,2) for s in tv],
                 "n_overlap_genes":ov}
        tag=(("+" if pc>=0 else "")+str(round(pc,1))+"%") if pc is not None else "N/A"
        print("    "+ax.ljust(15)+"  d="+str(round(d_val,1))+"  "+tag)
    res[cn]={"metadata":co["desc"],"tissue":co["tis"],"dose":co["dose"],"timepoint":co["tp"],"axes":asc}

outp=os.path.join(SCORES_DIR,"PSK_10axis_scores.json")
with open(outp,"w") as fh: json.dump(res,fh,indent=2)
print("Saved: "+outp)

print("")
print("SUMMARY TABLE (% calibration):")
ck=list(res.keys()); sh=["987_1uM","987_100nM","986_shoot","986_root"]
hd="Axis".ljust(18)+"".join(x.rjust(14) for x in sh)
print(hd); print("-"*len(hd))
for ax in AXES:
    vals=[]
    for k in ck:
        v=res[k]["axes"].get(ax,{}).get("pct_calibration")
        vals.append(((("+" if v>=0 else "")+str(round(v,1))+"%").rjust(14)) if v is not None else "N/A".rjust(14))
    print(ax.ljust(18)+"".join(vals))
print("Done.")
