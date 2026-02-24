import sys, gzip, tarfile, time, traceback, re
from datetime import datetime
from pathlib import Path
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import requests, numpy as np, pandas as pd
ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
SUPPL_DIR  = ROOT / "RNA_seq_data" / "suppl_files"
MATRIX_DIR = ROOT / "RNA_seq_data" / "series_matrices"
SIG_DIR    = ROOT / "RNA_seq_axes" / "signatures"
SCORES_DIR = ROOT / "RNA_seq_axes" / "scores"
LOG_FILE   = ROOT / "RNA_seq_axes" / "extend_pipeline_log.txt"
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
for d in [SUPPL_DIR, MATRIX_DIR, SIG_DIR, SCORES_DIR]: d.mkdir(parents=True, exist_ok=True)
log_fh = open(LOG_FILE, "w", encoding="utf-8", errors="replace")
def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    log_fh.write(line + chr(10))
    log_fh.flush()
log("=" * 70)
log("Biostimulant Scoring Pipeline -- Extended Run (v2)")
log("=" * 70)
def geo_prefix(acc): return acc[: len(acc) - 3] + "nnn"
def download_file(url, dest_path, retries=3, chunk=1 << 20):
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, stream=True, timeout=120)
            if r.status_code == 200:
                with open(dest_path, "wb") as f:
                    for block in r.iter_content(chunk): f.write(block)
                log(f"  Downloaded: {dest_path.name} ({dest_path.stat().st_size/1e6:.1f} MB)")
                return True
            else: log(f"  HTTP {r.status_code} for {url} (attempt {attempt})")
        except Exception as e_dl: log(f"  Error: {e_dl} (attempt {attempt})")
        time.sleep(5)
    return False
def get_filelist(acc):
    prefix = geo_prefix(acc)
    url = f"{GEO_FTP_BASE}/{prefix}/{acc}/suppl/filelist.txt"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200: return r.text
    except Exception as e_fl: log(f"  Filelist error {acc}: {e_fl}")
    return None
TARGET_GSES = ["GSE208223", "GSE138478", "GSE199501"]
log(chr(10) + "-" * 60)
log("STEP 1 -- Download series matrices")
log("-" * 60)
for acc in TARGET_GSES:
    prefix = geo_prefix(acc)
    dest2 = MATRIX_DIR / f"{acc}_series_matrix.txt.gz"
    if dest2.exists(): log(f"{acc}: series matrix already present."); continue
    url = f"{GEO_FTP_BASE}/{prefix}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    log(f"{acc}: downloading series matrix...")
    ok = download_file(url, dest2)
    if not ok: log(f"{acc}: WARNING -- series matrix download failed.")
log(chr(10) + "-" * 60)
log("STEP 2 -- Download count / supplementary files")
log("-" * 60)
def download_suppl(acc):
    out_dir = SUPPL_DIR / acc
    out_dir.mkdir(parents=True, exist_ok=True)
    filelist_text = get_filelist(acc)
    if filelist_text is None: log(f"  {acc}: no filelist.txt found."); return
    lines_fl = [l.strip() for l in filelist_text.splitlines() if l.strip()]
    files_to_get = []
    for line in lines_fl[1:]:
        parts = line.split(chr(9))
        if not parts: continue
        fname = parts[0].strip()
        if re.search(r"(RAW[.]tar|counts?[.]gz|count.*[.]gz|_count[.]gz)", fname, re.I): files_to_get.append(fname)
    if not files_to_get:
        for line in lines_fl[1:]:
            parts = line.split(chr(9))
            if parts and (parts[0].strip().endswith(".gz") or parts[0].strip().endswith(".tar")): files_to_get.append(parts[0].strip())
    if not files_to_get: log(f"  {acc}: no suppl files."); log(f"  {filelist_text[:400]}"); return
    pfx = geo_prefix(acc)
    for fname in files_to_get:
        d2 = out_dir / fname
        if d2.exists(): log(f"  {acc}/{fname}: already present."); continue
        u2 = f"{GEO_FTP_BASE}/{pfx}/{acc}/suppl/{fname}"
        log(f"  {acc}: downloading {fname}")
        download_file(u2, d2)
for acc in TARGET_GSES: log(f"{chr(10)}{acc}: fetching supplementary files..."); download_suppl(acc)
log(chr(10) + "-" * 60)
log("STEP 3 -- Search GEO for K-starvation RNA-seq datasets")
log("-" * 60)
ESEARCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
K_SEARCH_TERMS = ["Arabidopsis potassium starvation RNA-seq",
    "Arabidopsis thaliana potassium deficiency transcriptome",
    "Triticum aestivum potassium deficiency transcriptome",
    "potassium starvation Arabidopsis RNA-seq high-throughput sequencing"]
k_gse_hits = {}
def esearch_gds(term, retmax=20):
    params = dict(db="gds", term=term, retmax=retmax, retmode="json", email="moshe@research.com")
    try:
        r = requests.get(ESEARCH_URL, params=params, timeout=30)
        if r.status_code == 200:
            data = r.json()
            ids = data.get("esearchresult", {}).get("idlist", [])
            count = data.get("esearchresult", {}).get("count", "?")
            log(f"  Term: {repr(term)} => {count} results, {len(ids)} fetched")
            return ids
    except Exception as e_es: log(f"  esearch error: {e_es}")
    return []
def esummary_gds(uid_list):
    if not uid_list: return {}
    params = dict(db="gds", id=",".join(uid_list), retmode="json", email="moshe@research.com")
    try:
        r = requests.get(ESUMMARY_URL, params=params, timeout=30)
        if r.status_code == 200: return r.json().get("result", {})
    except Exception as e_sm: log(f"  esummary error: {e_sm}")
    return {}
all_k_uids = set()
for term in K_SEARCH_TERMS: all_k_uids.update(esearch_gds(term)); time.sleep(0.4)
log(f"{chr(10)}Total K-starvation GDS UIDs: {len(all_k_uids)}")
if all_k_uids:
    summaries = esummary_gds(list(all_k_uids))
    for uid, rec in summaries.items():
        if uid == "uids": continue
        a_k = rec.get("accession", ""); t_k = rec.get("title", "")
        g_k = rec.get("gdstype", ""); x_k = rec.get("taxon", "")
        n_k = rec.get("n_samples", "?")
        log(f"  {a_k} | {g_k} | {x_k} | n={n_k} | {t_k[:70]}")
        if a_k.startswith("GSE"): k_gse_hits[a_k] = rec
log(f"{chr(10)}K-starvation GSE entries: {list(k_gse_hits.keys())}")
log(chr(10) + "-" * 60)
log("STEP 4 -- Download K-starvation count files")
log("-" * 60)
K_ACCS_TO_TRY = list(k_gse_hits.keys())[:6]
for acc in K_ACCS_TO_TRY:
    log(f"{chr(10)}{acc}: fetching K-starvation suppl...")
    pfx = geo_prefix(acc)
    md = MATRIX_DIR / f"{acc}_series_matrix.txt.gz"
    if not md.exists(): download_file(f"{GEO_FTP_BASE}/{pfx}/{acc}/matrix/{acc}_series_matrix.txt.gz", md)
    download_suppl(acc)
log(chr(10) + "-" * 60)
log("STEP 5 -- Extract / parse count matrices")
log("-" * 60)
def extract_tar_for_acc(acc):
    out_dir = SUPPL_DIR / acc
    for f in list(out_dir.glob("*.tar")):
        log(f"  {acc}: extracting {f.name}...")
        try:
            with tarfile.open(f) as tf: tf.extractall(out_dir)
            log(f"  {acc}: extraction complete.")
        except Exception as e2: log(f"  {acc}: tar error: {e2}")
SKIP_COLS = {"chr", "start", "end", "strand", "length"}
def read_count_file(path):
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            df = pd.read_csv(fh, sep=chr(9), comment="#", index_col=0, header=0)
        df = df[~df.index.astype(str).str.startswith("__")]
        num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        count_cols = [c for c in num_cols if c.lower() not in SKIP_COLS]
        if not count_cols: return None
        sub = df[count_cols].copy()
        sub.columns = [c.replace(".bam", "").replace(".BAM", "") for c in sub.columns]
        return sub
    except Exception as er: log(f"    Cannot parse {path.name}: {er}"); return None
def read_xls_count_file(path):
    import tempfile, os
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rb") as fh:
            data = fh.read()
        tmp = Path(tempfile.gettempdir()) / (Path(path).stem + "_tmp.xls")
        with open(tmp, "wb") as f:
            f.write(data)
        df = pd.read_excel(str(tmp), index_col=0)
        try: os.remove(str(tmp))
        except: pass
        num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        skip_kw2 = {"cpm", "logfc", "log2fc", "log2cpm", "log2", "pval", "p-value", "fdr", "padj"}
        count_cols = [c for c in num_cols if " " not in str(c) and str(c).lower() not in skip_kw2 and not str(c).lower().startswith("log")]
        if not count_cols: return None
        sub = df[count_cols].fillna(0)
        if sub.min().min() >= 0: sub = sub.astype(int)
        return sub
    except Exception as ex:
        log(f"  xls read error {Path(path).name}: {ex}")
        return None

def load_counts_for_acc(acc):
    out_dir = SUPPL_DIR / acc
    if not out_dir.exists(): return None
    extract_tar_for_acc(acc)
    cands = []
    xls_cands = list(out_dir.rglob("*.xls.gz")) + list(out_dir.rglob("*.xls"))
    for ext in ["*.count.gz","*count*.gz","*counts*.gz","*.htseq*","*featureCounts*",
                "*.txt.gz","*.tsv.gz","*.csv.gz","*.count","*.counts","*.txt","*.tsv"]:
        cands.extend(out_dir.rglob(ext))
    skip_kw = ["series_matrix", "readme", "filelist", "md5"]
    cands = [c for c in cands if not any(kw in c.name.lower() for kw in skip_kw)]
    seen = set(); cands = [c for c in cands if str(c) not in seen and not seen.add(str(c))]
    if not cands and not xls_cands: log(f"  {acc}: no candidate count files."); return None
    if cands: log(f"  {acc}: {len(cands)} candidate count file(s).")
    dfs = []
    for c in cands[:50]:
        sub = read_count_file(c)
        if sub is not None and len(sub) > 100:
            stem = c.stem.split(".")[0]
            sub.columns = [f"{stem}__{col}" for col in sub.columns]
            dfs.append(sub)
    for xc in xls_cands:
        if str(xc) not in seen and not seen.add(str(xc)):
            xsub = read_xls_count_file(xc)
            if xsub is not None and len(xsub) > 100:
                stem2 = xc.stem.split(".")[0]
                xsub.columns = [f"{stem2}__{col}" for col in xsub.columns]
                dfs.append(xsub)
    if not dfs: log(f"  {acc}: no parseable count files."); return None
    counts = pd.concat(dfs, axis=1)
    if counts.min().min() >= 0:
        counts = counts.fillna(0).astype(int)
    else:
        counts = counts.fillna(0)
    log(f"  {acc}: counts matrix {counts.shape} (genes x samples)")
    return counts
log(chr(10) + "-" * 60)
log("STEP 6 -- Build K-starvation axis signature (pydeseq2)")
log("-" * 60)
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
K_SIG_PATH = SIG_DIR / "K_axis_logFC.csv"
def run_deseq2_sig(counts_df, cond_map, ref="control", aname="K"):
    samps = [s for s in counts_df.columns if s in cond_map]
    if len(samps) < 4: log(f"  Not enough labeled samples ({len(samps)}), skipping."); return None
    ct = counts_df[samps].T
    ct = ct.loc[:, ct.sum() > 0].astype(int)
    meta = pd.DataFrame({"condition": [cond_map[s] for s in samps]}, index=samps)
    nt = int((meta["condition"] != ref).sum())
    nc = int((meta["condition"] == ref).sum())
    log(f"  DESeq2: {ct.shape[1]} genes, {nc} ctrl + {nt} treat")
    if nc < 2 or nt < 2: log("  Need >=2 replicates per group."); return None
    try:
        dds = DeseqDataSet(counts=ct, metadata=meta, design_factors="condition",
                           ref_level=["condition", ref], refit_cooks=True)
        dds.deseq2()
        ds = DeseqStats(dds, contrast=["condition", "treatment", ref])
        ds.summary()
        res = ds.results_df[["log2FoldChange"]].rename(columns={"log2FoldChange": "logFC"})
        res = res.dropna()
        log(f"  {aname}: {len(res)} genes with logFC.")
        return res
    except Exception as ed:
        log(f"  DESeq2 failed: {ed}")
        traceback.print_exc()
        return None
def infer_cond_map(counts_df, acc):
    mat_path = MATRIX_DIR / f"{acc}_series_matrix.txt.gz"
    if not mat_path.exists(): return None
    try:
        with gzip.open(mat_path, "rt", encoding="utf-8", errors="replace") as fh:
            mlines = fh.readlines()
    except Exception: return None
    sline = tline = ""
    for ml in mlines:
        if ml.startswith("!Sample_geo_accession"): sline = ml
        if ml.startswith("!Sample_title"): tline = ml
    gsm_ids = re.findall(r"GSM[0-9]+", sline)
    titles = [f.strip(chr(34)).strip() for f in tline.split(chr(9))[1:] if f.strip(chr(34)).strip()]
    if not gsm_ids or not titles: return None
    cmap = {}
    for gsm, title in zip(gsm_ids, titles):
        tl = title.lower()
        if any(kw in tl for kw in ["starv", "deficien", "low k", "low-k", "without k", "k0", "potassium-free"]):
            cmap[gsm] = "treatment"
        elif any(kw in tl for kw in ["suffic", "control", "ctrl", "normal", "full", "replete"]):
            cmap[gsm] = "control"
        else: cmap[gsm] = "treatment" if "k" in tl else "control"
    mapped = {}
    for col in counts_df.columns:
        hit = False
        for gsm, cond in cmap.items():
            if gsm in col or col in gsm: mapped[col] = cond; hit = True; break
        if not hit:
            tl = col.lower()
            if any(kw in tl for kw in ["starv", "deficien", "low", "0mm", "k0"]): mapped[col] = "treatment"
            elif any(kw in tl for kw in ["control", "ctrl", "suff", "full", "replete"]): mapped[col] = "control"
    if not mapped:
        cols = list(counts_df.columns); half = len(cols)//2
        for c in cols[:half]: mapped[c] = "control"
        for c in cols[half:]: mapped[c] = "treatment"
        log(f"  {acc}: fallback condition split")
    nc2 = sum(1 for v in mapped.values() if v=="control")
    nt2 = sum(1 for v in mapped.values() if v=="treatment")
    log(f"  {acc}: {nc2} ctrl, {nt2} treat")
    return mapped
K_ALL_COUNTS = []
for acc in K_ACCS_TO_TRY:
    log(f"{chr(10)}  Loading counts for {acc}...")
    ct2 = load_counts_for_acc(acc)
    if ct2 is not None: K_ALL_COUNTS.append((acc, ct2))
K_sig_df = None
if K_SIG_PATH.exists():
    log(f"{chr(10)}K-axis exists, loading...")
    K_sig_df = pd.read_csv(K_SIG_PATH, index_col=0)
    log(f"  {len(K_sig_df)} genes.")
elif K_ALL_COUNTS:
    log(f"{chr(10)}Building K-axis from {len(K_ALL_COUNTS)} dataset(s)...")
    lfc_parts = []
    for acc, ct3 in K_ALL_COUNTS:
        log(f"{chr(10)}  DESeq2 on {acc}...")
        cm = infer_cond_map(ct3, acc)
        if cm is None: log(f"  {acc}: cannot infer conditions."); continue
        r2 = run_deseq2_sig(ct3, cm, aname=f"K({acc})")
        if r2 is not None: lfc_parts.append(r2)
    if lfc_parts:
        mg = pd.concat(lfc_parts, axis=1)
        mg.columns = [f"logFC_{i}" for i in range(len(mg.columns))]
        K_sig_df = mg.mean(axis=1).to_frame("logFC")
        K_sig_df.to_csv(K_SIG_PATH)
        log(f"{chr(10)}K-axis saved: {K_SIG_PATH} ({len(K_sig_df)} genes)")
    else: log(f"{chr(10)}No K-axis results -- signature not built.")
else: log(f"{chr(10)}No K data -- K-axis not built.")

log(chr(10) + "-" * 60)
log("STEP 7 -- Load existing axis signatures")
log("-" * 60)

SIGNATURES = {}
sig_files = {
    "P": SIG_DIR / "P_axis_logFC.csv",
    "N": SIG_DIR / "N_axis_logFC.csv",
    "ABA": SIG_DIR / "ABA_axis_logFC.csv",
    "K": SIG_DIR / "K_axis_logFC.csv",
}
for axis, spath in sig_files.items():
    if spath.exists():
        sdf = pd.read_csv(spath, index_col=0)
        SIGNATURES[axis] = sdf["logFC"]
        log(f"  {axis}-axis: {len(sdf)} genes from {spath.name}")
    else:
        log(f"  {axis}-axis: not found ({spath.name}), skipping.")
log(f"{chr(10)}Axes available: {list(SIGNATURES.keys())}")

log(chr(10) + "-" * 60)
log("STEP 8 -- Score biostimulant datasets on all axes")
log("-" * 60)

def score_on_axis(log2counts_s, logfc_s):
    common = log2counts_s.index.intersection(logfc_s.index)
    if len(common) < 50: return float("nan")
    return float((log2counts_s[common] * logfc_s[common]).sum())

BIO_ACCS = ["GSE208223", "GSE138478", "GSE199501"]
all_scores = []

for bacc in BIO_ACCS:
    log(f"{chr(10)}Scoring {bacc}...")
    cts = load_counts_for_acc(bacc)
    if cts is None:
        log(f"  {bacc}: no count data, skipping.")
        continue
    log(f"  {bacc}: {cts.shape[1]} samples, {cts.shape[0]} genes")
    for col in cts.columns:
        col_data = cts[col]
        if isinstance(col_data, pd.DataFrame): col_data = col_data.iloc[:, 0]
        raw = col_data.astype(float)
        if raw.min() < 0:
            log2c = raw
        else:
            log2c = np.log2(raw + 1)
        log2c.index = cts.index.astype(str).str.split(".").str[0]
        row = {"dataset": bacc, "sample": col}
        for axis, logfc_s in SIGNATURES.items():
            logfc_s.index = logfc_s.index.astype(str)
            sc = score_on_axis(log2c, logfc_s)
            row[f"{axis}_score"] = sc
            log(f"    {col} | {axis}: {sc:.2f}" if not (sc != sc) else f"    {col} | {axis}: nan")
        all_scores.append(row)

if all_scores:
    scores_df = pd.DataFrame(all_scores)
    score_cols = [c for c in scores_df.columns if c.endswith("_score")]
    out_csv = SCORES_DIR / "biostimulant_extended_scores.csv"
    scores_df.to_csv(out_csv, index=False)
    log(f"{chr(10)}Scores saved: {out_csv} ({len(scores_df)} rows)")

    summary_rows = []
    for bacc in scores_df["dataset"].unique():
        sub = scores_df[scores_df["dataset"] == bacc]
        srow = {"dataset": bacc, "n_samples": len(sub)}
        for sc in score_cols:
            srow[f"{sc}_mean"] = sub[sc].mean()
            srow[f"{sc}_std"] = sub[sc].std()
        summary_rows.append(srow)
    summary_df = pd.DataFrame(summary_rows)
    sum_csv = SCORES_DIR / "biostimulant_extended_summary.csv"
    summary_df.to_csv(sum_csv, index=False)
    log(f"Summary saved: {sum_csv}")
    log(f"{chr(10)}=== SUMMARY ===")
    log(summary_df.to_string())
else:
    log(f"{chr(10)}No scores computed.")

log(chr(10) + "=" * 60)
log("extend_pipeline.py complete.")
log("=" * 60)
log_fh.close()
log('\n' + '-' * 60)
log('STEP 6 -- Build K-starvation axis signature (pydeseq2)')
log('-' * 60)

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

K_SIG_PATH = SIG_DIR / 'K_axis_logFC.csv'

def run_deseq2_sig(counts_df, cond_map, ref='control', aname='K'):
    samps = [s for s in counts_df.columns if s in cond_map]
    if len(samps) < 4:
        log(f'  Not enough labeled samples ({len(samps)}), skipping.')
        return None
    ct = counts_df[samps].T
    ct = ct.loc[:, ct.sum() > 0].astype(int)
    meta = pd.DataFrame({'condition': [cond_map[s] for s in samps]}, index=samps)
    nt = int((meta['condition'] != ref).sum())
    nc = int((meta['condition'] == ref).sum())
    log(f'  DESeq2: {ct.shape[1]} genes, {nc} ctrl + {nt} treat')
    if nc < 2 or nt < 2:
        log('  Need >=2 replicates per group.')
        return None
    try:
        dds = DeseqDataSet(counts=ct, metadata=meta,
                           design_factors='condition',
                           ref_level=['condition', ref], refit_cooks=True)
        dds.deseq2()
        ds = DeseqStats(dds, contrast=['condition', 'treatment', ref])
        ds.summary()
        res = ds.results_df[['log2FoldChange']].rename(
            columns={'log2FoldChange': 'logFC'})
        res = res.dropna()
        log(f'  {aname}: {len(res)} genes with logFC.')
        return res
    except Exception as ed:
        log(f'  DESeq2 failed: {ed}')
        traceback.print_exc()
        return None

def infer_cond_map(counts_df, acc):
    mat_path = MATRIX_DIR / f'{acc}_series_matrix.txt.gz'
    if not mat_path.exists(): return None
    try:
        with gzip.open(mat_path, 'rt', encoding='utf-8', errors='replace') as fh:
            mlines = fh.readlines()
    except Exception: return None
    sline = tline = ''
    for ml in mlines:
        if ml.startswith('!Sample_geo_accession'): sline = ml
        if ml.startswith('!Sample_title'): tline = ml
    gsm_ids = re.findall(r'GSM\d+', sline)
    titles = [f.strip(chr(34)).strip() for f in tline.split(chr(9))[1:] if f.strip(chr(34)).strip()]
    if not gsm_ids or not titles: return None
    cmap = {}
    for gsm, title in zip(gsm_ids, titles):
        tl = title.lower()
        if any(kw in tl for kw in ['starv','deficien','low k',
                                    'low-k','0 mm k','without k',
                                    'k0','no k','potassium-free']):
            cmap[gsm] = 'treatment'
        elif any(kw in tl for kw in ['suffic','control','ctrl',
                                      'normal','full','replete']):
            cmap[gsm] = 'control'
        else:
            cmap[gsm] = 'treatment' if 'k' in tl else 'control'
    mapped = {}
    for col in counts_df.columns:
        hit = False
        for gsm, cond in cmap.items():
            if gsm in col or col in gsm:
                mapped[col] = cond; hit = True; break
        if not hit:
            tl = col.lower()
            if any(kw in tl for kw in ['starv','deficien','low','0mm','k0']):
                mapped[col] = 'treatment'
            elif any(kw in tl for kw in ['control','ctrl','suff','full','replete']):
                mapped[col] = 'control'
    if not mapped:
        cols = list(counts_df.columns)
        half = len(cols)//2
        for c in cols[:half]: mapped[c] = 'control'
        for c in cols[half:]: mapped[c] = 'treatment'
        log(f'  {acc}: fallback condition split')
    nc2 = sum(1 for v in mapped.values() if v=='control')
    nt2 = sum(1 for v in mapped.values() if v=='treatment')
    log(f'  {acc}: {nc2} ctrl, {nt2} treat')
    return mapped

K_ALL_COUNTS = []
for acc in K_ACCS_TO_TRY:
    log(f'\n  Loading counts for {acc}...')
    ct2 = load_counts_for_acc(acc)
    if ct2 is not None: K_ALL_COUNTS.append((acc, ct2))

K_sig_df = None
if K_SIG_PATH.exists():
    log(f'\nK-axis exists, loading from {K_SIG_PATH}...')
    K_sig_df = pd.read_csv(K_SIG_PATH, index_col=0)
    log(f'  {len(K_sig_df)} genes.')
elif K_ALL_COUNTS:
    log(f'\nBuilding K-axis from {len(K_ALL_COUNTS)} dataset(s)...')
    lfc_parts = []
    for acc, ct3 in K_ALL_COUNTS:
        log(f'\n  DESeq2 on {acc}...')
        cm = infer_cond_map(ct3, acc)
        if cm is None:
            log(f'  {acc}: cannot infer conditions.')
            continue
        r2 = run_deseq2_sig(ct3, cm, aname=f'K({acc})')
        if r2 is not None: lfc_parts.append(r2)
    if lfc_parts:
        mg = pd.concat(lfc_parts, axis=1)
        mg.columns = [f'logFC_{i}' for i in range(len(mg.columns))]
        K_sig_df = mg.mean(axis=1).to_frame('logFC')
        K_sig_df.to_csv(K_SIG_PATH)
        log(f'\nK-axis saved: {K_SIG_PATH} ({len(K_sig_df)} genes)')
    else:
        log('\nNo K-axis results -- signature not built.')
else:
    log('\nNo K data -- K-axis not built.')
log('\n' + '-' * 60)
log('STEP 6 -- Build K-starvation axis signature (pydeseq2)')
log('-' * 60)

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

K_SIG_PATH = SIG_DIR / 'K_axis_logFC.csv'

def run_deseq2_sig(counts_df, cond_map, ref='control', aname='K'):
    samps = [s for s in counts_df.columns if s in cond_map]
    if len(samps) < 4:
        log(f'  Not enough labeled samples ({len(samps)}), skipping.')
        return None
    ct = counts_df[samps].T
    ct = ct.loc[:, ct.sum() > 0].astype(int)
    meta = pd.DataFrame({'condition': [cond_map[s] for s in samps]}, index=samps)
    nt = int((meta['condition'] != ref).sum())
    nc = int((meta['condition'] == ref).sum())
    log(f'  DESeq2: {ct.shape[1]} genes, {nc} ctrl + {nt} treat')
    if nc < 2 or nt < 2:
        log('  Need >=2 replicates per group.')
        return None
    try:
        dds = DeseqDataSet(counts=ct, metadata=meta,
                           design_factors='condition',
                           ref_level=['condition', ref], refit_cooks=True)
        dds.deseq2()
        ds = DeseqStats(dds, contrast=['condition', 'treatment', ref])
        ds.summary()
        res = ds.results_df[['log2FoldChange']].rename(
            columns={'log2FoldChange': 'logFC'})
        res = res.dropna()
        log(f'  {aname}: {len(res)} genes with logFC.')
        return res
    except Exception as ed:
        log(f'  DESeq2 failed: {ed}')
        traceback.print_exc()
        return None

def infer_cond_map(counts_df, acc):
    mat_path = MATRIX_DIR / f'{acc}_series_matrix.txt.gz'
    if not mat_path.exists(): return None
    try:
        with gzip.open(mat_path, 'rt', encoding='utf-8', errors='replace') as fh:
            mlines = fh.readlines()
    except Exception: return None
    sline = tline = ''
    for ml in mlines:
        if ml.startswith('!Sample_geo_accession'): sline = ml
        if ml.startswith('!Sample_title'): tline = ml
    gsm_ids = re.findall(r'GSM\d+', sline)
    titles = [f.strip(chr(34)).strip() for f in tline.split(chr(9))[1:] if f.strip(chr(34)).strip()]
    if not gsm_ids or not titles: return None
    cmap = {}
    for gsm, title in zip(gsm_ids, titles):
        tl = title.lower()
        if any(kw in tl for kw in ['starv','deficien','low k',
                                    'low-k','0 mm k','without k',
                                    'k0','no k','potassium-free']):
            cmap[gsm] = 'treatment'
        elif any(kw in tl for kw in ['suffic','control','ctrl',
                                      'normal','full','replete']):
            cmap[gsm] = 'control'
        else:
            cmap[gsm] = 'treatment' if 'k' in tl else 'control'
    mapped = {}
    for col in counts_df.columns:
        hit = False
        for gsm, cond in cmap.items():
            if gsm in col or col in gsm:
                mapped[col] = cond; hit = True; break
        if not hit:
            tl = col.lower()
            if any(kw in tl for kw in ['starv','deficien','low','0mm','k0']):
                mapped[col] = 'treatment'
            elif any(kw in tl for kw in ['control','ctrl','suff','full','replete']):
                mapped[col] = 'control'
    if not mapped:
        cols = list(counts_df.columns)
        half = len(cols)//2
        for c in cols[:half]: mapped[c] = 'control'
        for c in cols[half:]: mapped[c] = 'treatment'
        log(f'  {acc}: fallback condition split')
    nc2 = sum(1 for v in mapped.values() if v=='control')
    nt2 = sum(1 for v in mapped.values() if v=='treatment')
    log(f'  {acc}: {nc2} ctrl, {nt2} treat')
    return mapped

K_ALL_COUNTS = []
for acc in K_ACCS_TO_TRY:
    log(f'\n  Loading counts for {acc}...')
    ct2 = load_counts_for_acc(acc)
    if ct2 is not None: K_ALL_COUNTS.append((acc, ct2))

K_sig_df = None
if K_SIG_PATH.exists():
    log(f'\nK-axis exists, loading from {K_SIG_PATH}...')
    K_sig_df = pd.read_csv(K_SIG_PATH, index_col=0)
    log(f'  {len(K_sig_df)} genes.')
elif K_ALL_COUNTS:
    log(f'\nBuilding K-axis from {len(K_ALL_COUNTS)} dataset(s)...')
    lfc_parts = []
    for acc, ct3 in K_ALL_COUNTS:
        log(f'\n  DESeq2 on {acc}...')
        cm = infer_cond_map(ct3, acc)
        if cm is None:
            log(f'  {acc}: cannot infer conditions.')
            continue
        r2 = run_deseq2_sig(ct3, cm, aname=f'K({acc})')
        if r2 is not None: lfc_parts.append(r2)
    if lfc_parts:
        mg = pd.concat(lfc_parts, axis=1)
        mg.columns = [f'logFC_{i}' for i in range(len(mg.columns))]
        K_sig_df = mg.mean(axis=1).to_frame('logFC')
        K_sig_df.to_csv(K_SIG_PATH)
        log(f'\nK-axis saved: {K_SIG_PATH} ({len(K_sig_df)} genes)')
    else:
        log('\nNo K-axis results -- signature not built.')
else:
    log('\nNo K data -- K-axis not built.')
log('\n' + '-' * 60)
log('STEP 6 -- Build K-starvation axis signature (pydeseq2)')
log('-' * 60)

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

K_SIG_PATH = SIG_DIR / 'K_axis_logFC.csv'

def run_deseq2_sig(counts_df, cond_map, ref='control', aname='K'):
    samps = [s for s in counts_df.columns if s in cond_map]
    if len(samps) < 4:
        log(f'  Not enough labeled samples ({len(samps)}), skipping.')
        return None
    ct = counts_df[samps].T
    ct = ct.loc[:, ct.sum() > 0].astype(int)
    meta = pd.DataFrame({'condition': [cond_map[s] for s in samps]}, index=samps)
    nt = int((meta['condition'] != ref).sum())
    nc = int((meta['condition'] == ref).sum())
    log(f'  DESeq2: {ct.shape[1]} genes, {nc} ctrl + {nt} treat')
    if nc < 2 or nt < 2:
        log('  Need >=2 replicates per group.')
        return None
    try:
        dds = DeseqDataSet(counts=ct, metadata=meta,
                           design_factors='condition',
                           ref_level=['condition', ref], refit_cooks=True)
        dds.deseq2()
        ds = DeseqStats(dds, contrast=['condition', 'treatment', ref])
        ds.summary()
        res = ds.results_df[['log2FoldChange']].rename(
            columns={'log2FoldChange': 'logFC'})
        res = res.dropna()
        log(f'  {aname}: {len(res)} genes with logFC.')
        return res
    except Exception as ed:
        log(f'  DESeq2 failed: {ed}')
        traceback.print_exc()
        return None

def infer_cond_map(counts_df, acc):
    mat_path = MATRIX_DIR / f'{acc}_series_matrix.txt.gz'
    if not mat_path.exists(): return None
    try:
        with gzip.open(mat_path, 'rt', encoding='utf-8', errors='replace') as fh:
            mlines = fh.readlines()
    except Exception: return None
    sline = tline = ''
    for ml in mlines:
        if ml.startswith('!Sample_geo_accession'): sline = ml
        if ml.startswith('!Sample_title'): tline = ml
    gsm_ids = re.findall(r'GSM\d+', sline)
    titles = [f.strip(chr(34)).strip() for f in tline.split(chr(9))[1:] if f.strip(chr(34)).strip()]
    if not gsm_ids or not titles: return None
    cmap = {}
    for gsm, title in zip(gsm_ids, titles):
        tl = title.lower()
        if any(kw in tl for kw in ['starv','deficien','low k',
                                    'low-k','0 mm k','without k',
                                    'k0','no k','potassium-free']):
            cmap[gsm] = 'treatment'
        elif any(kw in tl for kw in ['suffic','control','ctrl',
                                      'normal','full','replete']):
            cmap[gsm] = 'control'
        else:
            cmap[gsm] = 'treatment' if 'k' in tl else 'control'
    mapped = {}
    for col in counts_df.columns:
        hit = False
        for gsm, cond in cmap.items():
            if gsm in col or col in gsm:
                mapped[col] = cond; hit = True; break
        if not hit:
            tl = col.lower()
            if any(kw in tl for kw in ['starv','deficien','low','0mm','k0']):
                mapped[col] = 'treatment'
            elif any(kw in tl for kw in ['control','ctrl','suff','full','replete']):
                mapped[col] = 'control'
    if not mapped:
        cols = list(counts_df.columns)
        half = len(cols)//2
        for c in cols[:half]: mapped[c] = 'control'
        for c in cols[half:]: mapped[c] = 'treatment'
        log(f'  {acc}: fallback condition split')
    nc2 = sum(1 for v in mapped.values() if v=='control')
    nt2 = sum(1 for v in mapped.values() if v=='treatment')
    log(f'  {acc}: {nc2} ctrl, {nt2} treat')
    return mapped

K_ALL_COUNTS = []
for acc in K_ACCS_TO_TRY:
    log(f'\n  Loading counts for {acc}...')
    ct2 = load_counts_for_acc(acc)
    if ct2 is not None: K_ALL_COUNTS.append((acc, ct2))

K_sig_df = None
if K_SIG_PATH.exists():
    log(f'\nK-axis exists, loading from {K_SIG_PATH}...')
    K_sig_df = pd.read_csv(K_SIG_PATH, index_col=0)
    log(f'  {len(K_sig_df)} genes.')
elif K_ALL_COUNTS:
    log(f'\nBuilding K-axis from {len(K_ALL_COUNTS)} dataset(s)...')
    lfc_parts = []
    for acc, ct3 in K_ALL_COUNTS:
        log(f'\n  DESeq2 on {acc}...')
        cm = infer_cond_map(ct3, acc)
        if cm is None:
            log(f'  {acc}: cannot infer conditions.')
            continue
        r2 = run_deseq2_sig(ct3, cm, aname=f'K({acc})')
        if r2 is not None: lfc_parts.append(r2)
    if lfc_parts:
        mg = pd.concat(lfc_parts, axis=1)
        mg.columns = [f'logFC_{i}' for i in range(len(mg.columns))]
        K_sig_df = mg.mean(axis=1).to_frame('logFC')
        K_sig_df.to_csv(K_SIG_PATH)
        log(f'\nK-axis saved: {K_SIG_PATH} ({len(K_sig_df)} genes)')
    else:
        log('\nNo K-axis results -- signature not built.')
else:
    log('\nNo K data -- K-axis not built.')
