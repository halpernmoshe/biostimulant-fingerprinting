import gzip, numpy as np, pandas as pd
from pathlib import Path

BASE      = Path("C:/Users/moshe/Dropbox/ISF 2025")
COUNT_DIR = BASE / "RNA_seq_data/suppl_files/GSE297649"
SIG_DIR   = BASE / "RNA_seq_axes/signatures"
SCORE_DIR = BASE / "RNA_seq_axes/scores"

AXES = ["P_axis", "N_axis", "JA_axis", "SA_axis", "ABA_axis", "Fe_axis"]
TREATED_FILES = [
    "GSM8996165_htseqcount_WTaa_1.txt.gz",
    "GSM8996166_htseqcount_WTaa_2.txt.gz",
    "GSM8996167_htseqcount_WTaa_3.txt.gz",
]
CONTROL_FILES = [
    "GSM8996168_htseqcount_WTcontrol_1.txt.gz",
    "GSM8996169_htseqcount_WTcontrol_2.txt.gz",
    "GSM8996170_htseqcount_WTcontrol_3.txt.gz",
]
TREATED_LABELS = ["WTaa_1", "WTaa_2", "WTaa_3"]
CONTROL_LABELS = ["WTctl_1", "WTctl_2", "WTctl_3"]


def load_htseq(filepath):
    rows = {}
    with gzip.open(filepath, "rt") as fh:
        for line in fh:
            g, c = line.rstrip("\n").split("\t")
            if not g.startswith("__"):
                rows[g] = int(c)
    s = pd.Series(rows, name=str(filepath.name))
    s.index.name = "GeneID"
    return s


def to_log2cpm(s):
    return np.log2(s / s.sum() * 1e6 + 1)


def load_signature(ax):
    df  = pd.read_csv(SIG_DIR / (ax + "_logFC.csv"))
    sig = df.set_index(df.columns[0])[df.columns[1]]
    sig.index.name = "GeneID"
    return sig


def compute_score(tlist, clist, sig):
    mt = pd.concat(tlist, axis=1).mean(axis=1)
    mc = pd.concat(clist, axis=1).mean(axis=1)
    d  = mt - mc
    common = d.index.intersection(sig.index)
    return float((d[common] * sig[common]).sum()), len(common)


def main():
    SEP = "=" * 72
    print(SEP)
    print("JACKKNIFE SENSITIVITY ANALYSIS -- GSE297649 (amino acids vs control)")
    print(SEP)
    print()
    print("Loading count files...")
    tc = [load_htseq(COUNT_DIR / fn) for fn in TREATED_FILES]
    cc = [load_htseq(COUNT_DIR / fn) for fn in CONTROL_FILES]
    for i, s in enumerate(tc):
        print("  Treated  %s: %d total counts, %d genes" % (TREATED_FILES[i], s.sum(), (s>0).sum()))
    for i, s in enumerate(cc):
        print("  Control  %s: %d total counts, %d genes" % (CONTROL_FILES[i], s.sum(), (s>0).sum()))
    print()
    tl = [to_log2cpm(s) for s in tc]
    cl = [to_log2cpm(s) for s in cc]
    print("Loading axis signatures...")
    sigs = {ax: load_signature(ax) for ax in AXES}
    for ax, sig in sigs.items():
        print("  %s: %d genes" % (ax, len(sig)))
    print()
    sumdf   = pd.read_csv(SCORE_DIR / "biostimulant_6axis_summary.csv")
    aa_rows = sumdf[sumdf["treatment"] == "amino_acid"]
    kp      = dict(zip(aa_rows["axis"], aa_rows["delta_pct"]))
    print("Known full-sample %% scores (amino_acid):")
    for ax, pct in kp.items():
        print("  %s: %+.4f%%" % (ax, pct))
    print()
    print("Computing full-sample (n=3 vs n=3) scores...")
    fs = {}
    ng = {}
    for ax in AXES:
        score, n = compute_score(tl, cl, sigs[ax])
        fs[ax] = score
        ng[ax] = n
        print("  %-12s: raw = %+12.4f   shared genes = %d" % (ax, score, n))
    print()
    cal = {}
    for ax in AXES:
        if ax in kp and abs(kp[ax]) > 1e-6:
            cal[ax] = fs[ax] / (kp[ax] / 100.0)
        else:
            cal[ax] = fs[ax]
            print("  WARNING: no known %% for %s" % ax)
    print("Calibration divisors:")
    for ax in AXES:
        fp = fs[ax] / cal[ax] * 100
        print("  %-12s: divisor = %+12.1f   full%% = %+6.3f%%" % (ax, cal[ax], fp))
    print()
    results = []

    def do_jk(jt, jc, ax, lbl, jtype, nt, nc):
        jk_s, _ = compute_score(jt, jc, sigs[ax])
        fp = fs[ax]  / cal[ax] * 100
        jp = jk_s    / cal[ax] * 100
        pc = (jp - fp) / abs(fp) * 100 if abs(fp) > 1e-4 else float("nan")
        st = "YES" if abs(pc) < 20 else "NO "
        results.append({"experiment": "GSE297649_amino_acid", "jackknife_type": jtype,
            "dropped": lbl, "n_treated": nt, "n_control": nc, "axis": ax,
            "full_score_raw": fs[ax], "full_score_pct": fp,
            "jk_score_raw": jk_s, "jk_score_pct": jp,
            "pct_change_of_pct": pc, "stable_20pct": st})
        return fp, jp, pc, st

    def print_hdr():
        print("  %-12s %8s %8s %10s %7s" % ("Axis", "Full%", "JK%", "%change", "Stable"))
        print("  " + "-" * 50)

    def print_row(ax, fp, jp, pc, st):
        print("  %-12s %+8.3f %+8.3f %+10.1f %7s" % (ax, fp, jp, pc, st))

    # PART A: leave-one-treated-out
    print(SEP)
    print("PART A: Leave-one-TREATED-out  (n=2 treated vs n=3 control)")
    print(SEP)
    for di in range(3):
        kt = [j for j in range(3) if j != di]
        jt = [tl[j] for j in kt]
        kl = [TREATED_LABELS[j] for j in kt]
        print()
        print("  Drop treated %s | keep %s vs all 3 controls" % (TREATED_LABELS[di], kl))
        print_hdr()
        for ax in AXES:
            fp, jp, pc, st = do_jk(jt, cl, ax, TREATED_LABELS[di], "leave_one_treated_out", 2, 3)
            print_row(ax, fp, jp, pc, st)

    # PART B: leave-one-control-out
    print()
    print(SEP)
    print("PART B: Leave-one-CONTROL-out  (n=3 treated vs n=2 control)")
    print(SEP)
    for di in range(3):
        kc = [j for j in range(3) if j != di]
        jc = [cl[j] for j in kc]
        kl = [CONTROL_LABELS[j] for j in kc]
        print()
        print("  Drop control %s | keep all 3 treated vs %s" % (CONTROL_LABELS[di], kl))
        print_hdr()
        for ax in AXES:
            fp, jp, pc, st = do_jk(tl, jc, ax, CONTROL_LABELS[di], "leave_one_control_out", 3, 2)
            print_row(ax, fp, jp, pc, st)

    # PART C: leave-one-each-out
    print()
    print(SEP)
    print("PART C: Leave-one-EACH-out  (n=2 vs n=2, 9 combinations)")
    print(SEP)
    for dt in range(3):
        for dc in range(3):
            kt  = [j for j in range(3) if j != dt]
            kc  = [j for j in range(3) if j != dc]
            jt  = [tl[j] for j in kt]
            jc  = [cl[j] for j in kc]
            lbl = "T=%s, C=%s" % (TREATED_LABELS[dt], CONTROL_LABELS[dc])
            print()
            print("  Drop %s" % lbl)
            print_hdr()
            for ax in AXES:
                fp, jp, pc, st = do_jk(jt, jc, ax, lbl, "leave_one_each_out", 2, 2)
                print_row(ax, fp, jp, pc, st)

    # SUPPLEMENTARY: GSE138478
    print()
    print(SEP)
    print("SUPPLEMENTARY: GSE138478 (GMV & Diacetyl) -- per-sample replicate sensitivity")
    print("(n=2 per condition; leave-one-out leaves n=1 remaining)")
    print(SEP)
    gse138 = pd.read_csv(SCORE_DIR / "GSE138478_per_sample_scores_CPM.csv")
    aaxes  = gse138["axis"].unique()
    for tgrp, cgrp, ename in [
        ("GMV_treated",      "GMV_control",      "GMV_volatile"),
        ("diacetyl_treated", "diacetyl_control", "Pure_diacetyl"),
    ]:
        print()
        print("  Experiment: %s" % ename)
        for ax in aaxes:
            sub = gse138[gse138["axis"] == ax]
            ts  = sub[sub["group"] == tgrp]["cpm_pct"].values
            cs  = sub[sub["group"] == cgrp]["cpm_pct"].values
            if len(ts) == 0 or len(cs) == 0:
                continue
            fd = ts.mean() - cs.mean()
            print()
            print("    Axis: %s" % ax)
            print("    Treated per-sample%%: %s" % str(np.round(ts, 3).tolist()))
            print("    Control per-sample%%: %s" % str(np.round(cs, 3).tolist()))
            print("    Full delta (mean_T - mean_C): %+.4f%%" % fd)
            print("    %-25s %10s %10s" % ("Dropped rep", "JK delta%", "%change"))
            print("    " + "-" * 48)
            for i in range(len(ts)):
                rem = np.delete(ts, i)
                jd  = rem.mean() - cs.mean()
                pc  = (jd - fd) / abs(fd) * 100 if abs(fd) > 1e-4 else float("nan")
                print("    Drop treated rep %d:         %+10.4f %+10.1f%%" % (i+1, jd, pc))
                results.append({"experiment": "GSE138478_" + ename,
                    "jackknife_type": "leave_one_treated_out_per_sample",
                    "dropped": "treated_rep_%d" % (i+1),
                    "n_treated": len(ts)-1, "n_control": len(cs),
                    "axis": ax, "full_score_raw": float("nan"), "full_score_pct": fd,
                    "jk_score_raw": float("nan"), "jk_score_pct": jd,
                    "pct_change_of_pct": pc,
                    "stable_20pct": "YES" if abs(pc) < 20 else "NO "})
            for i in range(len(cs)):
                rem = np.delete(cs, i)
                jd  = ts.mean() - rem.mean()
                pc  = (jd - fd) / abs(fd) * 100 if abs(fd) > 1e-4 else float("nan")
                print("    Drop control rep %d:         %+10.4f %+10.1f%%" % (i+1, jd, pc))
                results.append({"experiment": "GSE138478_" + ename,
                    "jackknife_type": "leave_one_control_out_per_sample",
                    "dropped": "control_rep_%d" % (i+1),
                    "n_treated": len(ts), "n_control": len(cs)-1,
                    "axis": ax, "full_score_raw": float("nan"), "full_score_pct": fd,
                    "jk_score_raw": float("nan"), "jk_score_pct": jd,
                    "pct_change_of_pct": pc,
                    "stable_20pct": "YES" if abs(pc) < 20 else "NO "})

    # SUMMARY TABLES
    rdf = pd.DataFrame(results)
    print()
    print(SEP)
    print("SUMMARY TABLE 1: GSE297649 -- Per-axis jackknife statistics")
    print("(Across all 15 jackknife iterations per axis)")
    print(SEP)
    gr = rdf[rdf["experiment"] == "GSE297649_amino_acid"].copy()
    gr["abs_pc"] = gr["pct_change_of_pct"].abs()
    summ = gr.groupby("axis").agg(
        full_pct     = ("full_score_pct", "first"),
        jk_min_pct   = ("jk_score_pct",   "min"),
        jk_max_pct   = ("jk_score_pct",   "max"),
        max_abs_chg  = ("abs_pc",          "max"),
        mean_abs_chg = ("abs_pc",          "mean"),
    ).reset_index()
    def count_stable(x): return (x.str.strip() == "YES").sum()
    n_stable = gr.groupby("axis")["stable_20pct"].apply(count_stable).reset_index(name="n_stable")
    n_total  = gr.groupby("axis")["stable_20pct"].count().reset_index(name="n_total")
    summ = summ.merge(n_stable, on="axis").merge(n_total, on="axis")
    print()
    print("  %-12s %8s %9s %9s %9s %9s %9s" % ("Axis","Full%","JK_min%","JK_max%","MaxChg%","MeanChg%","Stable"))
    print("  " + "-" * 70)
    for _, row in summ.iterrows():
        stab = ("%d/%d" % (int(row["n_stable"]), int(row["n_total"]))).rjust(9)
        print("  %-12s %+8.3f %+9.3f %+9.3f %+9.1f %+9.1f %s" % (
            row["axis"], row["full_pct"], row["jk_min_pct"], row["jk_max_pct"],
            row["max_abs_chg"], row["mean_abs_chg"], stab))
    print()
    print(SEP)
    print("SUMMARY TABLE 2: Sign stability (did any jackknife flip the score direction?)")
    print(SEP)
    gr["full_sign"] = np.sign(gr["full_score_pct"])
    gr["jk_sign"]   = np.sign(gr["jk_score_pct"])
    gr["sign_flip"] = gr["full_sign"] != gr["jk_sign"]
    print()
    print("  %-12s %8s %7s %9s %14s" % ("Axis", "Full%", "Flips", "TotalJKs", "SignStable"))
    print("  " + "-" * 54)
    for ax in AXES:
        sub = gr[gr["axis"] == ax]
        fp  = sub["full_score_pct"].iloc[0]
        nf  = int(sub["sign_flip"].sum())
        nt  = len(sub)
        ok  = "YES" if nf == 0 else "NO (%d flips)" % nf
        print("  %-12s %+8.3f %7d %9d %14s" % (ax, fp, nf, nt, ok))
    print()
    print(SEP)
    print("SUMMARY TABLE 3: Max |%%change| by jackknife type")
    print(SEP)
    ts2  = gr.groupby(["jackknife_type", "axis"])["abs_pc"].max().reset_index()
    piv  = ts2.pivot(index="axis", columns="jackknife_type", values="abs_pc")
    cmap = {"leave_one_treated_out": "drop_1_treated",
            "leave_one_control_out": "drop_1_control",
            "leave_one_each_out":    "drop_1_each"}
    piv  = piv.rename(columns=cmap)
    def fmt_pct(x): return "%+.1f" % x
    print()
    print(piv.to_string(float_format=fmt_pct))
    print()
    print(SEP)
    print("INTERPRETATION GUIDE")
    print(SEP)
    print("  |%%change| < 20%% of full score = STABLE; sign flip = UNRELIABLE")
    print("  GSE297649 amino acid (n=3): 15 JK iterations (3+3+9 per axis)")
    print("  GSE138478 GMV/diacetyl (n=2): n=1 after dropout -- extreme sensitivity")
    print("  Amino acid scores <3%%: sign consistency is the key stability metric")
    out = SCORE_DIR / "jackknife_sensitivity_results.csv"
    rdf.to_csv(out, index=False)
    print("Full results saved to: %s" % out)
    print("Analysis complete.")


if __name__ == "__main__":
    main()
