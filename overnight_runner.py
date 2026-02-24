#!/usr/bin/env python3
"""
overnight_runner.py  —  Master controller for overnight pipeline
================================================================
Runs in sequence:
  1. geo_download.py          — extended search + download all target datasets
  2. build_axes_rnaseq.R      — build axis signatures + score biostimulants

Logs progress to overnight_run_log.txt
"""

import subprocess, sys, time, os
from pathlib import Path
from datetime import datetime

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT    = Path(__file__).parent
LOG_F   = ROOT / "overnight_run_log.txt"

def ts():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def log(msg):
    line = f"[{ts()}] {msg}"
    print(line)
    with open(LOG_F, "a", encoding="utf-8") as f:
        f.write(line + "\n")

def run_step(name, cmd, cwd=None, timeout_sec=7200):
    log(f"=== START: {name} ===")
    log(f"Command: {' '.join(str(c) for c in cmd)}")
    t0 = time.time()
    try:
        proc = subprocess.Popen(
            cmd, cwd=str(cwd or ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace"
        )
        step_log = ROOT / f"step_{name.replace(' ','_')}.log"
        with open(step_log, "w", encoding="utf-8") as lf:
            for line in proc.stdout:
                sys.stdout.write(line)
                lf.write(line)
        proc.wait(timeout=timeout_sec)
        elapsed = time.time() - t0
        if proc.returncode == 0:
            log(f"=== OK: {name} ({elapsed:.0f}s) ===")
            return True
        else:
            log(f"=== FAIL: {name} exit={proc.returncode} ({elapsed:.0f}s) ===")
            return False
    except subprocess.TimeoutExpired:
        proc.kill()
        log(f"=== TIMEOUT: {name} (>{timeout_sec}s) ===")
        return False
    except Exception as e:
        log(f"=== ERROR: {name}: {e} ===")
        return False

def find_rscript():
    """Find Rscript.exe on Windows."""
    candidates = [
        r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.3.3\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.3.2\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.3.1\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.4.0\bin\Rscript.exe",
    ]
    # Also try PATH
    import shutil
    from_path = shutil.which("Rscript")
    if from_path:
        return from_path
    # Glob R versions
    import glob
    for pattern in [r"C:\Program Files\R\R-4*\bin\Rscript.exe",
                    r"C:\Program Files\R\R-4*\bin\x64\Rscript.exe"]:
        hits = sorted(glob.glob(pattern), reverse=True)
        if hits:
            return hits[0]
    for c in candidates:
        if os.path.exists(c):
            return c
    return "Rscript"  # hope it's on PATH

def main():
    log("=" * 60)
    log("OVERNIGHT PIPELINE STARTED")
    log(f"Root: {ROOT}")
    log("=" * 60)

    python = sys.executable
    rscript = find_rscript()
    log(f"Python: {python}")
    log(f"Rscript: {rscript}")

    results = {}

    # ── Phase 1+2: Download ──────────────────────────────────────
    ok = run_step(
        "download",
        [python, str(ROOT / "Metaanalysis" / "geo_download.py")],
        cwd=ROOT,
        timeout_sec=7200   # 2 hours max
    )
    results["download"] = ok

    # ── Phase 3+4: Build axes + score ────────────────────────────
    ok = run_step(
        "build_axes",
        [rscript, "--vanilla", str(ROOT / "build_axes_rnaseq.R")],
        cwd=ROOT,
        timeout_sec=7200
    )
    results["build_axes"] = ok

    # ── Summary ────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("OVERNIGHT PIPELINE COMPLETE")
    log("=" * 60)
    for step, ok in results.items():
        log(f"  {step}: {'SUCCESS' if ok else 'FAILED'}")

    # Check what was produced
    log("\nOutput files produced:")
    check_files = [
        ROOT / "all_axes_scores_rnaseq.csv",
        ROOT / "RNA_seq_axes" / "signatures",
        ROOT / "RNA_seq_axes" / "scores",
        ROOT / "state_space_figures",
        ROOT / "RNA_seq_data" / "download_log.csv",
        ROOT / "Metaanalysis" / "geo_biostim_candidates_extended.csv",
    ]
    for f in check_files:
        if f.is_file():
            size = f.stat().st_size
            log(f"  [OK]   {f.name}  ({size//1024} KB)")
        elif f.is_dir():
            n = len(list(f.iterdir()))
            log(f"  [DIR]  {f.name}/  ({n} files)")
        else:
            log(f"  [MISS] {f.name}")

    log(f"\nFull log: {LOG_F}")
    log("Good morning!")

if __name__ == "__main__":
    main()
