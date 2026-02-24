#!/usr/bin/env python3
import sys, gzip
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

ROOT      = Path("C:/Users/moshe/Dropbox/ISF 2025")
SUPPL     = ROOT / "RNA_seq_data" / "suppl_files"
SIG_DIR   = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
PLOT_DIR  = ROOT / "state_space_figures"
CAL_FILE  = ROOT / "RNA_seq_axes" / "calibration.csv"
AGENT_DIR = ROOT / "agent_coordination"
