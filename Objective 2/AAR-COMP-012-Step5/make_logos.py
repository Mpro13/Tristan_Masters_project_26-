"""
make_logos.py — AAR-COMP-012 Step 5

Reads the active site CSVs produced by Step 4 and generates three
sequence logos using information content (bits):

  logo_amp_site.png    — MetRS AMP-binding site     (11 positions, Hao 2026)
  logo_noh_site.png    — MetRS N-OH substrate site   ( 6 positions, Hao 2026)
  logo_cupin_site.png  — Cupin ZN-coordination site  (5 Å cutoff from PyrN)

Outputs are written to this directory.

Usage:
    conda activate comp_analysis
    python make_logos.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import logomaker

# ── paths ──────────────────────────────────────────────────────────────────────
STEP4_RESULTS = Path(__file__).parent.parent / "AAR-COMP-012-Step4" / "results"
OUT_DIR       = Path(__file__).parent

METRS_CSV = STEP4_RESULTS / "metrs_active_site.csv"
CUPIN_CSV = STEP4_RESULTS / "cupin_active_site.csv"

# ── MetRS active site column groups (Hao 2026) ────────────────────────────────
AMP_COLS = ["S138", "F139", "P140", "T141", "V179", "Q182",
            "A384", "L387", "R390", "N421", "R425"]
NOH_COLS = ["N143", "E284", "E286", "D420", "N421", "L424"]

# ── amino acids ───────────────────────────────────────────────────────────────
AAS = sorted("ACDEFGHIKLMNPQRSTVWY")


# ── logo function ──────────────────────────────────────────────────────────────
def make_logo(df, columns, title, out_path):
    """
    df      : DataFrame read from Step 4 CSV
    columns : ordered list of column names (positions) to include
    title   : plot title
    out_path: output PNG path
    """
    # Keep only rows that have actual amino acid data (skip gaps)
    sub = df[columns].copy()
    n   = len(sub)

    counts = {col: {aa: 0 for aa in AAS} for col in columns}
    for _, row in sub.iterrows():
        for col in columns:
            aa = str(row[col]).strip()
            if aa in AAS:
                counts[col][aa] += 1

    freq_df = pd.DataFrame(
        {col: {aa: counts[col][aa] / n for aa in AAS} for col in columns}
    ).T
    freq_df.index = range(len(columns))

    info_df = logomaker.transform_matrix(
        freq_df, from_type="probability", to_type="information"
    )

    fig, ax = plt.subplots(figsize=(max(8, len(columns) * 1.05), 3.8))
    logomaker.Logo(info_df, ax=ax, color_scheme="chemistry",
                   stack_order="small_on_top")
    ax.set_ylabel("Information (bits)", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(columns, rotation=45, ha="right", fontsize=10)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path.name}  ({n} sequences, {len(columns)} positions)")


# ── main ───────────────────────────────────────────────────────────────────────
def main():
    if not METRS_CSV.exists():
        print(f"ERROR: {METRS_CSV} not found. Run Step 4 (analyse_metrs.py) first.")
        return
    if not CUPIN_CSV.exists():
        print(f"ERROR: {CUPIN_CSV} not found. Run Step 4 (analyse_cupin.py) first.")
        return

    metrs_df = pd.read_csv(METRS_CSV)
    cupin_df = pd.read_csv(CUPIN_CSV)

    # Cupin site columns = everything after the 4 metadata columns
    cupin_cols = [c for c in cupin_df.columns
                  if c not in ("Enzyme", "RMSD", "Best_model", "Confidence")]

    print("Generating sequence logos from Step 4 results...")
    print()

    make_logo(
        metrs_df, AMP_COLS,
        "MetRS AMP-binding site — Hao 2026",
        OUT_DIR / "logo_amp_site.png",
    )
    make_logo(
        metrs_df, NOH_COLS,
        "MetRS N-OH substrate site — Hao 2026",
        OUT_DIR / "logo_noh_site.png",
    )
    make_logo(
        cupin_df, cupin_cols,
        "Cupin ZN-coordination site (5 Å cutoff, PyrN reference)",
        OUT_DIR / "logo_cupin_site.png",
    )

    print(f"\nAll logos written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
