"""
analyse_metrs.py — AAR-COMP-012 Step 4

For each of the 17 MetRS domain Boltz2 predictions:
  1. Selects the highest-confidence model.
  2. Superimposes onto the PyrN crystal structure (PyrN.pdb) using
     Biopython Superimposer on Cα atoms from a pairwise sequence alignment.
  3. Maps Hao 2026 active site residues (AMP-binding site and N-OH
     substrate site) onto each enzyme via the alignment.

Output (results/):
  metrs_active_site.csv   — per-enzyme residue identities + RMSD

Logos are generated separately in Step 5 (make_logos.py).

Usage:
    conda activate comp_analysis
    python analyse_metrs.py
"""

import json
import csv
from pathlib import Path

from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.Align import PairwiseAligner

# ── paths ──────────────────────────────────────────────────────────────────────
STEP4_DIR   = Path(__file__).parent
STEP3_METRS = STEP4_DIR.parent / "AAR-COMP-012-Step3" / "metrs"
RESULTS_DIR = STEP4_DIR / "results"
PYRN_PDB    = STEP4_DIR / "PyrN.pdb"

# ── active site (Hao 2026) ─────────────────────────────────────────────────────
# PDB numbering = paper numbering − 119  (offset confirmed by SFPT motif search)
AMP_SITE = [
    (19, "S138"), (20, "F139"), (21, "P140"), (22, "T141"),
    (60, "V179"), (63, "Q182"),
    (265, "A384"), (268, "L387"), (271, "R390"),
    (302, "N421"), (306, "R425"),
]
NOH_SITE = [
    (24, "N143"),
    (165, "E284"), (167, "E286"),
    (301, "D420"), (302, "N421"), (305, "L424"),
]
ALL_ACTIVE_PDB = sorted({pdb for pdb, _ in AMP_SITE + NOH_SITE})

# ── amino acid helpers ─────────────────────────────────────────────────────────
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AAS = sorted(AA3TO1.values())


def get_ca_residues(chain):
    return [r for r in chain if r.id[0] == " " and "CA" in r]


def get_seq(residues):
    return "".join(AA3TO1.get(r.resname, "X") for r in residues)


# ── model selection ────────────────────────────────────────────────────────────
def find_best_model(pred_dir):
    """Return (cif_path, model_idx, confidence_score) for the best model."""
    best = (-1, None, -1)
    for i in range(3):
        name      = pred_dir.name
        conf_json = pred_dir / f"confidence_{name}_model_{i}.json"
        cif       = pred_dir / f"{name}_model_{i}.cif"
        if conf_json.exists() and cif.exists():
            score = json.loads(conf_json.read_text()).get("confidence_score", 0)
            if score > best[0]:
                best = (score, cif, i)
    return best[1], best[2], best[0]


# ── alignment + superimposition ────────────────────────────────────────────────
def align_and_superimpose(ref_residues, tgt_residues):
    """
    Global pairwise sequence alignment → extract Cα pairs →
    Biopython Superimposer → return (rmsd, [(ref_idx, tgt_idx), ...]).
    """
    aligner = PairwiseAligner()
    aligner.mode             = "global"
    aligner.match_score      = 2
    aligner.mismatch_score   = -1
    aligner.open_gap_score   = -10
    aligner.extend_gap_score = -0.5

    aln   = next(iter(aligner.align(get_seq(ref_residues), get_seq(tgt_residues))))
    pairs = []
    for (rs, re), (ts, te) in zip(aln.aligned[0], aln.aligned[1]):
        for d in range(re - rs):
            pairs.append((rs + d, ts + d))

    ref_atoms = [ref_residues[i]["CA"] for i, _ in pairs]
    tgt_atoms = [tgt_residues[j]["CA"] for _, j in pairs]
    sup = Superimposer()
    sup.set_atoms(ref_atoms, tgt_atoms)
    return sup.rms, pairs


# ── active site mapping ────────────────────────────────────────────────────────
def map_active_site(ref_residues, tgt_residues, pairs):
    """Return dict {pdb_num: single-letter aa} for all active site positions."""
    num_to_idx = {r.id[1]: i for i, r in enumerate(ref_residues)}
    ref_to_tgt = dict(pairs)
    result = {}
    for pdb_num in ALL_ACTIVE_PDB:
        ref_idx = num_to_idx.get(pdb_num)
        if ref_idx is None:
            result[pdb_num] = "-"
            continue
        tgt_idx = ref_to_tgt.get(ref_idx)
        result[pdb_num] = (
            AA3TO1.get(tgt_residues[tgt_idx].resname, "X")
            if tgt_idx is not None else "-"
        )
    return result


# ── main ───────────────────────────────────────────────────────────────────────
def main():
    RESULTS_DIR.mkdir(exist_ok=True)

    # Load PyrN crystal structure as reference
    print("Loading PyrN crystal structure...")
    ref_residues = get_ca_residues(
        PDBParser(QUIET=True).get_structure("PyrN", str(PYRN_PDB))[0]["A"]
    )
    print(f"  {len(ref_residues)} residues, chain A (res {ref_residues[0].id[1]}–{ref_residues[-1].id[1]})")

    # Add PyrN crystal as first row (self-reference, RMSD = 0)
    num_to_idx = {r.id[1]: i for i, r in enumerate(ref_residues)}
    pyrn_map   = {
        label: AA3TO1.get(ref_residues[num_to_idx[pdb]].resname, "X")
        for pdb, label in AMP_SITE + NOH_SITE if pdb in num_to_idx
    }
    results = [{"Enzyme": "PyrN_crystal", "RMSD": 0.0,
                "Best_model": "xtal", "Confidence": 1.0, **pyrn_map}]

    # Find all MetRS prediction directories
    pred_dirs = sorted(STEP3_METRS.glob("boltz_results_*/predictions/*_metrs"))
    print(f"\nFound {len(pred_dirs)} MetRS Boltz2 prediction directories")
    print(f"\n  {'Enzyme':<14} {'RMSD (Å)':>8}  {'Model':>5}  {'Conf':>6}")
    print("  " + "─" * 40)

    for pred_dir in pred_dirs:
        enzyme = pred_dir.name.replace("_metrs", "")
        best_cif, model_idx, conf = find_best_model(pred_dir)
        if best_cif is None:
            print(f"  {enzyme:<14}  no CIF found — skipping")
            continue

        tgt_residues = get_ca_residues(
            MMCIFParser(QUIET=True).get_structure("tgt", str(best_cif))[0]["A"]
        )
        rmsd, pairs  = align_and_superimpose(ref_residues, tgt_residues)
        active_map   = map_active_site(ref_residues, tgt_residues, pairs)

        label_map = {label: active_map[pdb] for pdb, label in AMP_SITE + NOH_SITE if pdb in active_map}
        print(f"  {enzyme:<14} {rmsd:>8.2f}  {model_idx:>5}  {conf:>6.3f}")
        results.append({
            "Enzyme": enzyme, "RMSD": round(rmsd, 3),
            "Best_model": model_idx, "Confidence": round(conf, 3),
            **label_map,
        })

    # Write CSV — AMP columns first, then extra N-OH columns
    amp_labels = [l for _, l in AMP_SITE]
    noh_extra  = [l for _, l in NOH_SITE if l not in amp_labels]
    fieldnames = ["Enzyme", "RMSD", "Best_model", "Confidence"] + amp_labels + noh_extra

    csv_out = RESULTS_DIR / "metrs_active_site.csv"
    with open(csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved: {csv_out.name}  ({len(results)} rows)")
    print(f"Run Step 5 (make_logos.py) to generate sequence logos.")


if __name__ == "__main__":
    main()
