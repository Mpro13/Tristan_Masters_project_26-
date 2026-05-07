"""
analyse_cupin.py — AAR-COMP-012 Step 4

Uses the best-confidence PyrN cupin Boltz2 model as reference.
Defines the active site as all residues with any heavy atom within
5 Å of the ZN ion in that reference structure.

For each of the other 16 cupin Boltz2 predictions:
  1. Selects the highest-confidence model.
  2. Superimposes onto PyrN cupin reference using Biopython Superimposer
     on Cα atoms from a pairwise sequence alignment.
  3. Maps each ZN-proximal position onto the target enzyme.

Outputs (results/):
  cupin_active_site.csv   — per-enzyme residue identities + RMSD
  cupin_zn_site.txt       — list of PyrN active site residues detected

Logos are generated separately in Step 5 (make_logos.py).

Usage:
    conda activate comp_analysis
    python analyse_cupin.py
"""

import json
import csv
from pathlib import Path

import numpy as np

from Bio.PDB import MMCIFParser, Superimposer
from Bio.Align import PairwiseAligner

# ── paths ──────────────────────────────────────────────────────────────────────
STEP4_DIR   = Path(__file__).parent
STEP3_CUPIN = STEP4_DIR.parent / "AAR-COMP-012-Step3" / "cupin"
RESULTS_DIR = STEP4_DIR / "results"

ZN_CUTOFF = 5.0   # Å

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


# ── ZN detection ───────────────────────────────────────────────────────────────
def get_zn_coords(struct_model):
    """Return numpy array (3,) of first ZN atom found in any chain."""
    for chain in struct_model:
        for res in chain:
            if res.resname in ("ZN", "ZN2"):
                for atom in res:
                    return atom.get_vector().get_array()
    return None


def find_zn_proximal(residues, zn_coords, cutoff):
    """Return sorted list of residue indices (into `residues`) within cutoff of ZN."""
    near = []
    for i, res in enumerate(residues):
        for atom in res:
            diff = atom.get_vector().get_array() - zn_coords
            if np.linalg.norm(diff) <= cutoff:
                near.append(i)
                break
    return sorted(near)


# ── alignment + superimposition ────────────────────────────────────────────────
def align_and_superimpose(ref_residues, tgt_residues):
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
def map_active_site(tgt_residues, pairs, active_indices):
    """Map PyrN active site indices onto target via alignment pairs."""
    ref_to_tgt = dict(pairs)
    result = {}
    for ref_idx in active_indices:
        tgt_idx = ref_to_tgt.get(ref_idx)
        result[ref_idx] = (
            AA3TO1.get(tgt_residues[tgt_idx].resname, "X")
            if tgt_idx is not None else "-"
        )
    return result


# ── main ───────────────────────────────────────────────────────────────────────
def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    parser = MMCIFParser(QUIET=True)

    # ── Load PyrN cupin reference (best-confidence model) ──────────────────────
    pyrn_pred_dir = next(
        STEP3_CUPIN.glob("boltz_results_*/predictions/PyrN_cupin"), None
    )
    if pyrn_pred_dir is None:
        raise FileNotFoundError("PyrN cupin prediction directory not found.")

    ref_cif, ref_model_idx, ref_conf = find_best_model(pyrn_pred_dir)
    print(f"PyrN cupin reference: model_{ref_model_idx}  (confidence {ref_conf:.3f})")

    ref_struct   = parser.get_structure("ref", str(ref_cif))
    ref_residues = get_ca_residues(ref_struct[0]["A"])
    print(f"  {len(ref_residues)} residues")

    # ── Find ZN-proximal active site residues in PyrN reference ───────────────
    zn_coords = get_zn_coords(ref_struct[0])
    if zn_coords is None:
        raise ValueError("No ZN found in PyrN cupin reference CIF.")

    active_indices = find_zn_proximal(ref_residues, zn_coords, ZN_CUTOFF)
    print(f"\nZN-proximal residues in PyrN cupin ({ZN_CUTOFF} Å cutoff): "
          f"{len(active_indices)} positions")

    # Build position labels: {aa}{1-based position in cupin sequence}
    labels = [
        f"{AA3TO1.get(ref_residues[i].resname,'X')}{ref_residues[i].id[1]}"
        for i in active_indices
    ]
    print(f"  Positions: {', '.join(labels)}")

    # Write site summary
    site_txt = RESULTS_DIR / "cupin_zn_site.txt"
    site_txt.write_text(
        f"PyrN cupin ZN-proximal active site ({ZN_CUTOFF} Å cutoff)\n"
        f"Reference: model_{ref_model_idx} (confidence {ref_conf:.3f})\n\n"
        + "\n".join(
            f"  idx {i:3d}  {AA3TO1.get(ref_residues[i].resname,'X')}"
            f"  res {ref_residues[i].id[1]}"
            for i in active_indices
        )
    )
    print(f"  Saved: {site_txt.name}")

    # PyrN self-row (residues from reference directly)
    pyrn_map = {i: AA3TO1.get(ref_residues[i].resname, "X") for i in active_indices}
    results  = [{"Enzyme": "PyrN_ref", "RMSD": 0.0,
                 "Best_model": ref_model_idx, "Confidence": ref_conf, **pyrn_map}]

    # ── Process all other cupin enzymes ───────────────────────────────────────
    all_pred_dirs = sorted(STEP3_CUPIN.glob("boltz_results_*/predictions/*_cupin"))
    other_dirs    = [d for d in all_pred_dirs if d.name != "PyrN_cupin"]
    print(f"\nFound {len(other_dirs)} other cupin prediction directories")
    print(f"\n  {'Enzyme':<14} {'RMSD (Å)':>8}  {'Model':>5}  {'Conf':>6}")
    print("  " + "─" * 40)

    for pred_dir in other_dirs:
        enzyme = pred_dir.name.replace("_cupin", "")
        best_cif, model_idx, conf = find_best_model(pred_dir)
        if best_cif is None:
            print(f"  {enzyme:<14}  no CIF found — skipping")
            continue

        tgt_residues = get_ca_residues(
            parser.get_structure("tgt", str(best_cif))[0]["A"]
        )
        rmsd, pairs  = align_and_superimpose(ref_residues, tgt_residues)
        active_map   = map_active_site(tgt_residues, pairs, active_indices)

        print(f"  {enzyme:<14} {rmsd:>8.2f}  {model_idx:>5}  {conf:>6.3f}")
        results.append({
            "Enzyme": enzyme, "RMSD": round(rmsd, 3),
            "Best_model": model_idx, "Confidence": round(conf, 3),
            **active_map,
        })

    # ── Write CSV ──────────────────────────────────────────────────────────────
    fieldnames = ["Enzyme", "RMSD", "Best_model", "Confidence"] + list(active_indices)
    csv_out = RESULTS_DIR / "cupin_active_site.csv"
    with open(csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        # Write with label header
        fh.write(",".join(["Enzyme", "RMSD", "Best_model", "Confidence"] + labels) + "\n")
        for row in results:
            fh.write(",".join([
                str(row["Enzyme"]), str(row["RMSD"]),
                str(row["Best_model"]), str(row["Confidence"]),
            ] + [str(row.get(i, "-")) for i in active_indices]) + "\n")
    print(f"\nSaved: {csv_out.name}  ({len(results)} rows)")
    print(f"Run Step 5 (make_logos.py) to generate sequence logos.")


if __name__ == "__main__":
    main()
