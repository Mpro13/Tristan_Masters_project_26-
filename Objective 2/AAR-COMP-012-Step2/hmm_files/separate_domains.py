"""
separate_domains.py  —  AAR-COMP-012 Step 2

Reads HS_database.csv, extracts cupin and MetRS domain sequences, and
writes two FASTA files. For entries where domain ranges are unknown,
runs hmmscan against Pfam HMMs to detect boundaries automatically.
Updates HS_database.csv with any newly detected ranges.

Usage:
    python separate_domains.py

Outputs (written to the same directory):
    cupin_domains.fasta
    metrs_domains.fasta

Requires:
    - HMMER loaded (module load HMMER/3.4-gompi-2023a)
    - Pfam HMMs downloaded to ./hmm/ (script downloads them if missing)
    - HS_database.csv from Step 1
"""

import csv
import subprocess
import urllib.request
import shutil
import sys
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
STEP1_CSV  = Path(__file__).parent.parent / "AAR-COMP-012-Step1" / "HS_database.csv"
OUT_DIR    = Path(__file__).parent
HMM_DIR    = OUT_DIR / "hmm"
CUPIN_HMM  = HMM_DIR / "cupin.hmm"
METRS_HMM  = HMM_DIR / "metrs.hmm"
CUPIN_FASTA = OUT_DIR / "cupin_domains.fasta"
METRS_FASTA = OUT_DIR / "metrs_domains.fasta"

# Pfam accessions
CUPIN_PFAM  = "PF00190"   # cupin_1
METRS_PFAM  = "PF09334"   # MetRS_core
METRS_PFAM2 = "PF00133"   # tRNA-synt_1 (class I catalytic domain — broader coverage)

CSV_FIELDS = [
    "Hydrazine synthase name", "GenBank ID", "Accession Code", "UniProt ID",
    "Organism", "Year Deposited", "Reference(s)", "Domain Architecture",
    "Sequence Length", "Cupin domain range", "MetRS domain range", "Sequence",
]


# ── HMM download ──────────────────────────────────────────────────────────────

def download_hmm(pfam_id, out_path):
    """Download a single Pfam HMM from InterPro."""
    url = f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{pfam_id}/?annotation=hmm"
    print(f"  Downloading {pfam_id} from InterPro ... ", end="", flush=True)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read()
        # Response is gzip-compressed
        import gzip, io
        with gzip.open(io.BytesIO(data)) as gz:
            hmm_text = gz.read().decode()
        out_path.write_text(hmm_text)
        print("done")
    except Exception as e:
        print(f"FAILED — {e}")
        sys.exit(1)


def ensure_hmms():
    HMM_DIR.mkdir(exist_ok=True)

    if not CUPIN_HMM.exists():
        download_hmm(CUPIN_PFAM, CUPIN_HMM)

    if not METRS_HMM.exists():
        # Concatenate both MetRS HMMs for broader coverage
        tmp1 = HMM_DIR / "metrs_core.hmm"
        tmp2 = HMM_DIR / "trna_synt1.hmm"
        download_hmm(METRS_PFAM,  tmp1)
        download_hmm(METRS_PFAM2, tmp2)
        METRS_HMM.write_text(tmp1.read_text() + tmp2.read_text())
        tmp1.unlink()
        tmp2.unlink()

    # Build combined pressed database (required by hmmscan)
    combined_hmm = HMM_DIR / "combined.hmm"
    combined_h3i = HMM_DIR / "combined.hmm.h3i"
    if not combined_h3i.exists():
        combined_hmm.write_text(CUPIN_HMM.read_text() + METRS_HMM.read_text())
        print(f"  Pressing HMM database ... ", end="", flush=True)
        result = subprocess.run(
            ["hmmpress", str(combined_hmm)],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"FAILED\n{result.stderr}")
            sys.exit(1)
        print("done")


# ── Hmmer ─────────────────────────────────────────────────────────────────────

def run_hmmscan(name, sequence):
    """
    Run hmmscan on a single sequence against both cupin and MetRS HMMs.
    Returns dict: {"cupin": "start-end" or None, "metrs": "start-end" or None}
    """
    tmp_fa  = OUT_DIR / f"_tmp_{name}.fasta"
    tmp_out = OUT_DIR / f"_tmp_{name}.hmmscan"

    tmp_fa.write_text(f">{name}\n{sequence}\n")

    combined_hmm = HMM_DIR / "combined.hmm"

    result = subprocess.run(
        ["hmmscan", "--domtblout", str(tmp_out),
         "--noali", "-E", "1e-3",
         str(combined_hmm), str(tmp_fa)],
        capture_output=True, text=True
    )

    hits = {"cupin": None, "metrs": None}

    if tmp_out.exists():
        for line in tmp_out.read_text().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 20:
                continue
            target_name = cols[0].lower()   # HMM name
            ali_from    = int(cols[17])     # alignment start on sequence (1-based)
            ali_to      = int(cols[18])     # alignment end on sequence (1-based)
            loc = f"{ali_from}-{ali_to}"
            if "cupin" in target_name or CUPIN_PFAM.lower() in target_name:
                hits["cupin"] = loc
            elif "metrs" in target_name or "mets" in target_name or METRS_PFAM.lower() in target_name:
                hits["metrs"] = loc

    # Clean up temp files
    tmp_fa.unlink(missing_ok=True)
    tmp_out.unlink(missing_ok=True)

    return hits


# ── domain extraction ─────────────────────────────────────────────────────────

def extract_domain(sequence, range_str):
    """Slice sequence using '28-101' range string (1-based inclusive)."""
    start, end = map(int, range_str.split("-"))
    return sequence[start - 1:end]


def resolve_entry(row):
    """
    Determine cupin and MetRS domain sequences for one CSV row.
    Returns (cupin_seq, metrs_seq, updated_cupin_range, updated_metrs_range, method)
    """
    seq        = row["Sequence"]
    arch       = row["Domain Architecture"]
    cupin_rng  = row["Cupin domain range"]
    metrs_rng  = row["MetRS domain range"]
    name       = row["Hydrazine synthase name"]

    cupin_seq  = None
    metrs_seq  = None
    method     = "annotated"

    # ── case 1: both ranges present (di-domain with annotations) ──────────────
    if cupin_rng != "NaN" and metrs_rng != "NaN":
        cupin_seq = extract_domain(seq, cupin_rng)
        metrs_seq = extract_domain(seq, metrs_rng)

    # ── case 2: cupin range present, no MetRS (cupin-only protein) ─────────────
    elif cupin_rng != "NaN" and metrs_rng == "NaN":
        cupin_seq = extract_domain(seq, cupin_rng)

    # ── case 3: MetRS range present, no cupin (MetRS-only protein) ─────────────
    elif metrs_rng != "NaN" and cupin_rng == "NaN":
        metrs_seq = extract_domain(seq, metrs_rng)

    # ── case 4: both NaN — run Hmmer ───────────────────────────────────────────
    elif cupin_rng == "NaN" and metrs_rng == "NaN":
        print(f"  [{name}] both ranges NaN — running Hmmer ...", end="", flush=True)
        hits = run_hmmscan(name, seq)
        method = "hmmer"

        if hits["cupin"]:
            cupin_rng = hits["cupin"]
            cupin_seq = extract_domain(seq, cupin_rng)
            print(f" cupin={cupin_rng}", end="")
        if hits["metrs"]:
            metrs_rng = hits["metrs"]
            metrs_seq = extract_domain(seq, metrs_rng)
            print(f" metrs={metrs_rng}", end="")

        # Hmmer found nothing — infer from sequence length as last resort
        if not hits["cupin"] and not hits["metrs"]:
            seq_len = len(seq)
            if seq_len <= 200:
                cupin_rng = f"1-{seq_len}"
                cupin_seq = seq
                arch      = "Cupin only"
                method    = "inferred (length)"
                print(f" no Hmmer hit — inferred cupin (full seq, {seq_len} aa)", end="")
            elif seq_len >= 300:
                metrs_rng = f"1-{seq_len}"
                metrs_seq = seq
                arch      = "MetRS only"
                method    = "inferred (length)"
                print(f" no Hmmer hit — inferred MetRS (full seq, {seq_len} aa)", end="")
            else:
                method = "failed"
                print(f" no hits — flagged for manual review", end="")
        print()

        # Update architecture based on Hmmer result (only if not already inferred)
        if method == "hmmer":
            if hits["cupin"] and hits["metrs"]:
                arch = "Di-domain"
            elif hits["cupin"]:
                arch = "Cupin only"
            elif hits["metrs"]:
                arch = "MetRS only"

    return cupin_seq, metrs_seq, cupin_rng, metrs_rng, arch, method


# ── FASTA writing ─────────────────────────────────────────────────────────────

def fasta_header(row, domain, range_str):
    """Format: >Name|Organism|range"""
    name     = row["Hydrazine synthase name"]
    organism = row["Organism"]
    return f">{name}|{organism}|{range_str}"


def write_fasta(path, entries):
    """entries: list of (header, sequence)"""
    with open(path, "w") as fh:
        for header, seq in entries:
            fh.write(f"{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if not STEP1_CSV.exists():
        print(f"ERROR: {STEP1_CSV} not found. Run Step 1 first.")
        sys.exit(1)

    if not shutil.which("hmmscan"):
        print("ERROR: hmmscan not found in PATH. Load HMMER first:")
        print("  module load HMMER/3.4-gompi-2023a")
        sys.exit(1)

    rows = list(csv.DictReader(open(STEP1_CSV)))
    print(f"Loaded {len(rows)} entries from {STEP1_CSV.name}")

    # Download Pfam HMMs if needed
    unknown = [r for r in rows if r["Domain Architecture"] == "Unknown"]
    if unknown:
        print(f"\n{len(unknown)} Unknown entries — downloading Pfam HMMs if needed...")
        ensure_hmms()

    # Process all entries
    print("\nExtracting domains...")
    cupin_entries = []
    metrs_entries = []
    summary       = []
    updated_rows  = []

    for row in rows:
        name = row["Hydrazine synthase name"]
        cupin_seq, metrs_seq, cupin_rng, metrs_rng, arch, method = resolve_entry(row)

        # Update row with any new ranges / architecture
        row["Cupin domain range"]  = cupin_rng
        row["MetRS domain range"]  = metrs_rng
        row["Domain Architecture"] = arch
        updated_rows.append(row)

        # Collect FASTA entries
        if cupin_seq:
            header = fasta_header(row, "cupin", cupin_rng)
            cupin_entries.append((header, cupin_seq))

        if metrs_seq:
            header = fasta_header(row, "metrs", metrs_rng)
            metrs_entries.append((header, metrs_seq))

        summary.append({
            "name":      name,
            "arch":      arch,
            "cupin_rng": cupin_rng,
            "metrs_rng": metrs_rng,
            "cupin_len": len(cupin_seq) if cupin_seq else 0,
            "metrs_len": len(metrs_seq) if metrs_seq else 0,
            "method":    method,
        })

    # Write FASTAs
    write_fasta(CUPIN_FASTA, cupin_entries)
    write_fasta(METRS_FASTA, metrs_entries)
    print(f"\nWrote: {CUPIN_FASTA.name}  ({len(cupin_entries)} sequences)")
    print(f"Wrote: {METRS_FASTA.name}  ({len(metrs_entries)} sequences)")

    # Update CSV
    with open(STEP1_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(updated_rows)
    print(f"Updated: {STEP1_CSV.name}")

    # Print summary table
    print("\n" + "─" * 90)
    print(f"  {'Name':<22} {'Architecture':<12} {'Cupin range':<14} {'MetRS range':<14} {'Method'}")
    print("─" * 90)
    flagged = []
    for s in summary:
        method_tag = "" if s["method"] == "annotated" else f"[{s['method']}]"
        print(f"  {s['name']:<22} {s['arch']:<12} {s['cupin_rng']:<14} {s['metrs_rng']:<14} {method_tag}")
        if s["method"] == "failed":
            flagged.append(s["name"])
    print("─" * 90)
    print(f"\nTotal cupin sequences : {len(cupin_entries)}")
    print(f"Total MetRS sequences : {len(metrs_entries)}")
    if flagged:
        print(f"\nFLAGGED (no Hmmer hit, not written to FASTA): {', '.join(flagged)}")
    else:
        print("\nAll entries resolved successfully.")


if __name__ == "__main__":
    main()
