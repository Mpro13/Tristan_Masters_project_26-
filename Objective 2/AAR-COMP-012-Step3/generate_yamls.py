"""
generate_yamls.py — AAR-COMP-012 Step 3

Reads cupin_domains.fasta and metrs_domains.fasta from Step 2,
generates Boltz2 YAML input files organised into 4 batches:

  yaml/cupin_batch1/   (sequences 1-9)
  yaml/cupin_batch2/   (sequences 10-17)
  yaml/metrs_batch1/   (sequences 1-9)
  yaml/metrs_batch2/   (sequences 10-17)

Cupin YAMLs: protein + 1 ZN cofactor.
MetRS YAMLs: protein + 2 ZN cofactors.
No affinity properties block — structure prediction only.

Usage:
    python generate_yamls.py
"""

from pathlib import Path

STEP2_DIR   = Path(__file__).parent.parent / "AAR-COMP-012-Step2"
CUPIN_FASTA = STEP2_DIR / "cupin_domains.fasta"
METRS_FASTA = STEP2_DIR / "metrs_domains.fasta"
YAML_DIR    = Path(__file__).parent / "yaml"

CUPIN_TEMPLATE = """\
sequences:
  - protein:
      id: A
      sequence: {sequence}
  - ligand:
      id: ZN1
      ccd: ZN
"""

METRS_TEMPLATE = """\
sequences:
  - protein:
      id: A
      sequence: {sequence}
  - ligand:
      id: ZN1
      ccd: ZN
  - ligand:
      id: ZN2
      ccd: ZN
"""


def parse_fasta(path):
    """Return list of (name, sequence) from a FASTA file."""
    entries = []
    name = None
    seq_lines = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                entries.append((name, "".join(seq_lines)))
            name = line[1:].split("|")[0].strip()
            seq_lines = []
        elif line.strip():
            seq_lines.append(line.strip())
    if name is not None:
        entries.append((name, "".join(seq_lines)))
    return entries


def write_yamls(entries, template, batch_size, domain_label):
    batches = [entries[i:i + batch_size] for i in range(0, len(entries), batch_size)]
    total = 0
    for batch_idx, batch in enumerate(batches, 1):
        batch_dir = YAML_DIR / f"{domain_label}_batch{batch_idx}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        for name, seq in batch:
            yaml_path = batch_dir / f"{name}_{domain_label}.yaml"
            yaml_path.write_text(template.format(sequence=seq))
            total += 1
        print(f"  {domain_label}_batch{batch_idx}: {len(batch)} YAMLs  →  {batch_dir.name}/")
        for name, _ in batch:
            print(f"    {name}_{domain_label}.yaml")
    return total


def main():
    YAML_DIR.mkdir(exist_ok=True)

    print("Parsing FASTAs...")
    cupin_entries = parse_fasta(CUPIN_FASTA)
    metrs_entries = parse_fasta(METRS_FASTA)
    print(f"  Cupin: {len(cupin_entries)} sequences")
    print(f"  MetRS: {len(metrs_entries)} sequences")

    print("\nGenerating cupin YAMLs (protein + 1 ZN)...")
    n_cupin = write_yamls(cupin_entries, CUPIN_TEMPLATE, batch_size=9, domain_label="cupin")

    print("\nGenerating MetRS YAMLs (protein + 2 ZN)...")
    n_metrs = write_yamls(metrs_entries, METRS_TEMPLATE, batch_size=9, domain_label="metrs")

    print(f"\nTotal: {n_cupin + n_metrs} YAMLs written to {YAML_DIR}/")


if __name__ == "__main__":
    main()
