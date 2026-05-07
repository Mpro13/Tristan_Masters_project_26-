"""
fetch_hs_entry.py  —  AAR-COMP-012 Step 1

Fetch hydrazine synthetase information from an NCBI accession and append
to HS_database.csv in the same directory.

Usage:
    python fetch_hs_entry.py <accession> [hs_name]

Examples:
    python fetch_hs_entry.py BAW27704.1 Spb40
    python fetch_hs_entry.py WP_184911086.1 Kit2
    python fetch_hs_entry.py UUL71891.1 Afn8
"""

import sys
import csv
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = "alexedargyrou@gmail.com"

CSV_OUT = Path(__file__).parent / "HS_database.csv"

FIELDS = [
    "Hydrazine synthase name",
    "GenBank ID",
    "Accession Code",
    "UniProt ID",
    "Organism",
    "Year Deposited",
    "Reference(s)",
    "Domain Architecture",
    "Sequence Length",
    "Cupin domain range",
    "MetRS domain range",
    "Sequence",
]

CUPIN_KEYWORDS  = ("cupin", "rmlc", "barrel")
METRS_KEYWORDS  = ("metg", "methionyl", "trna", "synthetase", "ligase")


def fetch_record(accession):
    """Fetch and parse a GenBank protein record using Biopython."""
    handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def get_uniprot_id(record):
    """Extract UniProt ID from record cross-references."""
    for xref in record.dbxrefs:
        if xref.startswith("UniProtKB"):
            return xref.split(":")[-1].strip()
    return "NaN"


def get_domain_ranges(record):
    """Extract cupin and MetRS domain ranges from Region features."""
    cupin_range = "NaN"
    metrs_range = "NaN"

    for feature in record.features:
        if feature.type != "Region":
            continue

        # Collect all text from this region's qualifiers for keyword matching
        qualifiers_text = " ".join(
            " ".join(v) for v in feature.qualifiers.values()
        ).lower()

        # Convert from 0-based half-open (Biopython) to 1-based inclusive
        start = int(feature.location.start) + 1
        end   = int(feature.location.end)
        loc   = f"{start}-{end}"

        if any(k in qualifiers_text for k in CUPIN_KEYWORDS):
            cupin_range = loc
        elif any(k in qualifiers_text for k in METRS_KEYWORDS):
            metrs_range = loc

    return cupin_range, metrs_range


def get_domain_architecture(cupin_range, metrs_range):
    """Infer domain architecture from detected ranges."""
    has_cupin = cupin_range != "NaN"
    has_metrs = metrs_range != "NaN"
    if has_cupin and has_metrs:
        return "Di-domain"
    elif has_cupin:
        return "Cupin only"
    elif has_metrs:
        return "MetRS only"
    else:
        return "Unknown"


def get_references(record):
    """Return a formatted string of literature references (excluding Direct Submission)."""
    refs = []
    for ref in record.annotations.get("references", []):
        if ref.title and ref.title != "Direct Submission":
            parts = []
            if ref.authors: parts.append(ref.authors)
            parts.append(ref.title)
            if ref.journal: parts.append(ref.journal)
            refs.append(". ".join(parts))
    return " | ".join(refs) if refs else "NaN"


def build_entry(accession, hs_name=""):
    print(f"  Fetching record ... ", end="", flush=True)
    record = fetch_record(accession)
    print(f"{len(record.seq)} aa")

    cupin_range, metrs_range = get_domain_ranges(record)
    architecture = get_domain_architecture(cupin_range, metrs_range)

    return {
        "Hydrazine synthase name": hs_name or "NaN",
        "GenBank ID":              record.id,
        "Accession Code":          accession,
        "UniProt ID":              get_uniprot_id(record),
        "Organism":                record.annotations.get("organism", "NaN"),
        "Year Deposited":          record.annotations.get("date", "NaN"),
        "Reference(s)":            get_references(record),
        "Domain Architecture":     architecture,
        "Sequence Length":         len(record.seq),
        "Cupin domain range":      cupin_range,
        "MetRS domain range":      metrs_range,
        "Sequence":                str(record.seq).upper(),
    }


def print_entry(entry):
    print("\n" + "─" * 72)
    for field in FIELDS:
        val = entry.get(field, "NaN")
        if field == "Sequence" and isinstance(val, str) and val != "NaN":
            display = f"{val[:60]}...  [{len(val)} aa]"
        elif field == "Reference(s)" and isinstance(val, str) and len(val) > 80:
            display = val[:77] + "..."
        else:
            display = str(val)
        print(f"  {field:<28} {display}")
    print("─" * 72)


def already_in_csv(accession):
    if not CSV_OUT.exists():
        return False
    with open(CSV_OUT, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row.get("Accession Code") == accession or row.get("GenBank ID") == accession:
                return True
    return False


def append_to_csv(entry):
    write_header = not CSV_OUT.exists()
    with open(CSV_OUT, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(entry)
    print(f"\nAppended to: {CSV_OUT}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    accession = sys.argv[1].strip()
    hs_name   = sys.argv[2].strip() if len(sys.argv) > 2 else ""

    print(f"\nAccession: {accession}  |  Name: {hs_name or '(not provided)'}")

    if already_in_csv(accession):
        print(f"WARNING: {accession} already in {CSV_OUT.name} — skipping.")
        sys.exit(0)

    entry = build_entry(accession, hs_name)
    print_entry(entry)
    append_to_csv(entry)
