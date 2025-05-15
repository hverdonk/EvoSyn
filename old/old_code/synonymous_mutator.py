from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import argparse
import sys
import re

# === INPUTS ===
parser = argparse.ArgumentParser()
parser.add_argument("--fasta_path", required=True, help="Path to input fasta")
parser.add_argument("--species_id", required=True, help="Fasta sequence id of the sequence to mutate")
parser.add_argument("--output_path", required=True, help="Path to output fasta file")
args = parser.parse_args()


# === CODON TABLE ===
# Standard table, can be adjusted if using a different genetic code
codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
amino_acid_to_codons = {}
for codon, aa in codon_table.forward_table.items():
    if aa not in amino_acid_to_codons:
        amino_acid_to_codons[aa] = []
    amino_acid_to_codons[aa].append(codon)

# === READ INPUT ===
species_seq = None
cleaned_species_id = None
for record in SeqIO.parse(args.fasta_path, "fasta"):
    if args.species_id in record.id:
        raw_seq = str(record.seq).upper()
        cleaned_seq = re.sub(r'[^ACGT]', '', raw_seq)
        if raw_seq != cleaned_seq:
            print(f"[Warning] Sequence for {args.species_id} contained gaps or ambiguous characters and was cleaned.")
        species_seq = cleaned_seq
        cleaned_species_id = record.id.replace(" ", "_")
        break

if species_seq is None:
    raise ValueError(f"Species ID '{args.species_id}' not found in alignment, or contained only gaps.")

# Ensure length is divisible by 3
if len(species_seq) % 3 != 0:
    raise ValueError("Sequence length is not divisible by 3 (not in-frame).")

# === EXTRACT CODONS ===
codons = [species_seq[i:i+3] for i in range(0, len(species_seq), 3)]

# Check for in-frame stop codons
for idx, codon in enumerate(codons):
    if codon in codon_table.stop_codons:
        raise ValueError(f"In-frame stop codon ({codon}) found at site {idx + 1}. Sequence processing aborted.")

# Validate start codon
if codons[0] != "ATG":
    raise ValueError(f"First codon is not Methionine (ATG): {codons[0]}")

# === CREATE OUTPUT RECORDS ===
output_records = []

# Add original sequence
output_records.append(SeqIO.SeqRecord(Seq("".join(codons)), id=cleaned_species_id + "_original", description=""))

# === PROCESS ALL SITES ===
for site_index, original_codon in enumerate(codons):
    aa = codon_table.forward_table.get(original_codon)
    if aa in ("M", "W"):
        print(f"[Info] Skipping site {site_index + 1}: codon {original_codon} encodes {aa} which has no synonymous codons.")
        continue
    if aa is None:
        print(f"[Warning] Skipping site {site_index + 1}: codon {original_codon} could not be translated.")
        continue  # not a valid coding codon, already filtered stop codons

    for alt_codon in amino_acid_to_codons[aa]:
        if alt_codon != original_codon:
            new_codons = codons.copy()
            new_codons[site_index] = alt_codon
            new_seq = Seq("".join(new_codons))
            output_records.append(
                SeqIO.SeqRecord(
                    new_seq,
                    id=f"Site{site_index + 1}_{original_codon}_{alt_codon}",
                    description=""
                )
            )

# === WRITE OUTPUT ===
SeqIO.write(output_records, args.output_path, "fasta")
print(f"Written {len(output_records)} sequences to {args.output_path}")
