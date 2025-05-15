import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser()
parser.add_argument('--scores', required=True, help='Path to tsv file containing Evo1 scored sequences')
parser.add_argument("--fasta", required=True, help="Path to fasta file containing sequence IDs")
parser.add_argument("--output", required=True, help="Path to output fasta file")
args = parser.parse_args()

# Load the scores TSV
df_scores = pd.read_csv(args.scores, sep="\t")

# Load and parse the FASTA file manually
# Ignore the first row, which contains the original sequence
fasta_records = []
with open(args.fasta, "r") as file:
    identifier = None
    sequence = []
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            if identifier:
                fasta_records.append((identifier, ''.join(sequence)))
            identifier = line[1:]
            sequence = []
        else:
            sequence.append(line)
    if identifier:
        fasta_records.append((identifier, ''.join(sequence)))

# Create a sequence-to-ID lookup
sequence_to_id = {seq: ident for ident, seq in fasta_records}

# Match fasta ID to scores
df_scores["fasta_id"] = df_scores["seqs"].map(sequence_to_id)

# Split fasta_id into site, wt_codon, mut_codon
def parse_fasta_id(fid):
    match = re.match(r"Site(\d+)_([ACGT]{3})_([ACGT]{3})", fid or "")
    if match:
        site = int(match.group(1))
        wt = match.group(2)
        mut = match.group(3)
        return pd.Series([site, wt, mut])
    return pd.Series([None, None, None])

df_scores[["site", "wt_codon", "mut_codon"]] = df_scores["fasta_id"].apply(parse_fasta_id)

# Calculate score_diff relative to the original sequence
original_score = df_scores[df_scores["fasta_id"].str.contains("_original", na=False)]["scores"].iloc[0]
df_scores["score_diff"] = df_scores["scores"] - original_score

# Select final columns and export to CSV
output_df = df_scores[["scores", "site", "wt_codon", "mut_codon", "score_diff"]]
output_df.to_csv(args.output, index=False)


