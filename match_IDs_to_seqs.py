import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--scores', type=int, default=250, help='Length of flanking region in base pairs')
parser.add_argument("--fasta", required=True, help="Path to fasta file containing sequence IDs")
parser.add_argument("--species_id", required=True, help="Fasta sequence id of the sequence to mutate")
parser.add_argument("--output_path", required=True, help="Path to output fasta file")
args = parser.parse_args()

# Create a lookup dictionary from sequences to identifiers
sequence_to_id = {seq: identifier for identifier, seq in fasta_records}

# Match each sequence in the TSV to its identifier from the FASTA
df_scores["fasta_id"] = df_scores["seqs"].map(sequence_to_id)
