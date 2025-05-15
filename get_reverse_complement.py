from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--fasta_path", type=str, required=True)
parser.add_argument("--output_path", type=str, required=True)
args = parser.parse_args()

output_records = []
for record in SeqIO.parse(args.fasta_path, "fasta"):
    rev_comp = record.seq.reverse_complement()
    record.seq = rev_comp
    output_records.append(record)

SeqIO.write(output_records, args.output_path, "fasta")
