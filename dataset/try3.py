from collections import Counter

protein_sequences = ["QCTGGADCTSCTGACTGCGNCPNAVTCTNSQHCVKANTCTGSTDCNTAQTCTNSKDCFEANTCTDSTNCYKATACTNSSGCPGH"]

frequency_rank = []
sequence_count = len(protein_sequences)

for position in range(len(protein_sequences[0])):
    amino_acids_at_position = [sequence[position] for sequence in protein_sequences]
    amino_acid_counts = Counter(amino_acids_at_position)
    rank = {amino_acid: count / sequence_count for amino_acid, count in amino_acid_counts.items()}
    frequency_rank.append(rank)

print(frequency_rank)
