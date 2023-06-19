import numpy as np
import matplotlib.pyplot as plt

# Example data for protein sequence and reference dataset
protein_sequence = 'MSRGDWDYYVDPFWRKPCKAGMEEWYQYAGPEVDVGFYRWSR'
reference_dataset = np.random.normal(loc=5, scale=2, size=50)  # Example reference dataset with normal distribution

# Example motif occurrence in protein sequence
motif = 'Y'
observed_occurrence = protein_sequence.count(motif)

# Calculate mean and standard deviation of occurrence in reference dataset
mean_reference = np.mean(reference_dataset)
std_reference = np.std(reference_dataset)

# Calculate z-score for protein sequence
z_score = (observed_occurrence - mean_reference) / std_reference

# Plot z-score against position in protein sequence
positions = range(1, len(protein_sequence) + 1)
z_scores = [z_score for _ in range(len(protein_sequence))]
plt.plot(positions, z_scores, 'o-')
plt.xlabel('Position in Protein Sequence')
plt.ylabel('Z-score')
plt.title('Z-score Plot for Motif Occurrence')
plt.show()
