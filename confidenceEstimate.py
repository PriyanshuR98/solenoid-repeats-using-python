import numpy as np
import matplotlib.pyplot as plt

def separate_lines(theta, solenoid, globular):
    best_line = None
    best_correct = -1

    for m in np.linspace(-1, 1, num=100):  # Range of slope values from -1 to 1
        for q in np.linspace(-10, 10, num=100):  # Range of intercept values from -10 to 10
            correct_predictions = 0

            for sequence in solenoid:
                zmax, rho = sequence  # Assuming each solenoid sequence is represented by (zmax, rho)

                if m * zmax - theta * rho + q > 0:
                    correct_predictions += 1

            for sequence in globular:
                zmax, rho = sequence  # Assuming each globular sequence is represented by (zmax, rho)

                if m * zmax - theta * rho + q <= 0:
                    correct_predictions += 1

            if correct_predictions > best_correct:
                best_correct = correct_predictions
                best_line = (m, q)

    return best_line

# Example usage
solenoid_sequences = [(3, 2), (4, 5), (2, 1)]  # Example solenoid sequences
globular_sequences = [(1, 4), (2, 3), (5, 6)]  # Example globular sequences
for i in range(1,6):
    best_line = separate_lines(i, solenoid_sequences, globular_sequences)
print("Best separating line:", best_line)

# # Plotting
# zmax_values = [sequence[0] for sequence in solenoid_sequences + globular_sequences]
# rho_values = [sequence[1] for sequence in solenoid_sequences + globular_sequences]
# colors = ['r'] * len(solenoid_sequences) + ['b'] * len(globular_sequences)

# plt.scatter(zmax_values, rho_values, c=colors)
# plt.xlabel('zmax')
# plt.ylabel('ρ')
# plt.title('Solenoid and Globular Sequences')

# # Plot the best separating line
# m, q = best_line
# zmax_line = np.linspace(min(zmax_values), max(zmax_values), num=100)
# rho_line = (m * zmax_line + q) / theta
# plt.plot(zmax_line, rho_line, 'g', label='Best Separating Line')

# plt.legend()
# plt.show()
