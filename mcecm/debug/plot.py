import os
import numpy as np
import matplotlib.pyplot as plt

# Temperatures from 1.0 to 0.1 with 10 sampling points
temperatures = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
folders = [f"{i}" for i in range(1, 10)]  # Folders named "run_1" to "run_8"
#folders = ["2","4","5","8"]  # Folders named "run_1" to "run_8"
lowest_energies_by_temp = {temp: [] for temp in temperatures}  # Store energies for each temperature

# Loop through all folders and read energy files
for folder in folders:
    for temp in temperatures:
        filename = os.path.join(folder, f"total_energy_T{temp:.2f}.txt")
        try:
            # Initialize an empty list to store energy values
            energies = []
            with open(filename, 'r') as file:
                for line in file:
                    # Split the line into columns
                    parts = line.split()
                    if len(parts) > 0:  # Ensure there is at least one column
                        energy = float(parts[0])  # Read the first column as energy
                        energies.append(energy)
            # Add the lowest energy for this file to the temperature-specific list
            if energies:
                lowest_energies_by_temp[temp].append(min(energies))
        except (FileNotFoundError, ValueError) as e:
            print(f"Error reading file {filename}: {e}")

# Compute mean and std for each temperature
mean_energies = []
std_energies = []
valid_temps = []

for temp in temperatures:
    if lowest_energies_by_temp[temp]:
        valid_temps.append(temp)
        mean_energies.append(np.mean(lowest_energies_by_temp[temp]))
        std_energies.append(np.std(lowest_energies_by_temp[temp]))

# Plot the lowest energy as a function of temperature
plt.figure(figsize=(8, 6))
plt.errorbar(valid_temps, mean_energies, yerr=std_energies, fmt='o-', capsize=5, label='Mean ± Std Dev')
plt.xlabel('Temperature (T*)')
plt.ylabel('Lowest Energy')
plt.title('Lowest Energy vs. Temperature (Across Runs)')
plt.grid(True)
plt.legend()
plt.savefig('lowest_energy_vs_temperature_runs.png')
plt.show()

