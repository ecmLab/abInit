import os
import matplotlib.pyplot as plt

# Directory containing the energy files
energy_files_directory = "./"  # Update this to your directory

# Initialize data storage
energy_data = {}

# Read the energy values
file_path = "totalE.txt"
with open(file_path, 'r') as f:
    energies = [float(line.strip()) for line in f.readlines()]

# Plot the data
plt.figure(figsize=(5, 4))
mc_steps = range(1, len(energies) + 1)
plt.plot(mc_steps, energies)

# Customize the plot
plt.xlabel("MC Steps")
plt.ylabel("Total Elastic Energy")
#plt.title("Total Energy as a Function of MC Steps")
plt.grid(True)
plt.show()
