import numpy as np
import matplotlib.pyplot as plt

# Temperatures from 1.0 to 0.1 with 10 sampling points
temperatures = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
lowest_energies = []

for temp in temperatures:
    filename = f"total_energy_T{temp:.2f}.txt"
   # Initialize an empty list to store energy values
    energies = []
    with open(filename, 'r') as file:
        for line in file:
         # Split the line into columns
           parts = line.split()
           energy = float(parts[0])
           energies.append(energy)

        # Compute the lowest energy if the list is not empty
    lowest_energies.append(min(energies))

# Filter out temperatures where no valid energy data was found
valid_temps = [t for t, e in zip(temperatures, lowest_energies) if e is not None]
valid_energies = [e for e in lowest_energies if e is not None]

# Plot the lowest energy as a function of temperature
plt.figure(figsize=(8, 6))
plt.plot(valid_temps, valid_energies, marker='o', label='Lowest Energy')
plt.xlabel('Temperature (T*)')
plt.ylabel('Lowest Energy')
plt.title('Lowest Energy vs. Temperature')
plt.grid(True)
plt.legend()
plt.savefig('lowest_energy_vs_temperature.png')
plt.show()

