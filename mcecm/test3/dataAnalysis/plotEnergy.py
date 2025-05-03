import os
import matplotlib.pyplot as plt

# Directory containing the energy files
energy_files_directory = "../5"  # Update this to your directory

# Initialize data storage
energy_data = {}

# Load all energy files
for filename in sorted(os.listdir(energy_files_directory)):
    if filename.startswith("total_energy_T") and filename.endswith(".txt"):
        # Extract temperature from the filename
        temp = filename.split("_T")[1].replace(".txt", "")
#        print(temp)
        temp = float(temp)  # Convert to float for sorting later

        # Read the energy values
        file_path = os.path.join(energy_files_directory, filename)
        with open(file_path, 'r') as f:
            energies = [float(line.strip()) for line in f.readlines()]

        # Store the energies with their corresponding temperature
        energy_data[temp] = energies

# Sort the data by temperature
sorted_temps = sorted(energy_data.keys())

# Plot the data
plt.figure(figsize=(5, 4))
for temp in sorted_temps:
    mc_steps = range(1, len(energy_data[temp]) + 1)
    plt.plot(mc_steps, energy_data[temp], label=f"T={temp:.2f}")

# Customize the plot
plt.xlabel("MC Steps")
plt.ylabel("Total Elastic Energy")
#plt.title("Total Energy as a Function of MC Steps")
plt.legend(fontsize=7)
plt.grid(True)

# Save and show the plot
plt.savefig("energy_5.png")
plt.show()
