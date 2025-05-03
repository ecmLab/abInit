import os
import numpy as np
import matplotlib.pyplot as plt

# Directory containing the strain files
strain_files_directory = "../results"  # Update this to your directory

# Initialize data storage for strain11
Nsample = 3000
strain11_data = {}
strain22_data = {}
strain33_data = {}
temperature_list = []
strain_deviation = []

# Load all strain files
for filename in sorted(os.listdir(strain_files_directory)):
    if filename.startswith("macro_strain_T") and filename.endswith(".txt"):
        # Extract temperature from the filename
        temp = filename.split("_T")[1].replace(".txt", "")
        temp = float(temp)  # Convert to float for sorting

        # Read the strain tensor values
        file_path = os.path.join(strain_files_directory, filename)
        strains = []
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for i in range(0, len(lines), 3):  # Read three lines at a time
                # Parse three lines to form a 3x3 matrix
                row1 = list(map(float, lines[i].strip().replace("[", "").replace("]", "").split()))
                row2 = list(map(float, lines[i + 1].strip().replace("[", "").replace("]", "").split()))
                row3 = list(map(float, lines[i + 2].strip().replace("[", "").replace("]", "").split()))
                matrix = np.array([row1, row2, row3])
                strains.append(matrix)

        # Extract strain11 component for each Monte Carlo step
        strain11_values = [strain[0, 0] for strain in strains]
        strain22_values = [strain[1, 1] for strain in strains]
        strain33_values = [strain[2, 2] for strain in strains]

        stdev11 = np.std(strain11_values[-Nsample:])
        stdev22 = np.std(strain22_values[-Nsample:])
        stdev33 = np.std(strain33_values[-Nsample:])

        # Store the strain11 data with their corresponding temperature
        strain11_data[temp] = strain11_values
        strain22_data[temp] = strain22_values
        strain33_data[temp] = strain33_values
        delta  = 1./(stdev11**2 + stdev22**2 + stdev33**2)
        strain_deviation.append(delta)
        temperature_list.append(temp)

# Sort the data by temperature
#sorted_temps = sorted(strain11_data.keys())
#print(strain_deviation)

# Plot the strain11 data
plt.figure(figsize=(8, 6))
nt = -10
plt.plot(temperature_list[nt:], strain_deviation[nt:])
#for temp in sorted_temps:
#    mc_steps = range(1, len(strain11_data[temp]) + 1)
#    plt.plot(mc_steps, strain11_data[temp], label=f"T={temp:.2f}")

# Customize the plot
plt.xlabel("MC Steps")
plt.ylabel("Strain Delta")
plt.legend(fontsize=7, loc="upper right")
plt.grid(True)

# Save and show the plot
#plt.savefig("strain11_vs_mc_steps_3.png")
plt.show()
