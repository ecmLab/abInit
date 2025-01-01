import matplotlib.pyplot as plt
import numpy as np

# Define Bragg peak positions and intensities for FCC Austenite
angles = [20, 29, 41, 50, 59, 69]  # Example 2θ values for (111), (200), (220), etc.
intensities = [100, 80, 60, 40, 30, 20]  # Relative intensities

# Plot the diffraction pattern
plt.figure(figsize=(8, 5))
plt.stem(angles, intensities, basefmt=" ")
plt.title("Bragg Reflection Peaks for Austenite (FCC)", fontsize=14)
plt.xlabel("2θ (Degrees)", fontsize=12)
plt.ylabel("Intensity (Arbitrary Units)", fontsize=12)
plt.xticks(np.arange(0, 81, 10))
plt.grid(alpha=0.5)
plt.tight_layout()
plt.show()

