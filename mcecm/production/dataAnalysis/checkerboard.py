import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Parameters for lattice size and spin states
lattice_size = 32  # 32x32x32 lattice
spin_states = [1, 2, 3]  # Possible spin states for checkerboard

# Initialize the lattice with checkerboard pattern
lattice_spins = np.zeros((lattice_size, lattice_size, lattice_size), dtype=int)

for x in range(lattice_size):
    for y in range(lattice_size):
        for z in range(lattice_size):
            lattice_spins[x, y, z] = spin_states[(x + y + z) % len(spin_states)]

# Define color map for spin states
colors = ['red', 'green', 'blue']  # Assign colors to spin states 1, 2, 3
cmap = ListedColormap(colors)

# Create a 3D plot
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Normalize spins to be between 0 and 2 (indices for the colormap)
normalized_spins = lattice_spins - 1

# Plot voxels
ax.voxels(
    np.ones_like(lattice_spins, dtype=bool),  # Show all voxels
    facecolors=cmap(normalized_spins),
    edgecolor='k'  # Add edges for clarity
)

# Set labels and title
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_title('3D Checkerboard Spin Configuration')

plt.show()
