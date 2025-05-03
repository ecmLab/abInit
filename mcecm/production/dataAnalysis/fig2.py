import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def generate_fig2_like_pattern(lattice_size, randomness=0.05):
    """
    Generate a spin configuration similar to Fig. 2 in the paper.

    Args:
        lattice_size (int): Size of the cubic lattice (NxNxN).
        randomness (float): Probability of introducing a random third spin.

    Returns:
        numpy.ndarray: 3D array of spins.
    """
    spins = np.zeros((lattice_size, lattice_size, lattice_size), dtype=int)

    for z in range(lattice_size):
        # Alternate spin pairs along the z-axis
        spin1, spin2 = (1, 2) if z % 2 == 0 else (2, 1)
        for x in range(lattice_size):
            for y in range(lattice_size):
                # Alternate spins on each face
                spins[x, y, z] = spin1 if (x + y) % 2 == 0 else spin2
                # Introduce randomness
                if np.random.rand() < randomness:
                    spins[x, y, z] = 3  # Add the third spin randomly

    return spins

def plot_3d_spin_configuration(spins, title):
    """
    Plot the 3D spin configuration using a voxel plot.

    Args:
        spins (numpy.ndarray): 3D array of spins.
        title (str): Title of the plot.
    """
    lattice_size = spins.shape[0]
    colors = ['red', 'green', 'blue']  # Map spin states to colors
    cmap = ListedColormap(colors)
    
    normalized_spins = spins - 1  # Normalize spins to range [0, 2] for colormap
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.voxels(np.ones_like(spins, dtype=bool), facecolors=cmap(normalized_spins), edgecolor='k', linewidth=0.2)
    
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

# Example Usage
lattice_size = 32  # Adjust lattice size as needed
spins = generate_fig2_like_pattern(lattice_size, randomness=0.05)
plot_3d_spin_configuration(spins, title="Fig. 2-Like Spin Pattern")
