import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def plot_3d_lattice_from_file(file_path, lattice_size, title, filename=None):
    """
    Plot the 3D lattice configuration from a file using voxels.

    Args:
        file_path (str): Path to the file containing the spin configuration.
        lattice_size (int): Size of the lattice (assumed cubic: lattice_size x lattice_size x lattice_size).
        title (str): Title of the plot.
        filename (str, optional): If provided, saves the plot to this file.
    """
    # Load the spins from the file
    spins = np.loadtxt(file_path, dtype=int).reshape((lattice_size, lattice_size, lattice_size))

    # Define color map for spin states
    colors = ['red', 'green', 'blue']  # Red, Green, Blue for spin states 1, 2, 3
    cmap = ListedColormap(colors)

    # Create a 3D plot
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Normalize spins to be between 0 and 2 (indices for the colormap)
    normalized_spins = spins - 1

    # Plot voxels
    ax.voxels(
        np.ones_like(spins, dtype=bool),  # Show all voxels
        facecolors=cmap(normalized_spins),
        edgecolor='k'  # Add edges for clarity
    )

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)

    if filename:
        plt.savefig(filename, dpi=300)
    plt.show()

# Example usage
file_path = "maxStep_spins.txt"  # Path to the file with the final spin configuration
lattice_size = 12  # Replace with the actual size of your lattice
plot_3d_lattice_from_file(file_path, lattice_size, title='3D Lattice Configuration', filename='lattice_3d.png')

