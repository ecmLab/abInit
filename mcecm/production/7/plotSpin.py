import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def plot_all_lattices_in_one_figure(directory, lattice_size, output_filename=None):
    """
    Plot all 3D lattice configurations in one figure, with each configuration as a subplot.

    Args:
        directory (str): Directory containing the lattice configuration files.
        lattice_size (int): Size of the lattice (assumed cubic: lattice_size x lattice_size x lattice_size).
        output_filename (str, optional): If provided, saves the figure to this file.
    """
    # Find all files starting with "final_spins" or "maxStep_spins"
    files = sorted([f for f in os.listdir(directory) if f.startswith("spins_T0.398")])

    if not files:
        print("No files found starting with 'final_spins' or 'maxStep_spins'.")
        return

    # Define color map for spin states
    colors = ['red', 'green', 'blue']  # Red, Green, Blue for spin states 1, 2, 3
    cmap = ListedColormap(colors)

    # Calculate the grid size for subplots
    num_files = len(files)
    cols = 5  # Number of columns
    rows = (num_files + cols - 1) // cols  # Number of rows, rounded up

    # Create a figure
    fig = plt.figure(figsize=(4 * cols, 4 * rows))

    for idx, file in enumerate(files):
        # Load the spins from the file
        file_path = os.path.join(directory, file)
        spins = np.loadtxt(file_path, dtype=int).reshape((lattice_size, lattice_size, lattice_size))

        # Normalize spins to be between 0 and 2 (indices for the colormap)
        normalized_spins = spins - 1

        # Add a subplot
        ax = fig.add_subplot(rows, cols, idx + 1, projection='3d')

        # Plot voxels
        ax.voxels(
            np.ones_like(spins, dtype=bool),  # Show all voxels
            facecolors=cmap(normalized_spins),
            edgecolor='k'  # Add edges for clarity
        )

        # Set title
        ax.set_title(file, fontsize=10)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

    # Adjust layout
    plt.tight_layout()

    # Save the figure if output filename is provided
    if output_filename:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')

    # Show the figure
    plt.show()

# Example usage
lattice_size = 32  # Replace with the actual size of your lattice
directory = "./"  # Replace with the actual path to your directory
output_filename = "lattices_spin.png"  # Optional: save the figure
plot_all_lattices_in_one_figure(directory, lattice_size, output_filename)
