import numpy as np
import random

# Define elastic constants
C11 = 200.0  # Elastic constant in GPa
C12 = 100.0  # Elastic constant in GPa
C44 = 50.0   # Elastic constant in GPa

# Define spontaneous strain tensors
epsilon_c = 0.01
epsilon_a = -0.005
epsilon_1 = np.array([[epsilon_c, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_a]])
epsilon_2 = np.array([[epsilon_a, 0, 0], [0, epsilon_c, 0], [0, 0, epsilon_a]])
epsilon_3 = np.array([[epsilon_a, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_c]])

# Function to read the initial lattice configuration from a LAMMPS dump file
def read_lattice_config(filename):
    """
    Read the lattice spin configuration from a LAMMPS dump file.
    Handles non-cubic lattices.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Extract the number of atoms
    num_atoms = int(lines[3].strip())
    
    # Extract box bounds to determine lattice dimensions
    x_bounds = list(map(float, lines[5].strip().split()))
    y_bounds = list(map(float, lines[6].strip().split()))
    z_bounds = list(map(float, lines[7].strip().split()))

    # Calculate lattice dimensions
    x_len = int(round((x_bounds[1] - x_bounds[0]) / 2.0))  # Lattice constant assumed 2.0
    y_len = int(round((y_bounds[1] - y_bounds[0]) / 2.0))
    z_len = int(round((z_bounds[1] - z_bounds[0]) / 2.0))

    # Ensure the total number of lattice sites matches the number of atoms
    if x_len * y_len * z_len != num_atoms:
        raise ValueError(f"Lattice dimensions ({x_len}, {y_len}, {z_len}) do not match number of atoms ({num_atoms}).")

    # Extract atom data (assumes x, y, z, and type are included in the dump)
    atom_data = []
    for line in lines[9:9 + num_atoms]:
        parts = line.strip().split()
        atom_type = int(parts[1])  # Atom type
        atom_data.append(atom_type)

    # Reshape the lattice spins into a 3D array
    lattice_spins = np.array(atom_data).reshape((x_len, y_len, z_len))

    return lattice_spins


# Stiffness tensor for cubic symmetry
def stiffness_tensor_cubic():
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if i == j and k == l:
                        C[i, j, k, l] = C11 if i == k else C12
                    elif i == k and j == l:
                        C[i, j, k, l] = C44 if i != j else 0
    return C

# Green's function tensor
def green_function_tensor(n):
    n_k = np.array(n)
    C = stiffness_tensor_cubic()
    omega_inv = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    omega_inv[i, j] += C[i, j, k, l] * n_k[k] * n_k[l]
    omega = np.linalg.inv(omega_inv) if np.linalg.det(omega_inv) != 0 else np.zeros((3, 3))
    return omega

# B_ijkl calculation
def B_ijkl(n):
    C = stiffness_tensor_cubic()
    omega = green_function_tensor(n)
    B = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    B[i, j, k, l] = C[i, j, k, l]
                    for p in range(3):
                        for q in range(3):
                            for r in range(3):
                                for s in range(3):
                                    B[i, j, k, l] -= (
                                        n[p]
                                        * C[p, q, i, j]
                                        * omega[q, r]
                                        * C[r, s, k, l]
                                        * n[s]
                                    )
    return B

# Fourier transform of the strain field
def fourier_transform_strain(strain_field):
    return np.fft.fftn(strain_field, axes=(0, 1, 2))

# Elastic energy calculation
def compute_elastic_energy(lattice_spins):
    N = lattice_spins.shape[0]
    strain_field = np.zeros((*lattice_spins.shape, 3, 3))
    for x in range(N):
        for y in range(N):
            for z in range(N):
                spin = lattice_spins[x, y, z]
                if spin == 1:
                    strain_field[x, y, z] = epsilon_1
                elif spin == 2:
                    strain_field[x, y, z] = epsilon_2
                elif spin == 3:
                    strain_field[x, y, z] = epsilon_3

    strain_ft = fourier_transform_strain(strain_field)
    q_range = np.fft.fftfreq(N, d=1 / (2 * np.pi))
    elastic_energy = 0
    for qx in range(N):
        for qy in range(N):
            for qz in range(N):
                k = np.array([q_range[qx], q_range[qy], q_range[qz]])
                n = k / np.linalg.norm(k) if np.linalg.norm(k) > 0 else np.zeros((3,))
                B = B_ijkl(n)
                epsilon_k = strain_ft[qx, qy, qz]
                epsilon_k_conj = np.conj(epsilon_k)
                elastic_energy += 0.5 * np.sum(B * epsilon_k * epsilon_k_conj)
    return np.real(elastic_energy)

# Monte Carlo step for spin updates
def monte_carlo_step(lattice_spins, temperature):
    N = lattice_spins.shape[0]
    x, y, z = np.random.randint(0, N, size=3)
    current_spin = lattice_spins[x, y, z]
    proposed_spin = random.choice([s for s in [1, 2, 3] if s != current_spin])
    lattice_spins[x, y, z] = proposed_spin
    E_new = compute_elastic_energy(lattice_spins)
    lattice_spins[x, y, z] = current_spin
    E_old = compute_elastic_energy(lattice_spins)
    delta_E = E_new - E_old
    if delta_E <= 0 or random.random() < np.exp(-delta_E / temperature):
        lattice_spins[x, y, z] = proposed_spin
        return True
    return False

def main():
    lattice_spins = read_lattice_config("init_config.dump")  # Load lattice from LAMMPS dump
    temperature = 300.0
    num_steps = 1000
    for step in range(num_steps):
        accepted_moves = 0
        for _ in range(lattice_spins.size):
            if monte_carlo_step(lattice_spins, temperature):
                accepted_moves += 1
        print(f"Step {step + 1}: Accepted moves = {accepted_moves}")
    np.savetxt("final_spins.txt", lattice_spins.reshape(-1), fmt='%d')

if __name__ == "__main__":
    main()
