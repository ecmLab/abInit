import numpy as np
import random
from joblib import Parallel, delayed

# Elastic constants (arbitrary units)
C11 = 200.0
C12 = 100.0
C44 = 50.0

# Strain tensors for spin states
epsilon_c = 0.01
epsilon_a = -0.005
epsilon_1 = np.array([[epsilon_c, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_a]])
epsilon_2 = np.array([[epsilon_a, 0, 0], [0, epsilon_c, 0], [0, 0, epsilon_a]])
epsilon_3 = np.array([[epsilon_a, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_c]])

# Precompute stiffness tensor for cubic symmetry
def precompute_stiffness_tensor():
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

# Compute reciprocal space contribution for a single wavevector
def reciprocal_energy_contribution(k, strain_ft, stiffness_tensor, idx):
    if np.linalg.norm(k) == 0:
        return 0  # Skip zero wavevector

    # Normalize the wavevector
    n = k / np.linalg.norm(k)

    # Compute Green's function tensor
    G_inv = np.einsum("ijmn,m,n->ij", stiffness_tensor, k, k)  # Efficient contraction
    G = np.linalg.inv(G_inv) if np.linalg.det(G_inv) != 0 else np.zeros((3, 3))

    # Access the Fourier-transformed strain field
    epsilon_k = strain_ft[idx[0], idx[1], idx[2]]
    epsilon_k_conj = np.conj(epsilon_k)

    # Compute elastic energy contribution
    energy = np.real(np.sum(epsilon_k * np.einsum("ij,jk,kl->il", G, epsilon_k, epsilon_k_conj)))
    return energy

# Compute the total elastic energy in parallel
def compute_elastic_energy(lattice_spins, stiffness_tensor, reciprocal_grid, n_jobs=-1):
    """
    Compute the total elastic energy using Fourier transforms and parallelization.
    """
    N = lattice_spins.shape[0]
    strain_field = np.zeros((*lattice_spins.shape, 3, 3))

    # Assign strain tensors
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

    # Fourier transform of strain field
    strain_ft = np.fft.fftn(strain_field, axes=(0, 1, 2))

    # Prepare arguments for parallelization
    tasks = [
        (reciprocal_grid[qx, qy, qz], strain_ft, stiffness_tensor, (qx, qy, qz))
        for qx in range(N) for qy in range(N) for qz in range(N)
    ]

    # Parallel computation of elastic energy contributions
    results = Parallel(n_jobs=n_jobs)(delayed(reciprocal_energy_contribution)(*task) for task in tasks)

    # Sum up contributions
    elastic_energy = sum(results)
    return elastic_energy

# Monte Carlo step
def monte_carlo_step(lattice_spins, temperature, stiffness_tensor, reciprocal_grid, n_jobs=-1):
    """
    Perform a single Monte Carlo step for spin updates.
    """
    N = lattice_spins.shape[0]
    x, y, z = np.random.randint(0, N, size=3)  # Randomly select a lattice site
    current_spin = lattice_spins[x, y, z]
    proposed_spin = random.choice([s for s in [1, 2, 3] if s != current_spin])
    lattice_spins[x, y, z] = proposed_spin  # Temporarily apply the proposed spin
    E_new = compute_elastic_energy(lattice_spins, stiffness_tensor, reciprocal_grid, n_jobs=n_jobs)
    lattice_spins[x, y, z] = current_spin  # Revert to the original spin
    E_old = compute_elastic_energy(lattice_spins, stiffness_tensor, reciprocal_grid, n_jobs=n_jobs)
    delta_E = E_new - E_old

    # Metropolis acceptance criterion
    if delta_E <= 0 or random.random() < np.exp(-delta_E / temperature):
        lattice_spins[x, y, z] = proposed_spin  # Accept the proposed spin
        return True
    return False

# Precompute reciprocal space grid
def precompute_reciprocal_space(N):
    """
    Precompute reciprocal space wavevector grid for a cubic lattice.
    """
    q_range = np.fft.fftfreq(N, d=1 / (2 * np.pi))
    qx, qy, qz = np.meshgrid(q_range, q_range, q_range, indexing="ij")
    return np.stack((qx, qy, qz), axis=-1)

# Main simulation
def main():
    lattice_size = 10  # Lattice dimensions
    lattice_spins = np.random.choice([1, 2, 3], size=(lattice_size, lattice_size, lattice_size))
    temperature = 300.0  # Simulation temperature
    num_steps = 1000  # Total number of Monte Carlo steps

    # Precompute stiffness tensor and reciprocal space wavevector grid
    stiffness_tensor = precompute_stiffness_tensor()
    reciprocal_grid = precompute_reciprocal_space(lattice_size)

    for step in range(1, num_steps + 1):
        accepted_moves = 0
        for _ in range(lattice_spins.size):  # Attempt one move per lattice site
            if monte_carlo_step(lattice_spins, temperature, stiffness_tensor, reciprocal_grid, n_jobs=-1):
                accepted_moves += 1

        # Print progress every 500 steps
        print(f"Timestep {step}: Accepted moves = {accepted_moves}")

    # Save final lattice configuration
    np.savetxt("final_spins.txt", lattice_spins.reshape(-1), fmt='%d')
    print("Simulation complete. Final spins saved to 'final_spins.txt'.")

if __name__ == "__main__":
    main()

