import numpy as np
import random

# Define elastic constants (in arbitrary units)
C11 = 200.0  # Elastic constant
C12 = 100.0  # Elastic constant
C44 = 50.0   # Elastic constant

# Define spontaneous strain tensors for the three spin states
epsilon_c = 0.01
epsilon_a = -0.005
epsilon_1 = np.array([[epsilon_c, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_a]])
epsilon_2 = np.array([[epsilon_a, 0, 0], [0, epsilon_c, 0], [0, 0, epsilon_a]])
epsilon_3 = np.array([[epsilon_a, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_c]])

# Compute the stiffness tensor for cubic symmetry
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

# Compute the elastic energy
def compute_elastic_energy(lattice_spins):
    """
    Compute the total elastic energy of the system based on strain tensors.
    """
    N = lattice_spins.shape[0]  # Lattice size
    strain_field = np.zeros((*lattice_spins.shape, 3, 3))  # 3x3 strain tensor at each site

    # Assign strain tensors based on spin state
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

    # Fourier transform of the strain field
    strain_ft = np.fft.fftn(strain_field, axes=(0, 1, 2))

    # Reciprocal space grid
    q_range = np.fft.fftfreq(N, d=1 / (2 * np.pi))  # Wavevectors in reciprocal space
    elastic_energy = 0
    stiffness_tensor = stiffness_tensor_cubic()

    # Loop over all reciprocal space points
    for qx in range(N):
        for qy in range(N):
            for qz in range(N):
                k = np.array([q_range[qx], q_range[qy], q_range[qz]])
                if np.linalg.norm(k) == 0:
                    continue  # Skip the zero wavevector

                # Normalize the wavevector
                n = k / np.linalg.norm(k)

                # Compute Green's function tensor
                G_inv = np.zeros((3, 3))
                for i in range(3):
                    for j in range(3):
                        for m in range(3):
                            for n in range(3):
                                G_inv[i, j] += stiffness_tensor[i, j, m, n] * k[m] * k[n]
                G = np.linalg.inv(G_inv) if np.linalg.det(G_inv) != 0 else np.zeros((3, 3))

                # Elastic energy calculation in reciprocal space
                epsilon_k = strain_ft[qx, qy, qz]
                epsilon_k_conj = np.conj(epsilon_k)
                for i in range(3):
                    for j in range(3):
                        elastic_energy += np.real(np.dot(epsilon_k[i, j], epsilon_k_conj[i, j]))

    return elastic_energy

# Monte Carlo step
def monte_carlo_step(lattice_spins, temperature):
    """
    Perform a single Monte Carlo step for spin updates.
    """
    N = lattice_spins.shape[0]
    x, y, z = np.random.randint(0, N, size=3)  # Randomly select a lattice site
    current_spin = lattice_spins[x, y, z]
    proposed_spin = random.choice([s for s in [1, 2, 3] if s != current_spin])
    lattice_spins[x, y, z] = proposed_spin  # Temporarily apply the proposed spin
    E_new = compute_elastic_energy(lattice_spins)
    lattice_spins[x, y, z] = current_spin  # Revert to the original spin
    E_old = compute_elastic_energy(lattice_spins)
    delta_E = E_new - E_old

    # Metropolis acceptance criterion
    if delta_E <= 0 or random.random() < np.exp(-delta_E / temperature):
        lattice_spins[x, y, z] = proposed_spin  # Accept the proposed spin
        return True
    return False

# Main simulation
def main():
    lattice_size = 10  # Lattice dimensions
    lattice_spins = np.random.choice([1, 2, 3], size=(lattice_size, lattice_size, lattice_size))
    temperature = 300.0  # Simulation temperature
    num_steps = 10000  # Total number of Monte Carlo steps

    for step in range(1, num_steps + 1):
        accepted_moves = 0
        for _ in range(lattice_spins.size):  # Attempt one move per lattice site
            if monte_carlo_step(lattice_spins, temperature):
                accepted_moves += 1

        # Print progress every 500 steps
        print(f"Timestep {step}: Accepted moves = {accepted_moves}")

    # Save final lattice configuration
    np.savetxt("final_spins.txt", lattice_spins.reshape(-1), fmt='%d')
    print("Simulation complete. Final spins saved to 'final_spins.txt'.")

if __name__ == "__main__":
    main()

