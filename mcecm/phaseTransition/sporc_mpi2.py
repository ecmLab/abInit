from mpi4py import MPI
import numpy as np
import random

# Define the new dimensionless lattice misfit strains
def define_dimensionless_strains(epsilon_a, gamma_0):
    factor = epsilon_a / gamma_0
    epsilon_1 = np.array([[1 + factor, 0, 0], [0, factor, 0], [0, 0, factor]])
    epsilon_2 = np.array([[factor, 0, 0], [0, 1 + factor, 0], [0, 0, factor]])
    epsilon_3 = np.array([[factor, 0, 0], [0, factor, 0], [0, 0, 1 + factor]])
    return epsilon_1, epsilon_2, epsilon_3

# Precompute stiffness tensor for cubic symmetry
def precompute_stiffness_tensor(anisoPar):
    C11 = 4 / anisoPar
    C12 = C11 / 2
    C44 = 1
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

# Precompute reciprocal space grid and kernel
def precompute_reciprocal_space_and_kernel(N, stiffness_tensor):
    q_range = np.fft.fftfreq(N, d=1 / (2 * np.pi))
    qx, qy, qz = np.meshgrid(q_range, q_range, q_range, indexing="ij")
    q_grid = np.stack((qx, qy, qz), axis=-1)
    B = np.zeros((N, N, N, 3, 3, 3, 3))
    for qx_idx in range(N):
        for qy_idx in range(N):
            for qz_idx in range(N):
                k = q_grid[qx_idx, qy_idx, qz_idx]
                if np.linalg.norm(k) == 0:
                    continue
                k_unit = k / np.linalg.norm(k)
                G_inv = np.einsum("ijmn,m,n->ij", stiffness_tensor, k_unit, k_unit)
                G = np.linalg.inv(G_inv) if np.linalg.det(G_inv) != 0 else np.zeros((3, 3))
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            for l in range(3):
                                B[qx_idx, qy_idx, qz_idx, i, j, k, l] = (
                                    stiffness_tensor[i, j, k, l]
                                    - np.einsum(
                                        "mn,im,nj->",
                                        G,
                                        stiffness_tensor[:, :, i, j],
                                        stiffness_tensor[:, :, k, l],
                                    )
                                )
    return q_grid, B

# Compute strain field for a single lattice site update
def update_strain_field(strain_field, x, y, z, spin, epsilon_tensors):
    strain_field[x, y, z] = epsilon_tensors[spin - 1]

# Compute elastic energy
def compute_elastic_energy(q_grid, B, strain_ft, rank, size):
    N = strain_ft.shape[0]
    local_energy = 0
    for qx_idx in range(rank, N, size):  # Divide work among ranks
        for qy_idx in range(N):
            for qz_idx in range(N):
                if np.allclose(q_grid[qx_idx, qy_idx, qz_idx], 0):
                    continue
                local_energy += np.real(
                    np.einsum(
                        "ij,ijkl,kl->",
                        strain_ft[qx_idx, qy_idx, qz_idx],
                        B[qx_idx, qy_idx, qz_idx],
                        strain_ft[qx_idx, qy_idx, qz_idx].conj(),
                    )
                )
    return local_energy

# Monte Carlo step with MPI
def monte_carlo_step(lattice_spins, temperature, q_grid, B, strain_field, epsilon_tensors, rank, size):
    N = lattice_spins.shape[0]
    x, y, z = np.random.randint(0, N, size=3)
    current_spin = lattice_spins[x, y, z]
    proposed_spin = random.choice([s for s in [1, 2, 3] if s != current_spin])

    # Compute new energy
    update_strain_field(strain_field, x, y, z, proposed_spin, epsilon_tensors)
    strain_ft_new = np.fft.fftn(strain_field, axes=(0, 1, 2))
    E_new = compute_elastic_energy(q_grid, B, strain_ft_new, rank, size)

    # Revert to original state
    update_strain_field(strain_field, x, y, z, current_spin, epsilon_tensors)
    strain_ft_original = np.fft.fftn(strain_field, axes=(0, 1, 2))
    E_old = compute_elastic_energy(q_grid, B, strain_ft_original, rank, size)

    # Calculate energy difference
    delta_E = E_new - E_old
    max_exp_arg = 700
    # Metropolis acceptance criterion
    if delta_E > 0:
       prob = np.exp(-np.clip(delta_E / temperature, None, max_exp_arg))
    else:
       prob = 1.0

    if random.random() < prob:
        lattice_spins[x, y, z] = proposed_spin
        update_strain_field(strain_field, x, y, z, proposed_spin, epsilon_tensors)
        return 1
    return 0

# Main simulation
def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    lattice_size = 32
    lattice_spins = np.random.choice([1, 2, 3], size=(lattice_size, lattice_size, lattice_size))
    epsilon_a = 0.1
    gamma_0 = 0.4
    anisoPar = 1
    temperature = 0.20
    num_steps = 4000

    epsilon_1, epsilon_2, epsilon_3 = define_dimensionless_strains(epsilon_a, gamma_0)
    stiffness_tensor = precompute_stiffness_tensor(anisoPar)
    q_grid, B = precompute_reciprocal_space_and_kernel(lattice_size, stiffness_tensor)
    strain_field = np.zeros((lattice_size, lattice_size, lattice_size, 3, 3))
    epsilon_tensors = [epsilon_1, epsilon_2, epsilon_3]

    zero_move_steps = 0  # Track successive steps with zero accepted moves
    terminate_flag = False  # Termination flag shared across all ranks
    for step in range(1, num_steps + 1):
        accepted_moves = 0

        if terminate_flag:  # Check if termination flag is set
           break

        for _ in range(lattice_spins.size // lattice_size):
            accepted_moves += monte_carlo_step(
                lattice_spins, temperature, q_grid, B, strain_field, epsilon_tensors, rank, size
            )

        total_accepted_moves = comm.reduce(accepted_moves, op=MPI.SUM, root=0)
        total_accepted_moves = int(total_accepted_moves or 0)
        if rank == 0:
           print(f"Timestep {step}: Accepted moves = {total_accepted_moves}")
           
           if total_accepted_moves <= 1:
                zero_move_steps += 1
           else:
                zero_move_steps = 0
           if zero_move_steps >= 5:
                print("Termination criteria met: Zero accepted moves for five successive steps.")
       # Save final lattice configuration
                np.savetxt("final_spins.txt", lattice_spins.reshape(-1), fmt='%d')
                print("Simulation complete. Final spins saved to 'final_spins.txt'.")
                terminate_flag = True
       
       # Broadcast termination flag to all ranks
        terminate_flag = comm.bcast(terminate_flag, root=0)

if __name__ == "__main__":
    main()

