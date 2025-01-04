from mpi4py import MPI
import torch
import random

# Use GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Define the new dimensionless lattice misfit strains
def define_dimensionless_strains(epsilon_a, gamma_0):
    factor = epsilon_a / gamma_0
    epsilon_1 = torch.tensor([[1 + factor, 0, 0], [0, factor, 0], [0, 0, factor]], device=device)
    epsilon_2 = torch.tensor([[factor, 0, 0], [0, 1 + factor, 0], [0, 0, factor]], device=device)
    epsilon_3 = torch.tensor([[factor, 0, 0], [0, factor, 0], [0, 0, 1 + factor]], device=device)
    return epsilon_1, epsilon_2, epsilon_3

# Precompute stiffness tensor for cubic symmetry
def precompute_stiffness_tensor(anisoPar):
    C11 = 4 / anisoPar
    C12 = C11 / 2
    C44 = 1
    C = torch.zeros((3, 3, 3, 3), device=device)
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
    q_range = torch.fft.fftfreq(N, d=1 / (2 * torch.pi)).to(device)
    qx, qy, qz = torch.meshgrid(q_range, q_range, q_range, indexing="ij")
    q_grid = torch.stack((qx, qy, qz), dim=-1)
    B = torch.zeros((N, N, N, 3, 3, 3, 3), dtype=torch.complex64, device=device)
    for qx_idx in range(N):
        for qy_idx in range(N):
            for qz_idx in range(N):
                k = q_grid[qx_idx, qy_idx, qz_idx]
                if torch.linalg.norm(k) == 0:
                    continue
                k_unit = k / torch.linalg.norm(k)
                G_inv = torch.einsum("ijmn,m,n->ij", stiffness_tensor, k_unit, k_unit)
                G = torch.linalg.inv(G_inv) if torch.linalg.det(G_inv) != 0 else torch.zeros((3, 3), dtype=torch.complex64, device=device)
                B[qx_idx, qy_idx, qz_idx] = stiffness_tensor - torch.einsum(
                    "ab,acbd,efcd->abef",
                    G,
                    stiffness_tensor,
                    stiffness_tensor,
                )
    return q_grid, B


# Ensure strain_ft and related tensors are created as complex types
def compute_strain_field(lattice_spins, epsilon_1, epsilon_2, epsilon_3):
    N = lattice_spins.shape[0]
    strain_field = torch.zeros((N, N, N, 3, 3), dtype=torch.complex64, device=device)
    for x in range(N):
        for y in range(N):
            for z in range(N):
                spin = lattice_spins[x, y, z].item()
                if spin == 1:
                    strain_field[x, y, z] = epsilon_1.to(dtype=torch.complex64)
                elif spin == 2:
                    strain_field[x, y, z] = epsilon_2.to(dtype=torch.complex64)
                elif spin == 3:
                    strain_field[x, y, z] = epsilon_3.to(dtype=torch.complex64)
    return strain_field

# Compute elastic energy
def compute_elastic_energy(lattice_spins, q_grid, B, strain_ft, rank, size):
    N = lattice_spins.shape[0]
    local_energy = 0
    for qx_idx in range(rank, N, size):
        for qy_idx in range(N):
            for qz_idx in range(N):
                if torch.allclose(q_grid[qx_idx, qy_idx, qz_idx], torch.tensor(0.0, device=device)):
                    continue
                local_energy += torch.real(
                    torch.einsum(
                        "ij,ijkl,kl->",
                        strain_ft[qx_idx, qy_idx, qz_idx],
                        B[qx_idx, qy_idx, qz_idx],
                        strain_ft[qx_idx, qy_idx, qz_idx].conj(),
                    )
                )
    return local_energy

# Monte Carlo step with MPI
def monte_carlo_step(lattice_spins, temperature, q_grid, B, strain_ft, epsilon_1, epsilon_2, epsilon_3, rank, size):
    N = lattice_spins.shape[0]
    x, y, z = torch.randint(0, N, (3,), device=device)
    current_spin = lattice_spins[x, y, z].item()
    proposed_spin = random.choice([s for s in [1, 2, 3] if s != current_spin])

    # Compute new energy
    lattice_spins[x, y, z] = proposed_spin
    new_strain_field = compute_strain_field(lattice_spins, epsilon_1, epsilon_2, epsilon_3)
    strain_ft_new = torch.fft.fftn(new_strain_field, dim=(0, 1, 2))
    E_new = compute_elastic_energy(lattice_spins, q_grid, B, strain_ft_new, rank, size)

    # Revert to original state
    lattice_spins[x, y, z] = current_spin
    original_strain_field = compute_strain_field(lattice_spins, epsilon_1, epsilon_2, epsilon_3)
    strain_ft_original = torch.fft.fftn(original_strain_field, dim=(0, 1, 2))
    E_old = compute_elastic_energy(lattice_spins, q_grid, B, strain_ft_original, rank, size)

    delta_E = E_new - E_old
    max_exp_arg = 700
    prob = torch.exp(-torch.clamp(delta_E / temperature, max=max_exp_arg)).item() if delta_E > 0 else 1.0

    if random.random() < prob:
        lattice_spins[x, y, z] = proposed_spin
        strain_ft[:] = strain_ft_new
        return True
    return False

# Main simulation
def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    lattice_size = 10
    lattice_spins = torch.randint(1, 4, (lattice_size, lattice_size, lattice_size), device=device)
    epsilon_a = 0.1
    gamma_0 = 0.4
    anisoPar = 1
    temperature = 0.20
    num_steps = 4000

    epsilon_1, epsilon_2, epsilon_3 = define_dimensionless_strains(epsilon_a, gamma_0)
    stiffness_tensor = precompute_stiffness_tensor(anisoPar)
    q_grid, B = precompute_reciprocal_space_and_kernel(lattice_size, stiffness_tensor)
    strain_field = compute_strain_field(lattice_spins, epsilon_1, epsilon_2, epsilon_3)
    strain_ft = torch.fft.fftn(strain_field, dim=(0, 1, 2))

    zero_move_steps = 0
    terminate_flag = False

    for step in range(1, num_steps + 1):
        if terminate_flag:
            break

        accepted_moves = 0
        for _ in range(lattice_spins.numel() // size):
            if monte_carlo_step(lattice_spins, temperature, q_grid, B, strain_ft, epsilon_1, epsilon_2, epsilon_3, rank, size):
                accepted_moves += 1

        total_accepted_moves = comm.reduce(accepted_moves, op=MPI.SUM, root=0)
        if rank == 0:
            if total_accepted_moves <= 1:
                zero_move_steps += 1
            else:
                zero_move_steps = 0
            if zero_move_steps >= 5:
                torch.save(lattice_spins.cpu(), "final_spins.pt")
                terminate_flag = True

        terminate_flag = comm.bcast(terminate_flag, root=0)

if __name__ == "__main__":
    main()

