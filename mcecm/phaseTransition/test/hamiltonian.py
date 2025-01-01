import numpy as np

# Define the cubic elastic constants (example values)
C11 = 200.0  # Elastic constant in GPa
C12 = 100.0  # Elastic constant in GPa
C44 = 50.0   # Elastic constant in GPa

# Define spontaneous strain tensors (example values for demonstration)
epsilon_c = 0.01  # Example value for epsilon_c
epsilon_a = -0.005  # Example value for epsilon_a
epsilon_1 = np.array([[epsilon_c, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_a]])
epsilon_2 = np.array([[epsilon_a, 0, 0], [0, epsilon_c, 0], [0, 0, epsilon_a]])
epsilon_3 = np.array([[epsilon_a, 0, 0], [0, epsilon_a, 0], [0, 0, epsilon_c]])

# Compute the stiffness tensor for cubic symmetry
def stiffness_tensor_cubic():
    """Construct the stiffness tensor for cubic symmetry."""
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

# Define the Green's function tensor
def green_function_tensor(n):
    """
    Compute the Green's function tensor Ω_ij(n) in reciprocal space.
    This is the inverse of C_ijkl n_k n_l.
    """
    n_k = np.array(n)  # Wavevector direction
    C = stiffness_tensor_cubic()  # Elastic modulus tensor (C_ijkl)

    # Compute Ω_ij^-1 = C_ijkl * n_k * n_l
    omega_inv = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    omega_inv[i, j] += C[i, j, k, l] * n[k] * n[l]

    # Invert Ω_ij^-1 to get Ω_ij
    omega = np.linalg.inv(omega_inv) if np.linalg.det(omega_inv) != 0 else np.zeros((3, 3))

    return omega

# Compute B_ijkl(n) using the Green's function tensor
def B_ijkl(n):
    """
    Compute B_ijkl(n) using the Green's function tensor.
    """
    C = stiffness_tensor_cubic()  # Full stiffness tensor
    omega = green_function_tensor(n)  # Green's function tensor Ω_ij(n)

    B = np.zeros((3, 3, 3, 3))

    # Compute B_ijkl
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    # First term: C_ijkl
                    B[i, j, k, l] = C[i, j, k, l]

                    # Second term: -n_p C_pqij Ω_qr(n) C_rskl n_s
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
    """Compute the Fourier transform of the strain field."""
    return np.fft.fftn(strain_field, axes=(0, 1, 2))

# Compute the elastic energy
def compute_elastic_energy(lattice_spins):
    """
    Compute the elastic energy based on the lattice spins and B_ijkl(n).
    Handles small lattices properly.
    """
    lattice_shape = lattice_spins.shape
    N = lattice_shape[0]  # Assume cubic lattice for simplicity
    strain_field = np.zeros((*lattice_shape, 3, 3))  # 3x3 strain tensor for each lattice site

    # Compute strain field in real space (from Equation 2)
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

    # Fourier transform the strain field
    strain_ft = fourier_transform_strain(strain_field)

    # Reciprocal-space grid
    q_range = np.fft.fftfreq(N, d=1 / (2 * np.pi))  # Wavevector grid in reciprocal space
    elastic_energy = 0

    # Loop over all wavevectors
    for qx in range(N):
        for qy in range(N):
            for qz in range(N):
                # Map indices to wavevector components
                k = np.array([q_range[qx], q_range[qy], q_range[qz]])
                n = k / np.linalg.norm(k) if np.linalg.norm(k) > 0 else np.zeros((3,))
                B = B_ijkl(n)

                # Access strain_ft using valid array indices
                epsilon_k = strain_ft[qx % N, qy % N, qz % N]  # Wrap indices if necessary
                epsilon_k_conj = np.conj(epsilon_k)
                elastic_energy += 0.5 * np.sum(B * epsilon_k * epsilon_k_conj)

    return np.real(elastic_energy)

# Example usage
#lattice_spins = np.random.choice([1, 2, 3], size=(4, 4, 4))  # Random lattice configuration
#elastic_energy = compute_elastic_energy(lattice_spins)
#print("Elastic Energy:", elastic_energy)
