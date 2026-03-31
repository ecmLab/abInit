"""
Set up NEB migration-barrier calculations for the three most promising phases:
  mp-560538_Na3AlS3    → 1×1×1 unit cell (24 Na)
  user_Na5AlS4_Na5AlS4 → 1×1×1 unit cell (40 Na)
  JVASP-12818_Na5InS4  → 1×2×1 supercell  (20 Na)

Strategy:
  1. Find the shortest Na-Na distance in the structure (primary hop).
  2. Remove the destination Na (vacancy), creating n-1 atom endpoint structures.
  3. Image 00: hopping Na at its equilibrium site A.
  4. Image 06: hopping Na displaced to site B (= removed Na position).
  5. Images 01-05: linear interpolation with PBC-aware shortest path.

Output tree:
  <phase>/neb/hop1/{00..06}/POSCAR   + INCAR, KPOINTS, POTCAR, submit.sh

Run from calculation/:  python setup_neb.py
"""

import shutil
import numpy as np
from pathlib import Path

CALC_DIR = Path(__file__).parent
N_IMAGES = 5          # interpolated images between endpoints (not counting 00 and 06)

PHASES = {
    "mp-560538_Na3AlS3":    {"supercell": (1, 1, 1)},
    "user_Na5AlS4_Na5AlS4": {"supercell": (1, 1, 1)},
    "JVASP-12818_Na5InS4":  {"supercell": (1, 2, 1)},
}


# ── POSCAR helpers ────────────────────────────────────────────────────────
def read_poscar(path):
    lines = Path(path).read_text().splitlines()
    scale   = float(lines[1])
    lattice = np.array([l.split() for l in lines[2:5]], dtype=float) * scale
    elements = lines[5].split()
    counts   = list(map(int, lines[6].split()))
    ctype    = lines[7].strip()
    n = sum(counts)
    pos = np.array([l.split()[:3] for l in lines[8:8+n]], dtype=float)
    if ctype.lower().startswith("c"):
        pos = np.linalg.solve(lattice.T, pos.T).T   # cart → frac
    return dict(lattice=lattice, elements=elements, counts=counts, positions=pos)


def make_supercell(p, sc):
    na, nb, nc = sc
    new_lat = p["lattice"].copy()
    new_lat[0] *= na; new_lat[1] *= nb; new_lat[2] *= nc
    new_pos, new_counts = [], []
    offset = 0
    for cnt in p["counts"]:
        atom_pos = p["positions"][offset:offset+cnt]
        expanded = []
        for ia in range(na):
            for ib in range(nb):
                for ic in range(nc):
                    shifted = (atom_pos + [ia, ib, ic]) / [na, nb, nc]
                    expanded.append(shifted)
        new_pos.extend(np.vstack(expanded))
        new_counts.append(cnt * na * nb * nc)
        offset += cnt
    return dict(lattice=new_lat, elements=p["elements"],
                counts=new_counts, positions=np.array(new_pos))


def write_poscar(p, path, comment="NEB image"):
    lines = [comment, "  1.0"]
    for row in p["lattice"]:
        lines.append("  " + "  ".join(f"{v:22.16f}" for v in row))
    lines.append("  " + "  ".join(p["elements"]))
    lines.append("  " + "  ".join(map(str, p["counts"])))
    lines.append("Direct")
    for pos in p["positions"]:
        lines.append("  " + "  ".join(f"{v:.16f}" for v in pos))
    Path(path).write_text("\n".join(lines) + "\n")


# ── Hop finder ────────────────────────────────────────────────────────────
def find_min_na_hop(p):
    """Return (idx_i, idx_j, distance_Å) for the shortest Na-Na pair."""
    # Build flat element list
    elem_list = []
    for el, cnt in zip(p["elements"], p["counts"]):
        elem_list.extend([el] * cnt)

    na_idx = [i for i, e in enumerate(elem_list) if e == "Na"]
    pos    = p["positions"]
    lat    = p["lattice"]

    best = (None, None, np.inf)
    for a in range(len(na_idx)):
        for b in range(a+1, len(na_idx)):
            i, j = na_idx[a], na_idx[b]
            diff = pos[j] - pos[i]
            diff -= np.round(diff)          # minimum image convention
            dist = np.linalg.norm(diff @ lat)
            if dist < best[2]:
                best = (i, j, dist)
    return best


# ── Image builder ─────────────────────────────────────────────────────────
def build_neb_images(p, idx_hop, idx_vac, n_images=5):
    """
    Create n_images+2 structures (00 … N+1).
    - idx_hop : index of the Na atom that moves (stays in structure)
    - idx_vac : index of the Na atom that is removed (becomes the vacancy / destination)
    """
    elem_list = []
    for el, cnt in zip(p["elements"], p["counts"]):
        elem_list.extend([el] * cnt)

    # base indices (all except the vacancy atom)
    base_idx = [i for i in range(len(elem_list)) if i != idx_vac]

    # position of hopping atom in the base list
    hop_in_base = base_idx.index(idx_hop)

    pos_A = p["positions"][idx_hop].copy()
    pos_B = p["positions"][idx_vac].copy()
    diff  = pos_B - pos_A
    diff -= np.round(diff)   # shortest PBC path

    # rebuild element/count arrays after removing idx_vac
    new_counts = []
    cursor = 0
    for el, cnt in zip(p["elements"], p["counts"]):
        indices_in_group = list(range(cursor, cursor+cnt))
        kept = [i for i in indices_in_group if i in base_idx]
        new_counts.append(len(kept))
        cursor += cnt

    base_positions = p["positions"][base_idx].copy()

    images = []
    for k in range(n_images + 2):            # 00 … N+1
        t = k / (n_images + 1)
        img_pos = base_positions.copy()
        img_pos[hop_in_base] = (pos_A + t * diff) % 1.0
        images.append(dict(
            lattice=p["lattice"],
            elements=p["elements"],
            counts=new_counts,
            positions=img_pos,
        ))
    return images


# ── File templates ────────────────────────────────────────────────────────
NEB_INCAR = """\
SYSTEM = NEB

! Electronic structure
ENCUT  = 400
EDIFF  = 1E-05
PREC   = Normal
ALGO   = Fast
LREAL  = Auto
ISYM   = 0

! NEB settings
IBRION = 3
POTIM  = 0
NSW    = 300
EDIFFG = -0.05
IMAGES = {n_images}
SPRING = -5
LCLIMB = .TRUE.

! Output
LWAVE  = .FALSE.
LCHARG = .FALSE.

! Parallelisation (NPAR = IMAGES for balanced load)
NPAR   = {n_images}
"""

KPOINTS = """\
Gamma-point only
 0
Gamma
 1 1 1
 0 0 0
"""


def submit_sh(job_name, n_images):
    # Each image uses 8 MPI tasks → total = n_images × 8
    ntasks = n_images * 8
    return f"""\
#!/bin/bash -l
#SBATCH --job-name={job_name}_neb
#SBATCH --account=membrane
#SBATCH --partition=tier3

#SBATCH --mail-user=slack:@qhteme

#SBATCH -t 10:00:00

#SBATCH --ntasks={ntasks}
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=500M

spack load vasp@6.3.2
srun vasp_std
"""


# ── Main ──────────────────────────────────────────────────────────────────
for phase_dir_name, cfg in PHASES.items():
    phase_dir = CALC_DIR / phase_dir_name
    contcar   = phase_dir / "relax" / "CONTCAR"
    potcar    = phase_dir / "relax" / "POTCAR"

    if not contcar.exists():
        print(f"WARN  missing CONTCAR: {phase_dir_name}")
        continue

    p  = read_poscar(contcar)
    sc = cfg["supercell"]
    sp = make_supercell(p, sc)

    n_atoms = sum(sp["counts"])
    n_na    = sp["counts"][sp["elements"].index("Na")]

    idx_i, idx_j, dist = find_min_na_hop(sp)
    print(f"{phase_dir_name}")
    print(f"  supercell {sc} → {n_atoms} atoms, {n_na} Na")
    print(f"  shortest Na-Na hop: atom {idx_i} → {idx_j}  ({dist:.3f} Å)")

    images = build_neb_images(sp, idx_i, idx_j, N_IMAGES)

    # Write files
    neb_dir = phase_dir / "neb" / "hop1"
    neb_dir.mkdir(parents=True, exist_ok=True)

    for k, img in enumerate(images):
        img_dir = neb_dir / f"{k:02d}"
        img_dir.mkdir(exist_ok=True)
        write_poscar(img, img_dir / "POSCAR",
                     comment=f"{phase_dir_name} NEB image {k:02d}")

    (neb_dir / "INCAR").write_text(NEB_INCAR.format(n_images=N_IMAGES))
    (neb_dir / "KPOINTS").write_text(KPOINTS)
    (neb_dir / "submit.sh").write_text(submit_sh(phase_dir_name.split("_")[0], N_IMAGES))

    if potcar.exists():
        shutil.copy(potcar, neb_dir / "POTCAR")

    print(f"  → {neb_dir.relative_to(CALC_DIR)}  "
          f"({N_IMAGES+2} images: 00…{N_IMAGES+1:02d})")

print(f"\nDone. Submit with: python batch_submit_neb.py")
print("NOTE: endpoints (00 and 06) should ideally be pre-relaxed with the")
print("      vacancy present before launching the NEB. See batch_relax_neb_endpoints.py")
