# DFT Plan: Na3SbS4 Hydration and Recyclability

This plan implements thermodynamic and transport calculations to support experiments on Na3SbS4 hydration/recovery and interfacial stabilization, aligned with Yaosen Tian et al. (Joule 2019) and Meng et al. (JES 2023).

## Scope
- Structures: Na3SbS4 (anhydrous) and Na3SbS4·8H2O (from Yaosen), with optional Na3SbS4·9H2O (modeled)
- Methods: VASP/PBE (+D3 optional), convex-hull reaction thermodynamics, gas mu(T,p) for H2O/H2S, NEB for Na migration, input scaffolding for reproducible runs

## Targets
- dE/dG for:
  - Na + Na3SbS4 -> 4 Na2S + Na3Sb (Reaction 1)
  - Na + Na3SbS4·8H2O -> 4 Na2S + Na3Sb + 16 NaH + 8 Na2O (Reaction 3)
  - Na + Na3SbS4·9H2O -> 4 Na2S + Na3Sb + 18 NaH + 9 Na2O (Reaction 2)
  - Hydration: Na3SbS4 + n H2O(g) -> Na3SbS4·nH2O (n = 8, 9)
- RH–T mapping for hydration favorability and dehydration recoverability at dry-room conditions
- Na migration barriers in Na3SbS4·8H2O vs 9H2O to rationalize conductivity trends

## Directory Map
- 10_structures: CIF inputs (copied from model/)
- 20_vasp_inputs: CIF->POSCAR converter, templates (INCAR, KPOINTS), and input prep
- 30_thermo: reactions, example energies, and dE/dG calculator with mu(T,p)
- 40_transport: NEB scaffolding

## DFT Settings (recommended)
- Code: VASP
- Functional: PBE; consider D3(BJ) (IVDW=12) for hydrates
- PAW: Na (Na or Na_pv), Sb, S, O, H
- ENCUT: 520 eV; PREC=Accurate; LASPH = .TRUE.
- Smearing: ISMEAR=0, SIGMA=0.05 eV
- Relax: IBRION=2, ISIF=3, NSW=100, EDIFF=1e-6, EDIFFG=-0.02 eV/Ang
- Static: NSW=0
- K-mesh: Gamma-centered, auto from ~0.2 A^-1 spacing (prepare_inputs.py)
- NEB: IMAGES=5, LCLIMB=.TRUE., IOPT=1

## Thermodynamics
- Solids: use static total energies (eV/f.u.) after relaxation
- Gases: mu(T,p) via Shomate equations (H2O, H2S) with ideal gas RT ln(p/p0)
- RH handling: p_H2O from RH and T using saturation vapor pressure fit
- Inputs: 30_thermo/energies.yaml and reactions.yaml

## Transport (NEB)
- Build supercells and identify Na vacancy hops
- Scaffold NEB folders (40_transport/make_neb.py)
- Use relaxed endpoints; interpolate images; run NEB

## Workflow
1) Prepare inputs
   - python3 20_vasp_inputs/prepare_inputs.py --cif 10_structures/Na3SbS4.cif --out 20_vasp_inputs/Na3SbS4
   - python3 20_vasp_inputs/prepare_inputs.py --cif 10_structures/Na3SbS4_yaosen.cif --out 20_vasp_inputs/Na3SbS4_8H2O
2) Relax -> Static for each structure and products (Na2S, Na3Sb, NaH, Na2O, Na)
3) Collect energies into 30_thermo/energies.yaml (see energies.example.yaml)
4) Compute dE/dG and RH–T trends
   - python3 30_thermo/compute_reaction_energies.py --energies 30_thermo/energies.yaml --reactions 30_thermo/reactions.yaml --T 298 --RH 0.68
5) NEB for Na migration in hydrates
   - python3 40_transport/make_neb.py --out 40_transport/neb_Na3SbS4_8H2O --images 5

## Notes & Assumptions
- CIF->POSCAR parser is minimal but tailored to provided CIFs
- POTCARs are user-supplied (not included)
- Gas mu uses 298–1000 K Shomate fits; sufficient for room-temperature to moderate anneal regimes
- For interface-specific modeling (Na/H2O at slab), extend with surface slab builders (pymatgen/ase) if available

## References
- Y. Tian et al., Joule 3, 1037–1050 (2019)
- Y.-T. Chen et al., J. Electrochem. Soc. 170, 080521 (2023)
- Y. Zhu & Y. Mo, Angew. Chem. Int. Ed. 59, 17472–17476 (2020)
