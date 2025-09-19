#!/usr/bin/env python3
"""
HCrO4- Ion Focused Analysis Script
Plots only axial and radial displacement using full simulation data
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def read_hcro4_positions(dump_file):
    """Read ALL HCrO4- ion positions from LAMMPS dump file"""
    print(f"Reading HCrO4- positions from {dump_file}...")

    # HCrO4- atom types: 2 (Cr), 3 (H), 10-13 (O1-O4)
    hcro4_types = [2, 3, 10, 11, 12, 13]

    data = {}
    timesteps = []

    with open(dump_file, 'r') as f:
        lines = f.readlines()

    i = 0
    frame_count = 0

    while i < len(lines):
        if 'ITEM: TIMESTEP' in lines[i]:
            timestep = int(lines[i+1].strip())
            timesteps.append(timestep)
            i += 2

            # Read number of atoms
            natoms = int(lines[i+1].strip())
            i += 6  # Skip box bounds and header

            # Read atom data and filter for HCrO4-
            hcro4_atoms = []
            for j in range(natoms):
                atom_line = lines[i+1+j].split()
                atom_type = int(atom_line[1])

                if atom_type in hcro4_types:
                    hcro4_atoms.append({
                        'id': int(atom_line[0]),
                        'type': atom_type,
                        'mol': int(atom_line[2]),
                        'x': float(atom_line[3]),
                        'y': float(atom_line[4]),
                        'z': float(atom_line[5])
                    })

            data[timestep] = hcro4_atoms
            i += natoms + 1
            frame_count += 1

            if frame_count % 100 == 0:
                print(f"Read {frame_count} frames...")

        else:
            i += 1

    print(f"Read {len(timesteps)} timesteps total")
    return data, timesteps

def calculate_hcro4_statistics(data, timesteps):
    """Calculate axial and radial position statistics for HCrO4- ions"""
    results = {
        'time_ps': [],
        'axial_mean': [],      # x position mean
        'axial_std': [],       # x position std
        'radial_mean': [],     # sqrt(y^2+z^2) mean
        'radial_std': [],      # sqrt(y^2+z^2) std
        'n_ions': []           # number of HCrO4- ions
    }

    print("Calculating statistics...")

    for idx, t in enumerate(timesteps):
        if idx % 100 == 0:
            print(f"Processing frame {idx}/{len(timesteps)}")

        atoms = data[t]
        if not atoms:
            continue

        # Group atoms by molecule ID to get complete HCrO4- ions
        molecules = {}
        for atom in atoms:
            mol_id = atom['mol']
            if mol_id not in molecules:
                molecules[mol_id] = []
            molecules[mol_id].append(atom)

        # Calculate center of mass for each HCrO4- ion
        ion_positions = []
        for mol_id, mol_atoms in molecules.items():
            if len(mol_atoms) >= 5:  # Complete HCrO4- should have 6 atoms, but allow some flexibility
                # Calculate center of mass
                x_com = np.mean([atom['x'] for atom in mol_atoms])
                y_com = np.mean([atom['y'] for atom in mol_atoms])
                z_com = np.mean([atom['z'] for atom in mol_atoms])

                # Calculate radial distance from tube axis (assuming tube along x)
                radial_dist = np.sqrt(y_com**2 + z_com**2)

                ion_positions.append({
                    'axial': x_com,
                    'radial': radial_dist
                })

        if ion_positions:
            axial_positions = [pos['axial'] for pos in ion_positions]
            radial_positions = [pos['radial'] for pos in ion_positions]

            results['time_ps'].append(t * 0.001)  # convert fs to ps
            results['axial_mean'].append(np.mean(axial_positions))
            results['axial_std'].append(np.std(axial_positions))
            results['radial_mean'].append(np.mean(radial_positions))
            results['radial_std'].append(np.std(radial_positions))
            results['n_ions'].append(len(ion_positions))

    return pd.DataFrame(results)

def plot_displacement_comparison(df_r10_1, df_r10_3, output_dir):
    """Plot axial and radial displacement comparison"""
    os.makedirs(output_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('HCrO4⁻ Ion Displacement Analysis', fontsize=16)

    # Calculate initial positions for displacement
    if len(df_r10_1) > 0:
        initial_axial_r10_1 = df_r10_1['axial_mean'].iloc[0]
        initial_radial_r10_1 = df_r10_1['radial_mean'].iloc[0]
        axial_disp_r10_1 = df_r10_1['axial_mean'] - initial_axial_r10_1
        radial_disp_r10_1 = df_r10_1['radial_mean'] - initial_radial_r10_1

    if len(df_r10_3) > 0:
        initial_axial_r10_3 = df_r10_3['axial_mean'].iloc[0]
        initial_radial_r10_3 = df_r10_3['radial_mean'].iloc[0]
        axial_disp_r10_3 = df_r10_3['axial_mean'] - initial_axial_r10_3
        radial_disp_r10_3 = df_r10_3['radial_mean'] - initial_radial_r10_3

    # Plot 1: Axial displacement
    if len(df_r10_1) > 0:
        axes[0].plot(df_r10_1['time_ps'], axial_disp_r10_1, 'blue', linewidth=2, label='r10_1')
        axes[0].fill_between(df_r10_1['time_ps'],
                            axial_disp_r10_1 - df_r10_1['axial_std'],
                            axial_disp_r10_1 + df_r10_1['axial_std'],
                            alpha=0.3, color='blue')

    if len(df_r10_3) > 0:
        axes[0].plot(df_r10_3['time_ps'], axial_disp_r10_3, 'red', linewidth=2, label='r10_3')
        axes[0].fill_between(df_r10_3['time_ps'],
                            axial_disp_r10_3 - df_r10_3['axial_std'],
                            axial_disp_r10_3 + df_r10_3['axial_std'],
                            alpha=0.3, color='red')

    axes[0].set_xlabel('Time (ps)', fontsize=12)
    axes[0].set_ylabel('Axial Displacement (Å)', fontsize=12)
    axes[0].set_title('Axial Displacement vs Time', fontsize=14)
    axes[0].legend(fontsize=12)
    axes[0].grid(True, alpha=0.3)

    # Plot 2: Radial displacement
    if len(df_r10_1) > 0:
        axes[1].plot(df_r10_1['time_ps'], radial_disp_r10_1, 'blue', linewidth=2, label='r10_1')
        axes[1].fill_between(df_r10_1['time_ps'],
                            radial_disp_r10_1 - df_r10_1['radial_std'],
                            radial_disp_r10_1 + df_r10_1['radial_std'],
                            alpha=0.3, color='blue')

    if len(df_r10_3) > 0:
        axes[1].plot(df_r10_3['time_ps'], radial_disp_r10_3, 'red', linewidth=2, label='r10_3')
        axes[1].fill_between(df_r10_3['time_ps'],
                            radial_disp_r10_3 - df_r10_3['radial_std'],
                            radial_disp_r10_3 + df_r10_3['radial_std'],
                            alpha=0.3, color='red')

    axes[1].set_xlabel('Time (ps)', fontsize=12)
    axes[1].set_ylabel('Radial Displacement (Å)', fontsize=12)
    axes[1].set_title('Radial Displacement vs Time', fontsize=14)
    axes[1].legend(fontsize=12)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/hcro4_displacement_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    # Analyze r10_1 and r10_3
    simulations = ['r10_1', 'r10_3']
    dataframes = {}

    for sim in simulations:
        dump_file = f'../{sim}/dump.lammpstrj'
        if os.path.exists(dump_file):
            print(f"\n=== Analyzing {sim} ===")

            # Read ALL data (no frame limit)
            data, timesteps = read_hcro4_positions(dump_file)

            # Calculate statistics
            df = calculate_hcro4_statistics(data, timesteps)

            if not df.empty:
                # Save results
                df.to_csv(f'hcro4_full_{sim}.csv', index=False)
                dataframes[sim] = df

                print(f"Results saved to hcro4_full_{sim}.csv")
                print(f"Total simulation time: {df['time_ps'].max():.1f} ps")
                print(f"Final axial displacement: {df['axial_mean'].iloc[-1] - df['axial_mean'].iloc[0]:.2f} Å")
                print(f"Final radial displacement: {df['radial_mean'].iloc[-1] - df['radial_mean'].iloc[0]:.2f} Å")
            else:
                print(f"No data found for {sim}")
        else:
            print(f"File ../{sim}/dump.lammpstrj not found")

    # Generate displacement plots
    if 'r10_1' in dataframes or 'r10_3' in dataframes:
        print(f"\nGenerating displacement comparison plots...")
        df_r10_1 = dataframes.get('r10_1', pd.DataFrame())
        df_r10_3 = dataframes.get('r10_3', pd.DataFrame())

        plot_displacement_comparison(df_r10_1, df_r10_3, '.')
        print("Analysis complete! Check hcro4_displacement_comparison.png")

if __name__ == "__main__":
    main()