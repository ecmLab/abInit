#!/usr/bin/env python3
"""
HCrO4- Ion Analysis Script
Analyzes axial (x) and radial (sqrt(y^2+z^2)) positions of HCrO4- ions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def read_hcro4_positions(dump_file, max_frames=None):
    """Read HCrO4- ion positions from LAMMPS dump file"""
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
        if max_frames and frame_count >= max_frames:
            break

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
        else:
            i += 1

    print(f"Read {len(timesteps)} timesteps")
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

    for t in timesteps:
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

def plot_hcro4_analysis(df_list, labels, output_dir):
    """Plot HCrO4- analysis comparing different simulations"""
    os.makedirs(output_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('HCrO4⁻ Ion Position Analysis', fontsize=16)

    colors = ['blue', 'red', 'green', 'orange', 'purple']

    # Plot 1: Axial position vs time
    for i, (df, label) in enumerate(zip(df_list, labels)):
        color = colors[i % len(colors)]
        axes[0,0].plot(df['time_ps'], df['axial_mean'], color=color, linewidth=2, label=label)
        axes[0,0].fill_between(df['time_ps'],
                              df['axial_mean'] - df['axial_std'],
                              df['axial_mean'] + df['axial_std'],
                              alpha=0.3, color=color)

    axes[0,0].set_xlabel('Time (ps)')
    axes[0,0].set_ylabel('Axial Position (Å)')
    axes[0,0].set_title('Axial Position vs Time')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)

    # Plot 2: Radial position vs time
    for i, (df, label) in enumerate(zip(df_list, labels)):
        color = colors[i % len(colors)]
        axes[0,1].plot(df['time_ps'], df['radial_mean'], color=color, linewidth=2, label=label)
        axes[0,1].fill_between(df['time_ps'],
                              df['radial_mean'] - df['radial_std'],
                              df['radial_mean'] + df['radial_std'],
                              alpha=0.3, color=color)

    axes[0,1].set_xlabel('Time (ps)')
    axes[0,1].set_ylabel('Radial Distance (Å)')
    axes[0,1].set_title('Radial Position vs Time')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)

    # Plot 3: Number of ions vs time
    for i, (df, label) in enumerate(zip(df_list, labels)):
        color = colors[i % len(colors)]
        axes[1,0].plot(df['time_ps'], df['n_ions'], color=color, linewidth=2, label=label, marker='o', markersize=3)

    axes[1,0].set_xlabel('Time (ps)')
    axes[1,0].set_ylabel('Number of HCrO4⁻ Ions')
    axes[1,0].set_title('Ion Count vs Time')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)

    # Plot 4: Axial displacement (relative to initial position)
    for i, (df, label) in enumerate(zip(df_list, labels)):
        if len(df) > 0:
            color = colors[i % len(colors)]
            initial_pos = df['axial_mean'].iloc[0]
            displacement = df['axial_mean'] - initial_pos
            axes[1,1].plot(df['time_ps'], displacement, color=color, linewidth=2, label=label)

    axes[1,1].set_xlabel('Time (ps)')
    axes[1,1].set_ylabel('Axial Displacement (Å)')
    axes[1,1].set_title('Axial Displacement vs Time')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/hcro4_analysis_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    # Simulation directories to analyze
    simulations = ['r10_1', 'r10_3']

    # Check which dump files exist
    available_sims = []
    for sim in simulations:
        dump_file = f'../{sim}/dump.lammpstrj'
        if os.path.exists(dump_file):
            available_sims.append(sim)
        else:
            print(f"Warning: {dump_file} not found, skipping...")

    if not available_sims:
        print("No dump files found!")
        return

    # Analyze each simulation
    dataframes = []
    labels = []

    for sim in available_sims:
        print(f"\n=== Analyzing {sim} ===")
        dump_file = f'../{sim}/dump.lammpstrj'

        # Read data (limit to first 20 frames for faster processing)
        data, timesteps = read_hcro4_positions(dump_file, max_frames=20)

        # Calculate statistics
        df = calculate_hcro4_statistics(data, timesteps)

        if not df.empty:
            # Save results
            df.to_csv(f'hcro4_positions_{sim}.csv', index=False)
            dataframes.append(df)
            labels.append(sim)

            print(f"Results saved to hcro4_positions_{sim}.csv")
            print(f"Average axial position: {df['axial_mean'].mean():.2f} ± {df['axial_std'].mean():.2f} Å")
            print(f"Average radial distance: {df['radial_mean'].mean():.2f} ± {df['radial_std'].mean():.2f} Å")
        else:
            print(f"No data found for {sim}")

    # Generate comparison plots
    if dataframes:
        print(f"\nGenerating comparison plots...")
        plot_hcro4_analysis(dataframes, labels, '.')
        print("Analysis complete! Check hcro4_analysis_comparison.png")

if __name__ == "__main__":
    main()