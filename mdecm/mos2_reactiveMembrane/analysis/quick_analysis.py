#!/usr/bin/env python3
"""
Quick LAMMPS Trajectory Analysis Script
Analyzes first 1000 timesteps for faster processing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def read_dump_quick(dump_file, max_timesteps=10):
    """Read first few timesteps from LAMMPS dump file"""
    print(f"Reading first {max_timesteps} timesteps from {dump_file}...")

    data = {}
    timesteps = []

    with open(dump_file, 'r') as f:
        lines = f.readlines()

    i = 0
    timestep_count = 0

    while i < len(lines) and timestep_count < max_timesteps:
        if 'ITEM: TIMESTEP' in lines[i]:
            timestep = int(lines[i+1].strip())
            timesteps.append(timestep)
            i += 2

            # Read number of atoms
            natoms = int(lines[i+1].strip())
            i += 6  # Skip box bounds

            # Read atom data
            atoms_data = []
            for j in range(natoms):
                atom_line = lines[i+1+j].split()
                atoms_data.append([
                    int(atom_line[0]),    # id
                    int(atom_line[1]),    # type
                    float(atom_line[3]),  # x
                    float(atom_line[4]),  # y
                    float(atom_line[5]),  # z
                ])

            data[timestep] = atoms_data
            i += natoms + 1
            timestep_count += 1
        else:
            i += 1

    return data, timesteps

def analyze_group_positions(data, timesteps, group_types):
    """Analyze positions for a specific group"""
    results = {
        'time': [],
        'x_mean': [], 'y_mean': [], 'z_mean': [],
        'x_std': [], 'y_std': [], 'z_std': []
    }

    for t in timesteps:
        # Filter atoms by group
        group_atoms = [atom for atom in data[t] if atom[1] in group_types]

        if group_atoms:
            positions = np.array([[atom[2], atom[3], atom[4]] for atom in group_atoms])

            results['time'].append(t * 0.001)  # convert to ps
            results['x_mean'].append(np.mean(positions[:, 0]))
            results['y_mean'].append(np.mean(positions[:, 1]))
            results['z_mean'].append(np.mean(positions[:, 2]))
            results['x_std'].append(np.std(positions[:, 0]))
            results['y_std'].append(np.std(positions[:, 1]))
            results['z_std'].append(np.std(positions[:, 2]))

    return pd.DataFrame(results)

def plot_comparison(data_r10_1, data_r10_3, group_name):
    """Plot comparison between r10_1 and r10_3"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'{group_name.upper()} Group Position Analysis', fontsize=16)

    # X position
    axes[0].plot(data_r10_1['time'], data_r10_1['x_mean'], 'b-', label='r10_1', linewidth=2)
    axes[0].plot(data_r10_3['time'], data_r10_3['x_mean'], 'r-', label='r10_3', linewidth=2)
    axes[0].fill_between(data_r10_1['time'],
                        data_r10_1['x_mean'] - data_r10_1['x_std'],
                        data_r10_1['x_mean'] + data_r10_1['x_std'], alpha=0.3, color='blue')
    axes[0].fill_between(data_r10_3['time'],
                        data_r10_3['x_mean'] - data_r10_3['x_std'],
                        data_r10_3['x_mean'] + data_r10_3['x_std'], alpha=0.3, color='red')
    axes[0].set_xlabel('Time (ps)')
    axes[0].set_ylabel('X Position (Å)')
    axes[0].set_title('X Position vs Time')
    axes[0].legend()
    axes[0].grid(True)

    # Y position
    axes[1].plot(data_r10_1['time'], data_r10_1['y_mean'], 'b-', label='r10_1', linewidth=2)
    axes[1].plot(data_r10_3['time'], data_r10_3['y_mean'], 'r-', label='r10_3', linewidth=2)
    axes[1].fill_between(data_r10_1['time'],
                        data_r10_1['y_mean'] - data_r10_1['y_std'],
                        data_r10_1['y_mean'] + data_r10_1['y_std'], alpha=0.3, color='blue')
    axes[1].fill_between(data_r10_3['time'],
                        data_r10_3['y_mean'] - data_r10_3['y_std'],
                        data_r10_3['y_mean'] + data_r10_3['y_std'], alpha=0.3, color='red')
    axes[1].set_xlabel('Time (ps)')
    axes[1].set_ylabel('Y Position (Å)')
    axes[1].set_title('Y Position vs Time')
    axes[1].legend()
    axes[1].grid(True)

    # Z position
    axes[2].plot(data_r10_1['time'], data_r10_1['z_mean'], 'b-', label='r10_1', linewidth=2)
    axes[2].plot(data_r10_3['time'], data_r10_3['z_mean'], 'r-', label='r10_3', linewidth=2)
    axes[2].fill_between(data_r10_1['time'],
                        data_r10_1['z_mean'] - data_r10_1['z_std'],
                        data_r10_1['z_mean'] + data_r10_1['z_std'], alpha=0.3, color='blue')
    axes[2].fill_between(data_r10_3['time'],
                        data_r10_3['z_mean'] - data_r10_3['z_std'],
                        data_r10_3['z_mean'] + data_r10_3['z_std'], alpha=0.3, color='red')
    axes[2].set_xlabel('Time (ps)')
    axes[2].set_ylabel('Z Position (Å)')
    axes[2].set_title('Z Position vs Time')
    axes[2].legend()
    axes[2].grid(True)

    plt.tight_layout()
    plt.savefig(f'comparison_{group_name}.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    # Group definitions
    groups = {
        'hcro': [2, 3, 10, 11, 12, 13],  # HCrO4- ions
        'water': [4, 5, 9],               # Water molecules
        'ions': [1, 8],                   # Cl- and Na+
    }

    # Read data from both simulations
    data_r10_1, timesteps_r10_1 = read_dump_quick('r10_1/dump.lammpstrj', max_timesteps=50)
    data_r10_3, timesteps_r10_3 = read_dump_quick('r10_3/dump.lammpstrj', max_timesteps=50)

    # Analyze each group
    for group_name, group_types in groups.items():
        print(f"\nAnalyzing {group_name} group...")

        # Analyze both simulations
        analysis_r10_1 = analyze_group_positions(data_r10_1, timesteps_r10_1, group_types)
        analysis_r10_3 = analyze_group_positions(data_r10_3, timesteps_r10_3, group_types)

        # Save results
        analysis_r10_1.to_csv(f'r10_1_quick_{group_name}.csv', index=False)
        analysis_r10_3.to_csv(f'r10_3_quick_{group_name}.csv', index=False)

        # Plot comparison
        plot_comparison(analysis_r10_1, analysis_r10_3, group_name)

        print(f"Analysis complete for {group_name}")

if __name__ == "__main__":
    main()