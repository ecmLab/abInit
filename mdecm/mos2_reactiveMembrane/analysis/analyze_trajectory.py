#!/usr/bin/env python3
"""
LAMMPS Trajectory Analysis Script
Analyzes position and velocity data from LAMMPS dump files
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from collections import defaultdict

class LAMMPSAnalyzer:
    def __init__(self, dump_file):
        self.dump_file = dump_file
        self.data = defaultdict(list)
        self.timesteps = []
        self.box_bounds = {}

        # Atom type mapping (based on your lmp.in)
        self.atom_types = {
            1: 'Cl', 2: 'Cr', 3: 'H', 4: 'H1', 5: 'H2',
            6: 'Mo1', 7: 'Mo2', 8: 'Na', 9: 'OW',
            10: 'O1', 11: 'O2', 12: 'O3', 13: 'O4',
            14: 'S1', 15: 'S2', 16: 'S3', 17: 'S4'
        }

        # Group definitions
        self.groups = {
            'water': [4, 5, 9],
            'mos2': [6, 7, 14, 15, 16, 17],
            'hcro': [2, 3, 10, 11, 12, 13],
            'ions': [1, 8],  # Cl- and Na+
            'all_mobile': [1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13]
        }

    def read_dump_file(self):
        """Read LAMMPS dump file and extract data"""
        print(f"Reading {self.dump_file}...")

        with open(self.dump_file, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            if 'ITEM: TIMESTEP' in lines[i]:
                timestep = int(lines[i+1].strip())
                self.timesteps.append(timestep)
                i += 2

                # Read number of atoms
                natoms = int(lines[i+1].strip())
                i += 2

                # Read box bounds
                box_data = []
                for j in range(3):
                    box_data.append([float(x) for x in lines[i+1+j].split()])
                self.box_bounds[timestep] = box_data
                i += 4

                # Read atom data
                atoms_data = []
                for j in range(natoms):
                    atom_line = lines[i+1+j].split()
                    atoms_data.append([
                        int(atom_line[0]),    # id
                        int(atom_line[1]),    # type
                        int(atom_line[2]),    # mol
                        float(atom_line[3]),  # x
                        float(atom_line[4]),  # y
                        float(atom_line[5]),  # z
                        float(atom_line[6])   # charge
                    ])

                self.data[timestep] = atoms_data
                i += natoms + 1
            else:
                i += 1

        print(f"Read {len(self.timesteps)} timesteps")

    def calculate_velocities(self):
        """Calculate velocities from position differences"""
        self.velocities = {}
        dt = 1.0  # fs (from your timestep)

        for i in range(len(self.timesteps)-1):
            t1, t2 = self.timesteps[i], self.timesteps[i+1]
            dt_actual = (t2 - t1) * dt / 1000.0  # convert to ps

            vel_data = []
            data1 = {atom[0]: atom for atom in self.data[t1]}
            data2 = {atom[0]: atom for atom in self.data[t2]}

            for atom_id in data1:
                if atom_id in data2:
                    atom1, atom2 = data1[atom_id], data2[atom_id]
                    vx = (atom2[3] - atom1[3]) / dt_actual
                    vy = (atom2[4] - atom1[4]) / dt_actual
                    vz = (atom2[5] - atom1[5]) / dt_actual

                    vel_data.append([
                        atom_id, atom1[1], atom1[2],  # id, type, mol
                        vx, vy, vz,  # velocities
                        np.sqrt(vx**2 + vy**2 + vz**2)  # speed
                    ])

            self.velocities[t1] = vel_data

    def filter_by_group(self, data, group_name):
        """Filter atoms by group"""
        if group_name not in self.groups:
            return data

        group_types = self.groups[group_name]
        return [atom for atom in data if atom[1] in group_types]

    def analyze_positions(self, group='all_mobile'):
        """Analyze position statistics"""
        results = {
            'time': [],
            'x_mean': [], 'x_std': [],
            'y_mean': [], 'y_std': [],
            'z_mean': [], 'z_std': [],
            'x_com': [], 'y_com': [], 'z_com': []  # center of mass
        }

        for t in self.timesteps:
            data = self.filter_by_group(self.data[t], group)
            if not data:
                continue

            positions = np.array([[atom[3], atom[4], atom[5]] for atom in data])

            results['time'].append(t * 0.001)  # convert to ps
            results['x_mean'].append(np.mean(positions[:, 0]))
            results['x_std'].append(np.std(positions[:, 0]))
            results['y_mean'].append(np.mean(positions[:, 1]))
            results['y_std'].append(np.std(positions[:, 1]))
            results['z_mean'].append(np.mean(positions[:, 2]))
            results['z_std'].append(np.std(positions[:, 2]))

            # Center of mass
            results['x_com'].append(np.mean(positions[:, 0]))
            results['y_com'].append(np.mean(positions[:, 1]))
            results['z_com'].append(np.mean(positions[:, 2]))

        return pd.DataFrame(results)

    def analyze_velocities(self, group='all_mobile'):
        """Analyze velocity statistics"""
        results = {
            'time': [],
            'vx_mean': [], 'vx_std': [],
            'vy_mean': [], 'vy_std': [],
            'vz_mean': [], 'vz_std': [],
            'speed_mean': [], 'speed_std': []
        }

        for t in sorted(self.velocities.keys()):
            data = self.filter_by_group(self.velocities[t], group)
            if not data:
                continue

            velocities = np.array([[atom[3], atom[4], atom[5], atom[6]] for atom in data])

            results['time'].append(t * 0.001)  # convert to ps
            results['vx_mean'].append(np.mean(velocities[:, 0]))
            results['vx_std'].append(np.std(velocities[:, 0]))
            results['vy_mean'].append(np.mean(velocities[:, 1]))
            results['vy_std'].append(np.std(velocities[:, 1]))
            results['vz_mean'].append(np.mean(velocities[:, 2]))
            results['vz_std'].append(np.std(velocities[:, 2]))
            results['speed_mean'].append(np.mean(velocities[:, 3]))
            results['speed_std'].append(np.std(velocities[:, 3]))

        return pd.DataFrame(results)

    def plot_analysis(self, pos_data, vel_data, group_name, save_dir):
        """Generate analysis plots"""
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'LAMMPS Analysis - {group_name.upper()} Group', fontsize=16)

        # Position plots
        axes[0,0].plot(pos_data['time'], pos_data['x_mean'], 'b-', label='X')
        axes[0,0].fill_between(pos_data['time'],
                              pos_data['x_mean'] - pos_data['x_std'],
                              pos_data['x_mean'] + pos_data['x_std'], alpha=0.3)
        axes[0,0].set_xlabel('Time (ps)')
        axes[0,0].set_ylabel('X Position (Å)')
        axes[0,0].set_title('X Position vs Time')
        axes[0,0].grid(True)

        axes[0,1].plot(pos_data['time'], pos_data['y_mean'], 'g-', label='Y')
        axes[0,1].fill_between(pos_data['time'],
                              pos_data['y_mean'] - pos_data['y_std'],
                              pos_data['y_mean'] + pos_data['y_std'], alpha=0.3)
        axes[0,1].set_xlabel('Time (ps)')
        axes[0,1].set_ylabel('Y Position (Å)')
        axes[0,1].set_title('Y Position vs Time')
        axes[0,1].grid(True)

        axes[0,2].plot(pos_data['time'], pos_data['z_mean'], 'r-', label='Z')
        axes[0,2].fill_between(pos_data['time'],
                              pos_data['z_mean'] - pos_data['z_std'],
                              pos_data['z_mean'] + pos_data['z_std'], alpha=0.3)
        axes[0,2].set_xlabel('Time (ps)')
        axes[0,2].set_ylabel('Z Position (Å)')
        axes[0,2].set_title('Z Position vs Time')
        axes[0,2].grid(True)

        # Velocity plots
        if not vel_data.empty:
            axes[1,0].plot(vel_data['time'], vel_data['vx_mean'], 'b-', label='Vx')
            axes[1,0].fill_between(vel_data['time'],
                                  vel_data['vx_mean'] - vel_data['vx_std'],
                                  vel_data['vx_mean'] + vel_data['vx_std'], alpha=0.3)
            axes[1,0].set_xlabel('Time (ps)')
            axes[1,0].set_ylabel('X Velocity (Å/ps)')
            axes[1,0].set_title('X Velocity vs Time')
            axes[1,0].grid(True)

            axes[1,1].plot(vel_data['time'], vel_data['vy_mean'], 'g-', label='Vy')
            axes[1,1].fill_between(vel_data['time'],
                                  vel_data['vy_mean'] - vel_data['vy_std'],
                                  vel_data['vy_mean'] + vel_data['vy_std'], alpha=0.3)
            axes[1,1].set_xlabel('Time (ps)')
            axes[1,1].set_ylabel('Y Velocity (Å/ps)')
            axes[1,1].set_title('Y Velocity vs Time')
            axes[1,1].grid(True)

            axes[1,2].plot(vel_data['time'], vel_data['speed_mean'], 'purple', label='Speed')
            axes[1,2].fill_between(vel_data['time'],
                                  vel_data['speed_mean'] - vel_data['speed_std'],
                                  vel_data['speed_mean'] + vel_data['speed_std'], alpha=0.3)
            axes[1,2].set_xlabel('Time (ps)')
            axes[1,2].set_ylabel('Speed (Å/ps)')
            axes[1,2].set_title('Speed vs Time')
            axes[1,2].grid(True)

        plt.tight_layout()
        plt.savefig(f'{save_dir}/analysis_{group_name}.png', dpi=300, bbox_inches='tight')
        plt.show()

def main():
    # Analyze both r10_1 and r10_3 simulations
    simulations = ['r10_1', 'r10_3']

    for sim in simulations:
        dump_file = f'{sim}/dump.lammpstrj'
        if not os.path.exists(dump_file):
            print(f"File {dump_file} not found, skipping...")
            continue

        print(f"\n=== Analyzing {sim} ===")
        analyzer = LAMMPSAnalyzer(dump_file)

        # Read data
        analyzer.read_dump_file()
        analyzer.calculate_velocities()

        # Create output directory
        os.makedirs(f'{sim}/analysis', exist_ok=True)

        # Analyze different groups
        groups_to_analyze = ['hcro', 'water', 'ions', 'all_mobile']

        for group in groups_to_analyze:
            print(f"Analyzing {group} group...")

            pos_data = analyzer.analyze_positions(group)
            vel_data = analyzer.analyze_velocities(group)

            # Save data
            pos_data.to_csv(f'{sim}/analysis/positions_{group}.csv', index=False)
            vel_data.to_csv(f'{sim}/analysis/velocities_{group}.csv', index=False)

            # Generate plots
            analyzer.plot_analysis(pos_data, vel_data, group, f'{sim}/analysis')

            print(f"Results saved to {sim}/analysis/")

if __name__ == "__main__":
    main()