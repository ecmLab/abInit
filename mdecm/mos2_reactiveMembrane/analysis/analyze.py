#!/usr/bin/env python3
"""
HCrO4- Ion Analysis for r10_3 and r10_5 with correct tube center
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def unwrap_coordinate(x_prev, x_curr, box_length):
    """Correctly unwrap coordinate accounting for periodic boundaries"""
    dx = x_curr - x_prev
    if dx > box_length / 2.0:
        dx -= box_length
    elif dx < -box_length / 2.0:
        dx += box_length
    return dx

def determine_tube_center(dump_file):
    """Determine the actual tube center from box bounds"""
    print(f"Determining tube center from {dump_file}...")

    with open(dump_file, 'r') as f:
        lines = f.readlines()

    # Find first box bounds
    for i, line in enumerate(lines):
        if 'ITEM: BOX BOUNDS' in line:
            # Read box bounds for Y and Z
            y_bounds = [float(x) for x in lines[i+2].split()]  # Y bounds
            z_bounds = [float(x) for x in lines[i+3].split()]  # Z bounds

            # Calculate center
            y_center = (y_bounds[0] + y_bounds[1]) / 2.0
            z_center = (z_bounds[0] + z_bounds[1]) / 2.0

            print(f"Y bounds: {y_bounds[0]:.3f} to {y_bounds[1]:.3f}, Center: {y_center:.3f}")
            print(f"Z bounds: {z_bounds[0]:.3f} to {z_bounds[1]:.3f}, Center: {z_center:.3f}")

            return y_center, z_center

    # Fallback
    return 17.548, 17.548

def read_dump_with_unwrapping(dump_file):
    """Read dump file and apply unwrapping frame by frame"""
    print(f"Reading and unwrapping {dump_file}...")

    hcro4_types = [2, 3, 10, 11, 12, 13]  # HCrO4- atom types
    data = {}
    timesteps = []
    box_bounds = {}

    with open(dump_file, 'r') as f:
        lines = f.readlines()

    i = 0
    frame_count = 0

    while i < len(lines):
        if 'ITEM: TIMESTEP' in lines[i]:
            timestep = int(lines[i+1].strip())
            timesteps.append(timestep)
            i += 2

            natoms = int(lines[i+1].strip())
            i += 2

            # Read box bounds
            box_data = []
            for j in range(3):
                bounds = [float(x) for x in lines[i+1+j].split()]
                box_data.append(bounds)
            box_bounds[timestep] = box_data
            i += 4

            # Read atom data
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

            if frame_count % 1000 == 0:
                print(f"Read {frame_count} frames...")

        else:
            i += 1

    print(f"Read {len(timesteps)} frames total")

    # Now unwrap trajectories
    return unwrap_trajectories_fixed(data, timesteps, box_bounds)

def unwrap_trajectories_fixed(data, timesteps, box_bounds):
    """Fixed unwrapping algorithm that tracks each molecule properly"""
    print("Unwrapping trajectories...")

    # Get box dimensions
    first_box = box_bounds[timesteps[0]]
    box_lengths = [bounds[1] - bounds[0] for bounds in first_box]
    print(f"Box lengths: X={box_lengths[0]:.2f}, Y={box_lengths[1]:.2f}, Z={box_lengths[2]:.2f} Å")

    # Group atoms by molecule
    molecules = {}
    for t in timesteps:
        for atom in data[t]:
            mol_id = atom['mol']
            if mol_id not in molecules:
                molecules[mol_id] = {}
            if t not in molecules[mol_id]:
                molecules[mol_id][t] = []
            molecules[mol_id][t].append(atom)

    # Unwrap each molecule's trajectory
    unwrapped_data = {}

    for mol_id in molecules:
        mol_timesteps = sorted(molecules[mol_id].keys())

        if len(mol_timesteps) < 2:
            continue

        # Track unwrapped center of mass for this molecule
        unwrapped_com = {}

        for i, t in enumerate(mol_timesteps):
            mol_atoms = molecules[mol_id][t]

            if len(mol_atoms) >= 5:  # Complete HCrO4- ion
                # Calculate wrapped center of mass
                x_com = np.mean([atom['x'] for atom in mol_atoms])
                y_com = np.mean([atom['y'] for atom in mol_atoms])
                z_com = np.mean([atom['z'] for atom in mol_atoms])

                if i == 0:
                    # First frame - no unwrapping needed
                    unwrapped_x = x_com
                    unwrapped_y = y_com
                    unwrapped_z = z_com
                else:
                    # Unwrap relative to previous frame
                    prev_t = mol_timesteps[i-1]
                    if prev_t in unwrapped_com:
                        prev_unwrapped = unwrapped_com[prev_t]

                        # Calculate unwrapped displacement
                        dx = unwrap_coordinate(prev_unwrapped[0], x_com, box_lengths[0])
                        dy = unwrap_coordinate(prev_unwrapped[1], y_com, box_lengths[1])
                        dz = unwrap_coordinate(prev_unwrapped[2], z_com, box_lengths[2])

                        # Update unwrapped position
                        unwrapped_x = prev_unwrapped[0] + dx
                        unwrapped_y = prev_unwrapped[1] + dy
                        unwrapped_z = prev_unwrapped[2] + dz
                    else:
                        unwrapped_x = x_com
                        unwrapped_y = y_com
                        unwrapped_z = z_com

                unwrapped_com[t] = (unwrapped_x, unwrapped_y, unwrapped_z)

        # Store in output format
        for t in unwrapped_com:
            if t not in unwrapped_data:
                unwrapped_data[t] = []

            x, y, z = unwrapped_com[t]
            unwrapped_data[t].append({
                'mol_id': mol_id,
                'x': x,
                'y': y,
                'z': z
            })

    print(f"Successfully unwrapped {len(molecules)} molecules")
    return unwrapped_data, timesteps

def calculate_statistics_from_unwrapped(unwrapped_data, timesteps, tube_center):
    """Calculate statistics from unwrapped data with correct tube center"""
    results = {
        'time_ps': [],
        'axial_mean': [],
        'axial_std': [],
        'radial_mean': [],
        'radial_std': [],
        'axial_velocity_mean': [],
        'axial_velocity_std': [],
        'radial_velocity_mean': [],
        'radial_velocity_std': [],
        'n_ions': []
    }

    ycenter, zcenter = tube_center

    # Store previous positions for each ion
    prev_positions = {}
    prev_time = None

    for t in timesteps:
        if t in unwrapped_data and unwrapped_data[t]:
            ions = unwrapped_data[t]

            axial_positions = [ion['x'] for ion in ions]
            radial_positions = [np.sqrt((ion['y'] - ycenter)**2 + (ion['z'] - zcenter)**2) for ion in ions]

            current_time = t * 0.001  # convert to ps

            # Calculate individual ion velocities
            axial_velocities = []
            radial_velocities = []

            if prev_time is not None:
                dt = current_time - prev_time
                if dt > 0:
                    for ion in ions:
                        mol_id = ion['mol_id']
                        if mol_id in prev_positions:
                            prev_x, prev_r = prev_positions[mol_id]
                            curr_x = ion['x']
                            curr_r = np.sqrt((ion['y'] - ycenter)**2 + (ion['z'] - zcenter)**2)

                            axial_vel = (curr_x - prev_x) / dt
                            radial_vel = (curr_r - prev_r) / dt

                            axial_velocities.append(axial_vel)
                            radial_velocities.append(radial_vel)

            # Store current positions for next iteration
            prev_positions = {}
            for ion in ions:
                mol_id = ion['mol_id']
                curr_r = np.sqrt((ion['y'] - ycenter)**2 + (ion['z'] - zcenter)**2)
                prev_positions[mol_id] = (ion['x'], curr_r)

            # Calculate statistics
            if axial_velocities:
                axial_vel_mean = np.mean(axial_velocities)
                axial_vel_std = np.std(axial_velocities)
                radial_vel_mean = np.mean(radial_velocities)
                radial_vel_std = np.std(radial_velocities)
            else:
                axial_vel_mean = 0.0
                axial_vel_std = 0.0
                radial_vel_mean = 0.0
                radial_vel_std = 0.0

            results['time_ps'].append(current_time)
            results['axial_mean'].append(np.mean(axial_positions))
            results['axial_std'].append(np.std(axial_positions))
            results['radial_mean'].append(np.mean(radial_positions))
            results['radial_std'].append(np.std(radial_positions))
            results['axial_velocity_mean'].append(axial_vel_mean)
            results['axial_velocity_std'].append(axial_vel_std)
            results['radial_velocity_mean'].append(radial_vel_mean)
            results['radial_velocity_std'].append(radial_vel_std)
            results['n_ions'].append(len(ions))

            prev_time = current_time

    return pd.DataFrame(results)

def plot_displacement_comparison(dataframes):
    """Plot displacement and velocity comparison for all simulations"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('HCrO4⁻ Ion Analysis: r10_1, r10_2, r10_3, r10_5', fontsize=16)

    # Define colors for each simulation
    colors = {'r10_1': 'blue', 'r10_2': 'orange', 'r10_3': 'green', 'r10_5': 'red'}
    datasets = [(sim, df, colors[sim]) for sim, df in dataframes.items()]

    for name, df, color in datasets:
        if len(df) > 0:
            initial_axial = df['axial_mean'].iloc[0]
            initial_radial = df['radial_mean'].iloc[0]
            axial_disp = df['axial_mean'] - initial_axial
            radial_disp = df['radial_mean'] - initial_radial

            # Axial displacement
            axes[0,0].plot(df['time_ps'], axial_disp, color=color, linewidth=2, label=name)
            axes[0,0].fill_between(df['time_ps'],
                                axial_disp - df['axial_std'],
                                axial_disp + df['axial_std'],
                                alpha=0.3, color=color)

            # Radial displacement
            axes[0,1].plot(df['time_ps'], radial_disp, color=color, linewidth=2, label=name)
            axes[0,1].fill_between(df['time_ps'],
                                radial_disp - df['radial_std'],
                                radial_disp + df['radial_std'],
                                alpha=0.3, color=color)

            # Axial velocity with shadows (same style as displacement)
            axes[1,0].plot(df['time_ps'], df['axial_velocity_mean'], color=color, linewidth=2, label=name)
            axes[1,0].fill_between(df['time_ps'],
                                df['axial_velocity_mean'] - df['axial_velocity_std'],
                                df['axial_velocity_mean'] + df['axial_velocity_std'],
                                alpha=0.3, color=color)

            # Radial velocity with shadows (same style as displacement)
            axes[1,1].plot(df['time_ps'], df['radial_velocity_mean'], color=color, linewidth=2, label=name)
            axes[1,1].fill_between(df['time_ps'],
                                df['radial_velocity_mean'] - df['radial_velocity_std'],
                                df['radial_velocity_mean'] + df['radial_velocity_std'],
                                alpha=0.3, color=color)

    # Format displacement plots
    axes[0,0].set_xlabel('Time (ps)', fontsize=12)
    axes[0,0].set_ylabel('Axial Displacement (Å)', fontsize=12)
    axes[0,0].set_title('Axial Displacement vs Time', fontsize=14)
    axes[0,0].legend(fontsize=12)
    axes[0,0].grid(True, alpha=0.3)

    axes[0,1].set_xlabel('Time (ps)', fontsize=12)
    axes[0,1].set_ylabel('Radial Displacement (Å)', fontsize=12)
    axes[0,1].set_title('Radial Displacement vs Time', fontsize=14)
    axes[0,1].legend(fontsize=12)
    axes[0,1].grid(True, alpha=0.3)

    # Format velocity plots
    axes[1,0].set_xlabel('Time (ps)', fontsize=12)
    axes[1,0].set_ylabel('Axial Velocity (Å/ps)', fontsize=12)
    axes[1,0].set_title('Axial Velocity vs Time', fontsize=14)
    axes[1,0].legend(fontsize=12)
    axes[1,0].grid(True, alpha=0.3)

    axes[1,1].set_xlabel('Time (ps)', fontsize=12)
    axes[1,1].set_ylabel('Radial Velocity (Å/ps)', fontsize=12)
    axes[1,1].set_title('Radial Velocity vs Time', fontsize=14)
    axes[1,1].legend(fontsize=12)
    axes[1,1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('r10_all_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    simulations = ['r10_1', 'r10_2', 'r10_3', 'r10_5']
    results = {}
    dataframes = {}

    for sim in simulations:
        dump_file = f'../{sim}/dump.lammpstrj'
        if os.path.exists(dump_file):
            print(f"\n=== Analyzing {sim} ===")

            # Determine tube center for this simulation
            tube_center = determine_tube_center(dump_file)
            print(f"Using tube center: Y={tube_center[0]:.3f}, Z={tube_center[1]:.3f}")

            # Read and unwrap data
            unwrapped_data, timesteps = read_dump_with_unwrapping(dump_file)

            # Calculate statistics
            df = calculate_statistics_from_unwrapped(unwrapped_data, timesteps, tube_center)

            if not df.empty:
                dataframes[sim] = df

                # Calculate final displacement and average velocities
                initial_axial = df['axial_mean'].iloc[0]
                final_axial = df['axial_mean'].iloc[-1]
                axial_displacement = final_axial - initial_axial

                initial_radial = df['radial_mean'].iloc[0]
                final_radial = df['radial_mean'].iloc[-1]
                radial_displacement = final_radial - initial_radial

                # Calculate average velocities (excluding first frame with zero velocity)
                avg_axial_velocity = df['axial_velocity_mean'].iloc[1:].mean()
                avg_radial_velocity = df['radial_velocity_mean'].iloc[1:].mean()

                results[sim] = {
                    'axial_displacement': axial_displacement,
                    'radial_displacement': radial_displacement,
                    'avg_axial_velocity': avg_axial_velocity,
                    'avg_radial_velocity': avg_radial_velocity,
                    'simulation_time': df['time_ps'].max()
                }

                print(f"Simulation time: {df['time_ps'].max():.1f} ps")
                print(f"Final axial displacement: {axial_displacement:.2f} Å")
                print(f"Final radial displacement: {radial_displacement:.2f} Å")
                print(f"Average axial velocity: {avg_axial_velocity:.4f} Å/ps")
                print(f"Average radial velocity: {avg_radial_velocity:.4f} Å/ps")

                # Save data
                df.to_csv(f'displacement_{sim}.csv', index=False)
                print(f"Data saved to displacement_{sim}.csv")
            else:
                print(f"No data found for {sim}")
        else:
            print(f"File {dump_file} not found")

    # Summary
    print(f"\n=== FINAL RESULTS ===")
    for sim in ['r10_1', 'r10_2', 'r10_3', 'r10_5']:
        if sim in results:
            data = results[sim]
            print(f"{sim}: Axial displacement = {data['axial_displacement']:+.2f} Å")
            print(f"{sim}: Radial displacement = {data['radial_displacement']:+.2f} Å")
            print(f"{sim}: Average axial velocity = {data['avg_axial_velocity']:+.4f} Å/ps")
            print(f"{sim}: Average radial velocity = {data['avg_radial_velocity']:+.4f} Å/ps")
            print()

    # Generate comparison plots if any simulations analyzed
    if dataframes:
        print(f"\nGenerating comparison plots...")
        plot_displacement_comparison(dataframes)
        print("Plots saved to r10_all_analysis.png")

if __name__ == "__main__":
    main()