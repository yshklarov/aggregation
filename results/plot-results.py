#!/usr/bin/env python3
# Generated by Claude AI, 2025-01-02.

import re
import matplotlib.pyplot as plt
import numpy as np

def parse_results(filename):
    sizes = []
    means = []
    stddevs = []
    densities = []
    
    with open(filename, 'r') as f:
        for line in f:
            size_match = re.search(r'size (\d+)×\d+', line)
            density_match = re.search(r'density ([\d.]+)', line)
            size = int(size_match.group(1))
            density = float(density_match.group(1))
            
            # Extract mean and stddev
            stats_match = re.search(r'mean±stddev\): \d+, \d+, ([\d.]+) ± ([\d.]+)', line)
            mean_val = float(stats_match.group(1))
            stddev_val = float(stats_match.group(2))
            
            sizes.append(size)
            means.append(mean_val)
            stddevs.append(stddev_val)
            densities.append(density)
    
    return sizes, means, stddevs, densities

def create_plot(sizes, means, stddevs, densities, errorbars=True):
    plt.figure(figsize=(10, 8))
    
    # Get unique densities and assign colors
    unique_densities = sorted(set(densities))
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_densities)))
    
    # Plot points for each density value
    for density, color in zip(unique_densities, colors):
        mask = [t == density for t in densities]
        if errorbars:
            plt.errorbar(
                [sizes[i] for i in range(len(sizes)) if mask[i]],
                [means[i] for i in range(len(means)) if mask[i]],
                yerr=[stddevs[i] for i in range(len(stddevs)) if mask[i]],
                fmt='o-', capsize=3, capthick=0.5, elinewidth=0.5, markersize=4,
                color=color, label=f'Density = {density:.2f}'
            )
        else:
            plt.plot(
                [sizes[i] for i in range(len(sizes)) if mask[i]],
                [means[i] for i in range(len(means)) if mask[i]],
                'o-',
                markersize=4, color=color, label=f'Density = {density:.2f}'
            )
    
    # Customize the plot
    plt.xlabel('N (Grid is N×N)', fontsize=12)
    plt.xticks(sizes, fontsize=7)
    plt.ylabel('Number of Steps', fontsize=12)
    plt.title('Mean steps to large cluster with 90% threshold', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('results.png', dpi=72, bbox_inches='tight')
    plt.show()

def main():
    # Parse the data
    sizes, means, stddevs, densities = parse_results('results.txt')
    
    # Create and save the plot
    create_plot(sizes, means, stddevs, densities, False)

if __name__ == "__main__":
    main()
