#!/usr/bin/env python3
"""Combine two dehydration-related figures into a single row (1x2 layout)."""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

# Set up the figure with 1 row x 2 columns
fig = plt.figure(figsize=(14, 5.5))
gs = GridSpec(1, 2, figure=fig, hspace=0.10, wspace=0.20,
              left=0.05, right=0.95, top=0.95, bottom=0.05)

# Load images (only 2 figures: dehydration thermodynamics and phase diagram)
# Figure 1 (hydration vs RH) removed - not relevant for solution-phase synthesis
img1 = mpimg.imread('../../docs/figure2_dehydration_vs_T.png')
img2 = mpimg.imread('../../docs/figure3_phase_diagram.png')

# Create 1x2 subplots
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])

# Display images
ax1.imshow(img1)
ax1.axis('off')
ax1.text(0.02, 0.98, '(a)', transform=ax1.transAxes,
         fontsize=16, fontweight='bold', va='top')

ax2.imshow(img2)
ax2.axis('off')
ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes,
         fontsize=16, fontweight='bold', va='top')

# Save combined figure
plt.savefig('../../docs/figure_combined.png', dpi=300, bbox_inches='tight')
plt.savefig('../../docs/figure_combined.pdf', bbox_inches='tight')
print("Combined figure saved: figure_combined.png/pdf (1x2 layout)")
print("  Panel (a): Dehydration free energy vs temperature")
print("  Panel (b): Temperature-pressure phase diagram")
plt.close()
