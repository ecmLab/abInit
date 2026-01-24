#!/usr/bin/env python3
"""Combine all figures into a single 2x2 multi-panel figure."""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

# Set up the figure with 2x2 subplots
fig = plt.figure(figsize=(14, 12))
gs = GridSpec(2, 2, figure=fig, hspace=0.15, wspace=0.15,
              left=0.05, right=0.95, top=0.95, bottom=0.05)

# Load images
img1 = mpimg.imread('../00_docs/figure1_hydration_vs_RH.png')
img2 = mpimg.imread('../00_docs/figure2_dehydration_vs_T.png')
img3 = mpimg.imread('../00_docs/figure3_phase_diagram.png')
img4 = mpimg.imread('../00_docs/figure4_energy_levels.png')

# Create 2x2 subplots
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

# Display images
ax1.imshow(img1)
ax1.axis('off')
ax1.text(0.02, 0.98, '(a)', transform=ax1.transAxes,
         fontsize=14, fontweight='bold', va='top')

ax2.imshow(img2)
ax2.axis('off')
ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes,
         fontsize=14, fontweight='bold', va='top')

ax3.imshow(img3)
ax3.axis('off')
ax3.text(0.02, 0.98, '(c)', transform=ax3.transAxes,
         fontsize=14, fontweight='bold', va='top')

ax4.imshow(img4)
ax4.axis('off')
ax4.text(0.02, 0.98, '(d)', transform=ax4.transAxes,
         fontsize=14, fontweight='bold', va='top')

# Save combined figure
plt.savefig('../00_docs/figure_combined.png', dpi=300, bbox_inches='tight')
plt.savefig('../00_docs/figure_combined.pdf', bbox_inches='tight')
print("Combined figure saved: figure_combined.png/pdf (2x2 layout)")
plt.close()
