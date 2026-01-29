#!/usr/bin/env python3
"""Combine crystal structure images into a single figure."""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

# Set up the figure with 1 row x 3 columns
fig = plt.figure(figsize=(15, 4))
gs = GridSpec(1, 3, figure=fig, hspace=0.05, wspace=0.10,
              left=0.05, right=0.95, top=0.90, bottom=0.15)

# Load images
img1 = mpimg.imread('../../docs/NAS.png')
img2 = mpimg.imread('../../docs/H2O.png')
img3 = mpimg.imread('../../docs/NAS8H2O.png')

# Create 1x3 subplots
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

# Display images
ax1.imshow(img1)
ax1.axis('off')
ax1.set_title('Na$_3$SbS$_4$', fontsize=14, fontweight='bold', pad=10)
ax1.text(0.02, 0.98, '(a)', transform=ax1.transAxes,
         fontsize=16, fontweight='bold', va='top', color='black',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax2.imshow(img2)
ax2.axis('off')
ax2.set_title('H$_2$O', fontsize=14, fontweight='bold', pad=10)
ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes,
         fontsize=16, fontweight='bold', va='top', color='black',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax3.imshow(img3)
ax3.axis('off')
ax3.set_title('Na$_3$SbS$_4$·8H$_2$O', fontsize=14, fontweight='bold', pad=10)
ax3.text(0.02, 0.98, '(c)', transform=ax3.transAxes,
         fontsize=16, fontweight='bold', va='top', color='black',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Save combined figure
plt.savefig('../../docs/figure_structures.png', dpi=300, bbox_inches='tight')
plt.savefig('../../docs/figure_structures.pdf', bbox_inches='tight')
print("Crystal structure figure saved: figure_structures.png/pdf")
plt.close()
