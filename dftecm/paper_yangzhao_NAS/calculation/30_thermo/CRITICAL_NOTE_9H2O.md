# CRITICAL ISSUE: Na₃SbS₄·9H₂O Calculation Results

## Problem Identified

The DFT energy for Na₃SbS₄·9H₂O appears anomalously high, resulting in a **positive hydration energy**.

## Calculated Energies (with E(H₂O) = -14.22 eV):

### 8H₂O Hydrate (NORMAL):
- E(Na₃SbS₄·8H₂O) = -148.508 eV
- E(Na₃SbS₄) = -29.673 eV
- 8×E(H₂O) = -113.76 eV
- **ΔE_hyd = -148.508 - (-29.673) - (-113.76) = -5.08 eV** ✓
- Per H₂O: -0.63 eV (favorable binding)

### 9H₂O Hydrate (ANOMALOUS):
- E(Na₃SbS₄·9H₂O) = -132.309 eV
- E(Na₃SbS₄) = -29.673 eV
- 9×E(H₂O) = -127.98 eV
- **ΔE_hyd = -132.309 - (-29.673) - (-127.98) = +25.34 eV** ❌
- Per H₂O: +2.82 eV (unfavorable!)

## Physical Interpretation

The positive hydration energy means:
```
Na₃SbS₄(s) + 9H₂O(g) → Na₃SbS₄·9H₂O(s)    ΔE = +25.34 eV (UNFAVORABLE)
```

This indicates the 9H₂O structure is **thermodynamically unstable** - it's energetically more favorable to have separated Na₃SbS₄ + 9 H₂O molecules than the "hydrate" structure.

## Possible Explanations

1. **Incorrect Structure**: The crystal structure used for 9H₂O DFT may not match the experimentally observed phase

2. **Convergence Issues**: The calculation may not have fully converged to the true minimum energy structure

3. **High-Energy Metastable State**: The structure may represent a kinetically accessible but thermodynamically unfavorable configuration

4. **Calculation Error**: OUTCAR energy extraction or unit cell formula unit counting may be incorrect

5. **Structure Collapse**: The 9H₂O structure may have collapsed or distorted during relaxation to a non-physical configuration

## Recommended Actions

1. **Verify OUTCAR energy**:
   - Check: `grep "free  energy   TOTEN" OUTCAR | tail -1`
   - Confirm: -529.236 eV for 4 formula units?

2. **Check structure convergence**:
   - Review OSZICAR for proper electronic and ionic convergence
   - Verify forces in final OUTCAR (should be < 0.02 eV/Å)

3. **Visualize final structure**:
   - Check if water molecules are properly bound or escaped from lattice
   - Verify crystal structure hasn't collapsed

4. **Re-run calculation**:
   - Start from experimental crystal structure (Tian et al., Joule 2019)
   - Use tighter convergence criteria
   - Check intermediate geometries

5. **Alternative approach**:
   - Focus thermodynamic analysis on 8H₂O (which shows reasonable behavior)
   - Note that experimental 9H₂O formation via solution precipitation may not be captured by gas-phase thermodynamics

## Current Status in Report

The report now:
- Presents corrected hydration energies using E(H₂O) = -14.22 eV
- Notes the anomalous 9H₂O result
- Focuses analysis primarily on well-behaved 8H₂O system
- Acknowledges that solution-phase precipitation may access different pathways

## Impact on Conclusions

- **Hydration to 8H₂O**: Thermodynamically favorable (-5.08 eV)
- **Hydration to 9H₂O**: DFT result is unphysical, requires verification
- **Recycling viability**: Still supported based on 8H₂O thermodynamics
- **Experimental 9H₂O formation**: May be kinetically controlled or involve solvent effects not captured in gas-phase model

---

**Date**: 2026-01-24
**Identified by**: Analysis of corrected hydration energies with proper H₂O reference
