# Critical Issues Found in Na₃SbS₄·9H₂O DFT Calculations

## Date: 2026-01-24

## Issue 1: D3 Dispersion Correction NOT Enabled

**Problem**: All INCAR files have VDW/IVDW commented out:
```
#VDW    = 12   # D3(BJ), if compiled as IVDW
#IVDW   = 12
```

**Affected calculations**:
- Na₃SbS₄ (anhydrous)
- Na₃SbS₄·8H₂O
- Na₃SbS₄·9H₂O

**Impact**:
- Calculations were run with **PBE only**, NOT PBE+D3
- Van der Waals interactions (critical for water-framework binding) are NOT captured
- Hydration energies are likely severely underestimated
- This explains why hydration energies seem too small!

**Expected effect**:
- Without D3: Water binding is too weak (what we calculated)
- With D3: Water binding should be ~2-4 eV stronger per H₂O molecule

---

## Issue 2: Na₃SbS₄·9H₂O Relaxation Did NOT Converge

**Problem**: The relax calculation stopped at NSW=100 (maximum limit)

**Evidence**:
```
OSZICAR line 100: F= -.52039720E+03 E0= -.52039042E+03  d E =-.733191E-01
```
The "d E" (energy change) is still -0.073 eV, indicating optimization was still progressing.

**Impact**:
- Final structure is NOT at the true energy minimum
- Forces on atoms have not fully converged
- True equilibrium energy could be significantly different

---

## Issue 3: Large Energy Discrepancy Between Relax and Static

**Findings**:
- **Relax final energy**: -520.397 eV (for 4 formula units)
- **Static energy**: -529.236 eV (for 4 formula units)
- **Difference**: 8.84 eV (static is LOWER!)

**Analysis**:
- Both calculations used identical POSCAR (verified by MD5 checksum)
- Static should give nearly identical energy to relaxed structure
- 8.84 eV difference is physically unreasonable
- Suggests possible calculation error or corrupted data

**Possible explanations**:
1. VASP version difference between relax and static runs
2. Different POTCAR files used
3. Electronic convergence issue
4. VASP bug or numerical instability
5. The static calculation may have been restarted from wrong WAVECAR

---

## Issue 4: Incorrect Hydration Energy Calculation

**Current formula used**:
```
ΔE_hyd = E(Na₃SbS₄·nH₂O) - E(Na₃SbS₄) - n×E(H₂O, gas)
```

**Results with E(H₂O) = -14.22 eV**:
- 8H₂O: ΔE = -5.08 eV (-0.63 eV per H₂O) - too weak!
- 9H₂O: ΔE = +25.34 eV (POSITIVE - unphysical!)

**Root cause**: Combination of:
1. No D3 dispersion (underbinding water)
2. Using wrong static energy for 9H₂O
3. 9H₂O structure not fully relaxed

---

## Corrected Analysis Using RELAX Energies

If we use the **relax** energies instead of static:

**Per formula unit energies**:
- Na₃SbS₄: -29.673 eV (from static - appears okay)
- Na₃SbS₄·8H₂O: -148.508 eV ÷ 4 = -37.127 eV (from static)
- Na₃SbS₄·9H₂O: -520.397 eV ÷ 4 = **-130.099 eV** (from RELAX, not static!)

**Corrected hydration energies** (still without D3):
- ΔE(9H₂O) = -130.099 - (-29.673) - 9×(-14.22) = **+27.55 eV** (still positive!)

Even using the relax energy, 9H₂O is thermodynamically unstable.

---

## Recommended Actions

### Immediate (Critical):

1. **Re-run ALL calculations with D3 enabled**:
   - Edit all INCARs to uncomment: `IVDW = 12`
   - This will add ~2-4 eV binding per water molecule

2. **Extend 9H₂O relaxation**:
   - Increase NSW to 200 or 300
   - Copy CONTCAR → POSCAR and continue relaxation
   - Check if structure is physically reasonable

3. **Verify H₂O molecule energy**:
   - Run isolated H₂O calculation WITH D3
   - Should get E(H₂O) ≈ -14.22 eV (PBE) or ~-14.8 eV (PBE+D3)

4. **Re-run static calculations**:
   - Use final converged CONTCAR from relaxation
   - Enable D3 in INCAR
   - Verify energies are consistent

### Medium-term:

5. **Validate 9H₂O crystal structure**:
   - Check against experimental structure from Tian et al.
   - Visualize to ensure water molecules are properly bound
   - Consider alternative structures if current one is unphysical

6. **Benchmark against experimental hydration enthalpy**:
   - Compare DFT hydration energies with calorimetry data (if available)
   - Typical hydrogen bonding: 0.2-0.5 eV per H₂O

---

## Current Status of Thermodynamic Report

The report currently uses:
- **E(H₂O) = -14.22 eV** (assumed, needs verification)
- **Energies from static calculations** (which may be incorrect for 9H₂O)
- **No D3 correction** (missing critical physics)

**Conclusion**: The thermodynamic analysis needs to be completely re-done once:
1. D3 is properly enabled
2. All structures are fully converged
3. H₂O reference energy is confirmed
4. 9H₂O structure is validated

---

## Expected Results After Corrections

With D3 enabled, we expect:
- 8H₂O: ΔE ≈ -20 to -30 eV total (-2.5 to -3.8 eV per H₂O)
- 9H₂O: ΔE should be negative if structure is valid
- Much stronger water binding
- Results more consistent with experimental reversibility

---

**Investigator**: Claude (AI assistant)
**User**: Howard Tu
**Date**: January 24, 2026
