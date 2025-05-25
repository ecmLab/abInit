##### This file is used to remove water and ions that we don't want and
##    ensure the net charge is zero by adding additional ions to balance the mos2 charge density #####
#####  By Howard Tu 2025.4.20  #####
### Notice ###
# There are 4 steps in this file:
  ## step1: Align to center and remove water molecules that we don't want---
  ## step2: set partial charge to atoms ---
  ## step4: Delete some ions so that the net charge will be zero---
## Variables used in the file
## Constant Value (in armstrong)
set  pi    3.1415926535897932
set  rr    20
set xyz [list [expr 2.5*$rr] [expr 2.5*$rr] 12.5]
package require topotools

# Align the mos2 tube to the center of the system
mol delete all
mol load pdb sol_r${rr}.pdb
pbc set $xyz
set all [atomselect top all]
set mos2 [atomselect top "resname MOL"]
set vec [measure center $mos2]
$all moveby [vecinvert $vec ]
$all moveby [vecscale 0.5 $xyz]
$all writepdb sol_r${rr}_fix.pdb

##Set partial charge to atoms
mol delete all
mol load pdb sol_r${rr}_fix.pdb
## partial charge in H2O
set H2O_O [atomselect top "name O"]
$H2O_O set charge -0.834
set H2O_H [atomselect top "name H1 H2"]
$H2O_H set charge 0.417
## partial charge in (HCrO4)-
set Cr [atomselect top "name CR"]
$Cr set charge 0.813
set O  [atomselect top "name O1 O2 O3 O4"]
$O set charge -0.703
set H  [atomselect top "name H"]
$H set charge 1.0
## partial charge in MoS2 with charge density 0.006 e/A^2
set Mo  [atomselect top "name Mo"]
$Mo set charge 0.734
set S  [atomselect top "name S"]
$S set charge -0.37425
## charge for Na and Cl
set Na  [atomselect top "name NA"]
$Na set charge 1.0
set Cl  [atomselect top "name CL"]
$Cl set charge -1.0

topo clearbonds
topo clearangles
# Add Cr–O bonds in HCrO4 molecules
set nCr [$Cr num]
set crid [$Cr get index]
# Iterate through each Cr and find bonded oxygen
for {set i 0} {$i < $nCr} {incr i} {
    set iCr  [lindex $crid $i]
    set aCr [atomselect top "index $iCr"]
    set resid  [$aCr get resid]
    # Find O and H atoms belong to this Cr
    set O1 [atomselect top "name O1 and resid $resid"]
    set O1id [lindex [$O1 get index] 0]
    set O2 [atomselect top "name O2 and resid $resid"]
    set O2id [lindex [$O2 get index] 0]
    set O3 [atomselect top "name O3 and resid $resid"]
    set O3id [lindex [$O3 get index] 0]
    set O4 [atomselect top "name O4 and resid $resid"]
    set O4id [lindex [$O4 get index] 0]
    set H [atomselect top "name H and resid $resid"]
    set hid [lindex [$H get index] 0]

    # Add Cr–O bonds
    topo addbond $iCr $O1id -bondtype 1
    topo addbond $iCr $O2id -bondtype 1
    topo addbond $iCr $O3id -bondtype 1
    topo addbond $iCr $O4id -bondtype 1
    # Add O–H bonds
    topo addbond $O1id $hid -bondtype 2

    # Add O–Cr–O angle
    topo addangle $O1id $iCr $O2id 1
    topo addangle $O1id $iCr $O3id 1
    topo addangle $O1id $iCr $O4id 1
    topo addangle $O2id $iCr $O3id 1
    topo addangle $O2id $iCr $O4id 1
    topo addangle $O3id $iCr $O4id 1
}

# Add O-H bonds in water molecules
set nO [$H2O_O num]
set oid [$H2O_O get index]
# Iterate through each oxygen and find bonded hydrogens
for {set i 0} {$i < $nO} {incr i} {
    set io  [lindex $oid $i]
    set iwo [atomselect top "index $io"]
    set resid  [$iwo get resid]
    # Find H atoms bonded to this oxygen
    set H1 [atomselect top "name H1 and resid $resid"]
    set h1id [lindex [$H1 get index] 0]
    set H2 [atomselect top "name H2 and resid $resid"]
    set h2id [lindex [$H2 get index] 0]

    # Add O–H bonds
    topo addbond $io $h1id -bondtype 3
    topo addbond $io $h2id -bondtype 3

    # Add H–O–H angle
    topo addangle $h1id $io $h2id 2
}

# Guess bond and angle types if needed
#topo guessbondtypes
#topo guessangletypes

set all [atomselect top all]
puts " TOTAL CHARGE ON THE SYSTEM: [eval vecadd [$all get charge]]"
$all writepdb sol_r${rr}.pdb

# Export LAMMPS data file with charges
package require topotools
topo writelammpsdata sol_r${rr}.data full

mol delete all
mol load pdb sol_r${rr}.pdb
