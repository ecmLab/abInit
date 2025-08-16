##### This file is used to remove water and ions that we don't want and
##    ensure the net charge is zero by adding additional ions to balance the mos2 charge density #####
#####  By Howard Tu 2025.7.30  #####
### Notice ###
# There are 4 steps in this file:
  ## step1: Align to center---
  ## step2: set partial charge to atoms ---
  ## step3: Build correct topoloy
  ## step4: Delete some ions so that the net charge will be zero---
## Variables used in the file
## Constant Value (in armstrong)
set  pi    3.1415926535897932
set  a0    3.1828

set  nn    100
## Radius of Mo layer
# when nn = 10, rr = 8.7738 ;  when nn = 20, rr  = 17.5477;
# when nn = 50, rr = 43.8692;  when nn = 100, rr = 87.7385;
if {$nn == 10} {
   set  rr    8.7738
} elseif {$nn == 20} {
   set  rr    17.5477
} elseif {$nn == 50} {
   set  rr    43.8692
} elseif {$nn == 100} {
   set  rr    87.7385
}

set xyz [list [expr 7*$a0] [expr 4.0*$rr] [expr 4.0*$rr]]
package require topotools

# Align the mos2 tube to the center of the system
mol delete all
mol load pdb sol_N${nn}.pdb
pbc set $xyz
set all [atomselect top all]
set mos2 [atomselect top "resname MOL"]
set vec [measure center $mos2]
$all moveby [vecinvert $vec ]
$all moveby [vecscale 0.5 $xyz]
$all writepdb sol_N${nn}_fix.pdb

##Set partial charge to atoms
mol delete all
mol load pdb sol_N${nn}_fix.pdb
## partial charge in H2O
set H2O_O [atomselect top "name O"]
$H2O_O set charge -0.834
set H2O_H [atomselect top "name H1 H2"]
$H2O_H set charge 0.417
## partial charge in (HCrO4)-
set Cr [atomselect top "name CR"]
$Cr set charge 0.812
set O  [atomselect top "name O1 O2 O3 O4"]
$O set charge -0.703
set H  [atomselect top "name H"]
$H set charge 1.0
## partial charge in MoS2 with charge density 1/140=0.0071 e/A^2
set Mo  [atomselect top "name Mo1 Mo2"]
$Mo set charge 0.7340000
set S  [atomselect top "name S1 S2 S3 S4"]
$S set charge -0.3705714
## charge for Na and Cl
set Na  [atomselect top "name NA"]
$Na set charge 1.0
set Cl  [atomselect top "name CL"]
$Cl set charge -1.0
set all [atomselect top all]

##Build topo relations
topo clearbonds
topo clearangles
set Mo1  [atomselect top "name Mo1"]
set nMo1 [$Mo1 num]
set moid [$Mo1 get index]
# Iterate through all Mo1 and find bonded sulfur
for {set i 0} {$i < $nMo1} {incr i} {
    set iMo1  [lindex $moid $i]
    set aMo1 [atomselect top "index $iMo1"]
    set resid  [$aMo1 get resid]
 # Find the other Mo atom having the same resid
    set aMo2 [atomselect top "name Mo2 and resid $resid"]
    set iMo2 [lindex [$aMo2 get index] 0]
 # Find S atoms having the same resid
    set S1 [atomselect top "name S1 and resid $resid"]
    set iS1 [lindex [$S1 get index] 0]
    set S2 [atomselect top "name S2 and resid $resid"]
    set iS2 [lindex [$S2 get index] 0]
    set S3 [atomselect top "name S3 and resid $resid"]
    set iS3 [lindex [$S3 get index] 0]
    set S4 [atomselect top "name S4 and resid $resid"]
    set iS4 [lindex [$S4 get index] 0]

## Add Mo1-S bonds, with bond length=2.4127A
  # Add Mo1-S bonds in the same resid
    topo addbond $iMo1 $iS1 -bondtype 1
    topo addbond $iMo1 $iS2 -bondtype 1
  # Add Mo1-S bonds with its neighing resid along x-1 direction id(x-1) = id - 6*nn
    set iS1_xm1 [expr $iS1-$nn*6]
    set iS2_xm1 [expr $iS2-$nn*6]
    if {$resid <= $nn} {
       set iS1_xm1 [expr $iS1+$nn*36]
       set iS2_xm1 [expr $iS2+$nn*36]
    }
    topo addbond $iMo1 $iS1_xm1 -bondtype 1
    topo addbond $iMo1 $iS2_xm1 -bondtype 1
  # Add Mo1-S bonds with its neighing resid along y-1 direction id(y-1) = id - 6
    set iS3_ym1 [expr $iS3-6]
    set iS4_ym1 [expr $iS4-6]
    if {[ expr $resid % $nn ] == 1} {
       set iS3_ym1 [expr $iS3+($nn-1)*6]
       set iS4_ym1 [expr $iS4+($nn-1)*6]
    }
    topo addbond $iMo1 $iS3_ym1 -bondtype 1
    topo addbond $iMo1 $iS4_ym1 -bondtype 1

## Add Mo2-S bonds, with bond length=2.4127A
  # Add Mo2-S bonds in the same resid
    topo addbond $iMo2 $iS1 -bondtype 1
    topo addbond $iMo2 $iS2 -bondtype 1
    topo addbond $iMo2 $iS3 -bondtype 1
    topo addbond $iMo2 $iS4 -bondtype 1
  # Add Mo2-S bonds with its neighing resid along x+1 direction id(x+1) = id + 6*nn
    set iS3_xp1 [expr $iS3+$nn*6]
    set iS4_xp1 [expr $iS4+$nn*6]
    if {$resid > $nn*6} {
       set iS3_xp1 [expr $iS3-$nn*36]
       set iS4_xp1 [expr $iS4-$nn*36]
    }
    topo addbond $iMo2 $iS3_xp1 -bondtype 1
    topo addbond $iMo2 $iS4_xp1 -bondtype 1

## Add Mo-S1-Mo angle, and Mo-S2-Mo angle
  # Add Mo1 along x+1 direction id(x+1) = id + 6*nn
    set iMo1_xp1 [expr $iMo1+$nn*6]
    if {$resid > $nn*6} {
       set iMo1_xp1 [expr $iMo1-$nn*36]
    }
    topo addangle $iMo1 $iS1 $iMo2 1
    topo addangle $iMo1 $iS1 $iMo1_xp1 1
    topo addangle $iMo2 $iS1 $iMo1_xp1 1
    topo addangle $iMo1 $iS2 $iMo2 1
    topo addangle $iMo1 $iS2 $iMo1_xp1 1
    topo addangle $iMo2 $iS2 $iMo1_xp1 1

## Add Mo-S3-Mo angle, and Mo-S4-Mo angle
  # Add Mo1 along y+1 direction id(y+1) = id - 6
    set iMo1_yp1 [expr $iMo1+6]
    if {[ expr $resid % $nn ] == 0} {
       set iMo1_yp1 [expr $iMo1-($nn-1)*6]
    }
  # Add Mo2 along x-1 direction id(x-1) = id - 6*nn
    set iMo2_xm1 [expr $iMo2-$nn*6]
    if {$resid <= $nn} {
       set iMo2_xm1 [expr $iMo2_xm1+$nn*36]
    }

    topo addangle $iMo2 $iS3 $iMo2_xm1 1
    topo addangle $iMo2 $iS3 $iMo1_yp1 1
    topo addangle $iMo1_yp1 $iS3 $iMo2_xm1 1
    topo addangle $iMo2 $iS4 $iMo2_xm1 1
    topo addangle $iMo2 $iS4 $iMo1_yp1 1
    topo addangle $iMo1_yp1 $iS4 $iMo2_xm1 1

## Add Mo-S plane impropers

}

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
    topo addbond $iCr $O1id -bondtype 2
    topo addbond $iCr $O2id -bondtype 2
    topo addbond $iCr $O3id -bondtype 2
    topo addbond $iCr $O4id -bondtype 2
    # Add O–H bonds
    topo addbond $O1id $hid -bondtype 3

    # Add O–Cr–O angle
    topo addangle $O1id $iCr $O2id 2
    topo addangle $O1id $iCr $O3id 2
    topo addangle $O1id $iCr $O4id 2
    topo addangle $O2id $iCr $O3id 2
    topo addangle $O2id $iCr $O4id 2
    topo addangle $O3id $iCr $O4id 2
}

# Add O-H bonds in water molecules
set nO [$H2O_O num]
set oid [$H2O_O get index]
# Iterate through each oxygen and find bonded hydrogens
for {set i 0} {$i < $nO} {incr i} {
    set io  [lindex $oid $i]
    set iwo [atomselect top "index $io"]
    set resid  [$iwo get resid]
    set chain  [$iwo get chain]
    # Find H atoms bonded to this oxygen
    set H1 [atomselect top "name H1 and resid $resid and chain $chain"]
    set h1id [lindex [$H1 get index] 0]
    set H2 [atomselect top "name H2 and resid $resid and chain $chain"]
    set h2id [lindex [$H2 get index] 0]

    # Add O–H bonds
    topo addbond $io $h1id -bondtype 4
    topo addbond $io $h2id -bondtype 4

    # Add H–O–H angle
    topo addangle $h1id $io $h2id 3
}

# Guess bond and angle types if needed
#topo guessbondtypes
#topo guessangletypes

set all [atomselect top all]
puts " TOTAL CHARGE ON THE SYSTEM: [eval vecadd [$all get charge]]"
$all writepdb sol_N${nn}.pdb

# Export LAMMPS data file with charges
package require topotools
topo writelammpsdata sol_N${nn}.data full

mol delete all
mol load pdb sol_N${nn}.pdb
