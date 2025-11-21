#!/bin/bash
# Helper script to generate POTCAR files
# Modify VASP_PP_PATH to your pseudopotential directory

VASP_PP_PATH="/Users/howardtu/Documents/modeling/localPackage/potpaw_PBE64"

# Recommended PAW potentials for NMC:
# Li_sv (1s2s2p), Ni_pv (3p4s3d), Mn_pv (3p4s3d), Co (4s3d), O (2s2p)

for dir in */Li_*/relax */Li_*/elastic; do
    if [ -d "$dir" ]; then
        echo "Generating POTCAR for $dir"
        cd "$dir"

        # Read elements from POSCAR line 6
        elements=$(sed -n '6p' POSCAR)

        # Generate POTCAR
        rm -f POTCAR
        for el in $elements; do
            case $el in
                Li) cat "$VASP_PP_PATH/Li_sv/POTCAR" >> POTCAR ;;
                Ni) cat "$VASP_PP_PATH/Ni_pv/POTCAR" >> POTCAR ;;
                Mn) cat "$VASP_PP_PATH/Mn_pv/POTCAR" >> POTCAR ;;
                Co) cat "$VASP_PP_PATH/Co/POTCAR" >> POTCAR ;;
                O)  cat "$VASP_PP_PATH/O/POTCAR" >> POTCAR ;;
            esac
        done

        cd - > /dev/null
    fi
done

echo "POTCAR generation complete!"
