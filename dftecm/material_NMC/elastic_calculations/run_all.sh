#!/bin/bash
# Master script to submit all calculations
# First run relaxations, then elastic calculations

echo "Submitting relaxation calculations..."
for dir in */Li_*/relax; do
    if [ -d "$dir" ]; then
        echo "Submitting $dir"
        cd "$dir"
        # Uncomment the line below for your scheduler:
        # sbatch submit.sh
        cd - > /dev/null
    fi
done

echo ""
echo "After relaxations complete, copy CONTCAR to elastic/POSCAR:"
echo "for dir in */Li_*/relax; do"
echo '    elastic_dir="${dir/relax/elastic}"'
echo '    cp "$dir/CONTCAR" "$elastic_dir/POSCAR"'
echo "done"
echo ""
echo "Then submit elastic calculations..."
