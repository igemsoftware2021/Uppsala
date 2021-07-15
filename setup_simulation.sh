#!/usr/bin/env bash

protprep="/home/patrhenn/Downloads/qligfep_env/qligfep/protprep.py"
qresfep="/home/patrhenn/Downloads/qligfep_env/qligfep/QresFEP.py"
qprep="/home/patrhenn/Downloads/qligfep_env/Q6/bin/qprep"

current_wd=$(pwd)

folder=$1 # where the protein file below is located, which is also where all of the setup files will be placed
protein_file=$2 # in pdb format
position=$3 # e.g. 88
mutation=$4 # eg. V88A

echo "residue: ${position}"
echo "mutation: ${mutation}"

cd $folder
echo "current wd:" $(pwd)

# run protprep
nextcmd="python $protprep -p $protein_file -r 25 -c RESN:$position"
echo $nextcmd
$($nextcmd)

# shift third column, save protein.pdb output from protprep, and move aligned file into its place
nextcmd="python3 $current_wd/shift_3rd_column_pymol_pdb.py protein.pdb protein_aligned.pdb"
echo $nextcmd
$nextcmd > /dev/null

# overwrite protprep protein output with aligned protein (step is really fast, so no backup necessary)
#mv protein.pdb protein.pdb.nonaligned
mv protein_aligned.pdb protein.pdb

# run qresfep
nextcmd="python $qresfep -m $mutation -f OPLS2015 -s linear -w 20 -S protein -T 298 -r 10 -C SNOWY -l 1 -mc A"
echo $nextcmd
$nextcmd

cd FEP_$mutation/inputfiles
python ${current_wd}/add_ter.py complex.pdb
mv complex_with_term.pdb complex.pdb
$qprep < qprep.inp > qprep.out
cat qprep.out | grep "^molecule "


cd $current_wd

# add cluster option to sbatch script
sed -i '2s/#.*/#SBATCH -M snowy/' $folder/FEP_$mutation/inputfiles/runSNOWY.sh
# replace srun with mpirun runtime, and specify number of cores to use
sed -i 's/srun/$job_runtime/g' $folder/FEP_$mutation/inputfiles/runSNOWY.sh
