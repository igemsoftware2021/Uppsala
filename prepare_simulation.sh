protprep="/home/patrhenn/Downloads/qligfep_env/qligfep/protprep.py"
qresfep="/home/patrhenn/Downloads/qligfep_env/qligfep/QresFEP.py"

position=$1
from_aa=$2
to_aa=$3
mutation="${from_aa}${position}${to_aa}"
echo "preparing mutation: " $mutation

function echo_cmd ( ) {
	nextcmd=$1
	echo $nextcmd
	$nextcmd > /dev/null
}
function echo_cmd_noredirect ( ) {
	nextcmd=$1
	echo $nextcmd
	$nextcmd
}

# run the column shift script on the input file for protprep compatibility
echo_cmd "python3 shift_3rd_column_pymol_pdb.py 3-complex-ppi/WT/FGF+R2_Q.pdb 3-complex-ppi/WT/$mutation/FGF+R2_Q_aligned.pdb"
echo_cmd_noredirect "./setup_simulation.sh 3-complex-ppi/WT/$mutation FGF+R2_Q_aligned.pdb $position ${from_aa}${position}A"

# align the third column of the mutated file after pymol
echo_cmd "python3 align_3rd_column_pymol_pdb.py 3-complex-ppi/MUT/$mutation/FGF+R2_Q_${mutation}.pdb 3-complex-ppi/MUT/$mutation/FGF+R2_Q_${mutation}_aligned.pdb"
echo_cmd_noredirect "./setup_simulation.sh 3-complex-ppi/MUT/$mutation FGF+R2_Q_${mutation}_aligned.pdb  $position ${to_aa}${position}A"

# run the column shift script on the input file for protprep compatibility
echo_cmd "python3 shift_3rd_column_pymol_pdb.py 2-single-protein/WT/FGF2_Q.pdb 2-single-protein/WT/${mutation}/FGF2_Q_aligned.pdb"
echo_cmd_noredirect "./setup_simulation.sh 2-single-protein/WT/$mutation FGF2_Q_aligned.pdb $position ${from_aa}${position}A"

# align the third column of the mutated file after pymol
echo_cmd "python3 align_3rd_column_pymol_pdb.py 2-single-protein/MUT/${mutation}/FGF2_Q_${mutation}.pdb 2-single-protein/MUT/${mutation}/FGF2_Q_${mutation}_aligned.pdb"
echo_cmd_noredirect "./setup_simulation.sh 2-single-protein/MUT/$mutation FGF2_Q_${mutation}_aligned.pdb $position ${to_aa}${position}A"
