import sys
assert sys.version_info[0]>2,"python3 required. (try running with 'python3 myscript.py' instead of 'python myscript.py')"
import os
import subprocess
import shutil

program_base_path="/home/patrhenn/Downloads/qligfep_env"

protprep_path=os.path.join(program_base_path,"qligfep/protprep.py")
qresfep_path=os.path.join(program_base_path,"qligfep/QresFEP.py")
qprep_path=os.path.join(program_base_path,"Q6/bin/qprep")
qligfep_analysis_path=os.path.join(program_base_path,"qligfep/analyze_FEP.py")

#script location on uppmax
qligfep_analysis_path="/proj/uppmax2021-2-11/qligfep_env/qligfep/analyze_FEP.py"

def is_float(num_as_string):
    try:
        f=float(num_as_string)
        return True
    except:
        return False

# return <does directory already exist?>
def create_dir_if_not_exists(dir):
    if not os.path.exists(dir):
        print(f"folder '{dir}' does not exist. creating")
        os.makedirs(dir)
        return False
    return True

class Wrapper:
    def __init__(self,directory,complex_pdb,single_pdb,mutations):
        self.directory=directory
        self.complex_pdb=complex_pdb
        self.single_pdb=single_pdb
        self.mutations=mutations

    def __str__(self):
        return f"dir:'{self.directory}',complex:'{self.complex_pdb}',single:'{self.single_pdb}',mutations:\"{self.mutations}\"\n"

    def build(self):
        #download and build the q stuff, if not already present
        # git clone q
        # git clone qligfep
        #build q
        #prepare qligfep/settings.py
        #prepare settings file for the q stuff as well, including runtime, number of cores etc.
        full_work_path=os.path.join(os.getcwd(),self.directory)
        create_dir_if_not_exists(full_work_path)
        #and save the current configuration into a file in the local directory
        with open(os.path.join(full_work_path,"setup.config"),"w+") as file:
            file.write(str(self))
        pass

    def prepare(self):
        #import here to not required a pymol installation on e.g. the cluster where the simulations are run, which usually is headless
        from perform_mutagenesis import prepare
        
        #actually prepare the simulations
        full_work_path=os.path.join(os.getcwd(),self.directory)
        #create one folder to contain the simulations for the complex
        create_dir_if_not_exists(os.path.join(full_work_path,"complex"))
        #and one folder for the single protein simulations
        create_dir_if_not_exists(os.path.join(full_work_path,"single"))
        #and copy the WT pdb files there
        single_pdb_file=os.path.split(self.single_pdb)[1]
        complex_pdb_file=os.path.split(self.complex_pdb)[1]

        single_reference_full_path=os.path.join(full_work_path,'single',single_pdb_file)
        complex_reference_full_path=os.path.join(full_work_path,'complex',complex_pdb_file)

        shutil.copyfile(self.single_pdb,single_reference_full_path)
        shutil.copyfile(self.complex_pdb,complex_reference_full_path)

        #then create folders and setup the simulations for each simulation that needs to run, based on the list of mutations given as argument
        for (upper_dir,reference_pdb,reference_full_path) in [("single",single_pdb_file,single_reference_full_path),("complex",complex_pdb_file,complex_reference_full_path)]:
            for mutation in self.mutations:
                mutation_from=mutation[0]
                mutation_to=mutation[-1]
                mutation_position=mutation[1:-1]

                for (simulation_name,is_wt) in [(mutation_to,False),(mutation_from,True)]:
                    simulation_name=simulation_name+mutation_position+"A"
                    mutation_folder=os.path.join(full_work_path,upper_dir,simulation_name)

                    #create folder to contain stuff related to (mutant to A) and (WT to A)
                    if not create_dir_if_not_exists(mutation_folder):
                        if not is_wt:
                            preparation_pdb_out_name=os.path.join(mutation_folder,f"{reference_pdb[:-4]}_mutated_{mutation_from}{mutation_position}{mutation_to}.pdb")
                            prepare(
                                in_file=os.path.join(os.path.split(mutation_folder)[0],reference_pdb),
                                out_file=preparation_pdb_out_name,
                                molecule_name="irrelevant_molecule_name",
                                chain_id="A",
                                mutation={
                                    "from":mutation_from,
                                    "to":mutation_to,
                                    "position":mutation_position
                                }
                            )
                        else:
                            preparation_pdb_out_name=os.path.join(mutation_folder,reference_pdb)
                            shutil.copyfile(reference_full_path,preparation_pdb_out_name)

                        # run the column shift script on the input file for protprep compatibility (shift_3rd_column_pymol_pdb.py and align_3rd_column_pymol_pdb.py)
                        fix_3rd_column(preparation_pdb_out_name,preparation_pdb_out_name[:-4]+"_aligned.pdb")

                        #basically setup_simulation.sh below:
                        #change directory to mutation folder
                        os.chdir(mutation_folder)
                        #run protprep
                        nextcmd=f"python3 {protprep_path} -p {os.path.split(preparation_pdb_out_name)[1][:-4]}_aligned.pdb -r 25 -c RESN:{mutation_position}"
                        assert subprocess.run(nextcmd.split()).returncode==0, "protprep likely not found"
                        #fix the 3rd column again
                        fix_3rd_column("protein.pdb","protein_fixed.pdb")
                        #move the fixed result into protein.pdb (exactly that name is required for qresfep)
                        shutil.move("protein_fixed.pdb","protein.pdb")
                        #run qresfep
                        assert subprocess.run(f"python3 {qresfep_path} -m {simulation_name} -f OPLS2015 -s linear -w 20 -S protein -T 298 -r 10 -C SNOWY -l 1 -mc A".split()).returncode==0, "qresfep likely not found"
                        #go into 'inputfiles' folder, where the qresfep output is
                        os.chdir(os.path.join(mutation_folder,f"FEP_{simulation_name}","inputfiles"))
                        #add TER line between chains/proteins in pdb file containing the protein-protein complex (qresfep wrongly strips that line)
                        if upper_dir=="complex":
                            protein_lengths=[]
                            #count amino acids per protein in pdb file
                            with open(preparation_pdb_out_name,"r") as prepared_pdb:
                                current_protein_length=0
                                for line in prepared_pdb.readlines():
                                    #count C alpha atoms, where each AA should have exactly one
                                    if line[:4]=="ATOM":
                                        if line[13:15]=="CA":
                                            current_protein_length+=1
                                    elif line[:3]=="TER":
                                        protein_lengths.append(current_protein_length)
                                        current_protein_length=0
                            #and add the terminator line between the two proteins, seperated by the knowledge of their length
                            with open("complex_with_ter.pdb","w+") as complex_with_terminator:
                                with open("complex.pdb","r") as complex:
                                    current_protein_length=0
                                    current_aa=""
                                    for line in complex.readlines():
                                        #again, count C alpha atoms as number of AAs
                                        if line[13:15]=="CA":
                                            current_protein_length+=1
                                        #save residue name from current line if the current line is the Ca atom of the last amino acid
                                        if current_protein_length==protein_lengths[0]:
                                            if line[13:15]=="CA":
                                                current_aa=line[17:26]
                                        #if residue name in current line does not match the saved residue name, it should be the first amino acid of the second chain/protein of the file, therefore a terminator line should be added before the current line
                                        if line[17:26]!=current_aa and current_aa!="":
                                            complex_with_terminator.write("TER\n")
                                            current_aa=""
                                        complex_with_terminator.write(line)
                            #rename file to fixed naming scheme
                            shutil.move("complex_with_ter.pdb","complex.pdb")
                            #rerun qprep after the ter line has been added, using this command: qprep < qprep.inp > qprep.out (but pythonified)
                            cat=subprocess.Popen(["cat","qprep.inp"],stdout=subprocess.PIPE)
                            cat.wait()
                            assert subprocess.run([qprep_path],stdin=cat.stdout,stdout=open("qprep.out","w+")).returncode==0, "qprep likely not found"
                            #check if output actually contains the relevant lines (e.g. if pdb file should include 2 chains/proteins, check if the output contains two lines matching ^molecule)
                            #grep "^molecule " qprep.out
                            molecules_in_file=[False,False]
                            with open("qprep.out","r") as file:
                                for line in file.readlines():
                                    mol="molecule "
                                    if line[:len(mol)]==mol:
                                        first_molecule="molecule    1:"
                                        if line[:len(first_molecule)]==first_molecule:
                                            molecules_in_file[0]=True
                                        second_molecule="molecule    2:"
                                        if line[:len(second_molecule)]==second_molecule:
                                            molecules_in_file[1]=True
                            if False in molecules_in_file:
                                print("could not add terminator to pdb file to seperate the molecules inside the complex. please contact a dev to fix this.")
                                sys.exit(1)
                        #add the SBATCH -M snowy line (specify computing cluster)
                        #and replace srun with $jobruntime inside the runSNOWY.sh file (change mpi runtime)
                        with open("runSNOWY.sh","r") as cluster_start_job:
                            with open("runSNOWY.sh.replacement","w+") as cluster_fixed:
                                for line_number,line in enumerate(cluster_start_job.readlines()):
                                    #add this flag in second line, which only every contains '#\n'
                                    if line_number==1:
                                        cluster_fixed.write("#SBATCH -M snowy\n")
                                    else:
                                        #and replace srun with an environment variable in all other lines
                                        cluster_fixed.write(line.replace("srun","$job_runtime"))
                        #rename file to fixed naming scheme
                        shutil.move("runSNOWY.sh.replacement","runSNOWY.sh")

    def run(self):
        #somehow run the code on a supercomputer cluster?
        #probably just prepare everything here, where a schrodinger license is available, and then instruct the user to copy everything to the supercomputer cluster and run the 'run' step there
        full_work_path=os.path.join(os.getcwd(),self.directory)
        for upper_dir in ["single","complex"]:
            for mutation in self.mutations:
                mutation_from=mutation[0]
                mutation_to=mutation[-1]
                mutation_position=mutation[1:-1]

                for (simulation_name,is_wt) in [(mutation_to,False),(mutation_from,True)]:
                    simulation_name=simulation_name+mutation_position+"A"
                    mutation_folder=os.path.join(full_work_path,upper_dir,simulation_name)

                    #create folder to contain stuff related to (mutant to A) and (WT to A)
                    os.chdir(os.path.join(mutation_folder,f"FEP_{simulation_name}"))
                    assert subprocess.run("bash FEP_submit.sh".split()).returncode==0, f"could not submit simulation {upper_dir}/{simulation_name}"

    def analyze(self):
        #just import numpy here for version checking (numpy.nanmean is used inside the external qligfep analysis script, which was introduced in numpy 1.8)
        import numpy
        assert int(numpy.__version__.split(".")[1])>=15,"numpy version >=1.15.0 required for analysis"
        #create a folder to contain the analysis results within each simulation set folder
        full_work_path=os.path.join(os.getcwd(),self.directory)
        
        with open(os.path.join(full_work_path,"results.out"),"w+") as result_out_file:
            for mutation in self.mutations:
                mutation_from=mutation[0]
                mutation_to=mutation[-1]
                mutation_position=mutation[1:-1]

                delta_g={
                    "single_wt":0.0,
                    "single_mut":0.0,
                    "complex_wt":0.0,
                    "complex_mut":0.0,
                }   
                sem=0.0

                for upper_dir in ["single","complex"]:
                    for (simulation_name,is_wt) in [(mutation_to,False),(mutation_from,True)]:
                        simulation_name=simulation_name+mutation_position+"A"

                        mutation_folder=os.path.join(full_work_path,upper_dir,simulation_name)

                        os.chdir(mutation_folder)

                        if True or not os.path.exists("analyze.out"):
                            with open("analyze.out","w+") as analysis_out_file:
                                #and run the analysis scipt in there
                                assert subprocess.run(f"python3 {qligfep_analysis_path} -F FEP_{simulation_name} -C SNOWY".split(),stdout=analysis_out_file).returncode==0
                        with open("analyze.out","r") as analyis_out_file:
                            #line:=FEP_mutation deltag deltag_sem deltag_forward deltag_forward_sem deltag_reverse deltag_reverse_sem deltag_overlap_sampling deltag_overlap_sampling_sem deltag_bennet_acceptance_ratio deltag_bennet_acceptance_ratio_sem
                            #e.g. "FEP_E58A 217.05  0.82 216.79  0.73 -217.32  0.92 217.03  0.82 217.04  0.82"
                            for line in analyis_out_file.readlines():
                                if line[:4]=="FEP_":
                                    cells=[]
                                    for cell in line.split():
                                        #cells are seperated by at least one space
                                        if len(cell)>0 and cell[:4]!="FEP_":
                                            assert is_float(cell),cell
                                            cells.append(float(cell))

                            sem+=cells[1]
                            if upper_dir=="complex":
                                if is_wt:
                                    delta_g["complex_wt"]=cells[0]
                                else:
                                    delta_g["complex_mut"]=cells[0]
                            else:
                                if is_wt:
                                    delta_g["single_wt"]=cells[0]
                                else:
                                    delta_g["single_mut"]=cells[0]

                assert delta_g["single_wt"]!=0.0
                assert delta_g["single_mut"]!=0.0
                assert delta_g["complex_wt"]!=0.0
                assert delta_g["complex_mut"]!=0.0
                                        
                #then combine the output from the output files into a small set of relevant numbers, and print these here
                result=(delta_g["complex_wt"]-delta_g["single_wt"])-(delta_g["complex_mut"]-delta_g["single_mut"])
                result_string='delta delta g for {}: {:.2f} kcal/mol (cw:{:.2f}-sw:{:.2f})-(cm:{:.2f}-sm:{:.2f}), sum of sems: {:.2f} kcal/mol'.format(
                    mutation,
                    result,
                    delta_g["complex_wt"],delta_g["single_wt"],delta_g["complex_mut"],delta_g["single_mut"],
                    sem
                )
                print(result_string)
                
                result_out_file.write(result_string+"\n")

def fix_3rd_column(in_filename,out_filename=None):
    if not out_filename:
        out_filename=in_filename[:-4]+"_column_3_fixed.pdb"

    with open(in_filename) as infile:
        with open(out_filename,"w+") as outfile:
            for line_number,line in enumerate(infile.readlines()):
                    line_name=line[:4]
                    if line_name.strip() not in {"ATOM","END","TER"}:
                        print(f"line type is not atom/ter/end: {line} (line:{line_number+1})")
                        outfile.write(line)
                        continue

                    if len(line)>4 and line_name.strip()=="ATOM":
                        assert line[4]==" "

                        if len(line)>11:
                            atom_number=line[5:11]
                            atom_number=int(atom_number)
                        
                            assert line[11]==" "

                            atom_type=line[12:17]

                            try:
                                new_atom_type=atom_type.strip()
                                if new_atom_type[0].isnumeric():
                                    new_atom_type=new_atom_type[1:]+new_atom_type[0]
                                if new_atom_type[0].isnumeric():
                                    new_atom_type=new_atom_type[1:]+new_atom_type[0]
                            except:
                                print("something went wrong",line)
                                sys.exit(1)
                            assert not new_atom_type[0].isnumeric()

                            assert len(new_atom_type)<5

                            new_atom_type=" "+new_atom_type+" "*(4-len(new_atom_type))

                            outfile.write(line.replace(atom_type,new_atom_type))
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)
