#!/bin/env python3

# https://pymolwiki.org/index.php/Launching_From_a_Script # how to use pymol as cli
# https://pymolwiki.org/index.php/Iterate # how to get information from selected atoms
# https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/wizard/mutagenesis.py # pymol mutagenesis internals, including bump_scores (strain value for rotamers), and an example on pymol.storage for pymol interal->external communication
# http://www.mayachemtools.org/docs/scripts/html/code/PyMOLMutateAminoAcids.html # the inspiration for this script, though not much (none?) of that code is actually included in this file
# http://www.mayachemtools.org/docs/modules/html/PyMOLUtil.py.html#GetChains # reference for pymol chain retrieval, though not much (none?) of that code is actually included in this file
# http://www.endmemo.com/bio/codon.php # codon table because i cannot remember them myself

import sys

try:
    import pymol
    # wait for pymol to finish launching cli (-c), in quiet mode (-q), and disable loading of pymolrc and plugins (-k)
    pymol.finish_launching(['pymol', '-ckq'])
except ImportError as ErrMsg:
    sys.stderr.write(f"Failed to import PyMOL package, most likely because pymol is not installed: {ErrMsg}\n")
    sys.exit(1)

aa_table_one_to_three={
    "A":"ALA",
    "C":"CYS",
    "D":"ASP",
    "E":"GLU",
    "F":"PHE",
    "G":"GLY",
    "H":"HIS",
    "I":"ILE",
    "K":"LYS",
    "L":"LEU",
    "M":"MET",
    "N":"ASN",
    "P":"PRO",
    "Q":"GLN",
    "R":"ARG",
    "S":"SER",
    "T":"THR",
    "V":"VAL",
    "W":"TRP",
    "Y":"TYR",
}

from contextlib import redirect_stdout, redirect_stderr

def prepare(in_file,out_file,molecule_name,chain_id,mutation):
    with open("redirect_output.out","a+") as file:
        with redirect_stdout(file):
            """Mutate specified residues across chains and generate an output file."""

            #reinitialize to restore program to default state after possible previous calls of this function
            pymol.cmd.reinitialize()
            pymol.cmd.load(in_file, molecule_name)

            if "from_long" not in mutation:
                mutation["from_long"]=aa_table_one_to_three[mutation["from"]]

            if "to_long" not in mutation:
                mutation["to_long"]=aa_table_one_to_three[mutation["to"]]

            mutation["position"]=int(mutation["position"])

            assert mutation["from_long"]==aa_table_one_to_three[mutation["from"]]
            assert mutation["to_long"]==aa_table_one_to_three[mutation["to"]]

            """Retrieve chain IDs. """    
            try:
                ChainIDs = pymol.cmd.get_chains(f"model {molecule_name}")
            except pymol.CmdException as ErrMsg:
                print(f"PyMOLUtil.GetChains: Invalid molecule name: {molecule_name}")
                sys.exit(1)

            if len(ChainIDs)==0:
                print(f"no chains found for molecule {molecule_name}")
                sys.exit(1)
            if not chain_id in ChainIDs:
                print(f"chain {chain_id} not found in chain ids of molecule {molecule_name}")
                sys.exit(1)
            
            # Apply mutations
            
            # Setup mutagenesis wizard
            pymol.cmd.wizard("mutagenesis")
            pymol.cmd.refresh_wizard()

            # select residue to be mutated
            select_residue = "/%s//%s/%s" % (molecule_name, chain_id, mutation["position"])
            pymol.cmd.iterate(select_residue, "(stored.resi,stored.resn)=(resi,oneletter)")
            if int(pymol.stored.resi)!=mutation["position"] or pymol.stored.resn!=mutation["from"]:
                print(f"WT residue or position wrong: expected {mutation}, but structure actually contains {pymol.stored.resn} at {pymol.stored.resi}")
                sys.exit(1)

            pymol.cmd.get_wizard().do_select(select_residue)

            # setup mutation (create mutation and open list of rotamers)
            pymol.cmd.get_wizard().set_mode(mutation["to_long"])

            # get rotamer with lowest strain
            strain=pymol.cmd.get_wizard().bump_scores[0]
            least_strain_rotamer_index=0
            for i in range(0,len(pymol.cmd.get_wizard().bump_scores)):
                if pymol.cmd.get_wizard().bump_scores[i]<strain:
                    strain=pymol.cmd.get_wizard().bump_scores[i]
                    least_strain_rotamer_index=i
            
            # apply mutation
            pymol.cmd.get_wizard().do_state(least_strain_rotamer_index+1)
            pymol.cmd.get_wizard().apply()
            
            # quit wizard
            pymol.cmd.set_wizard("done")

            select_residue = "/%s//%s/%s" % (molecule_name, chain_id, mutation["position"])
            pymol.cmd.iterate(select_residue, "(stored.resi,stored.resn)=(resi,oneletter)")
            if int(pymol.stored.resi)!=mutation["position"] or pymol.stored.resn!=mutation["to"]:
                print(f"WT residue or position wrong: expected mutation {mutation}, but got {pymol.stored.resn} in position {pymol.stored.resi}")
                sys.exit(1)

            #  save to output file
            pymol.cmd.save(out_file, molecule_name)
            
            # pymol.cmd.delete(MolName)

if __name__=="__main__":
    mutation={
        "from":"Q",
        "from_long":"GLN",
        "to":"R",
        "to_long":"ARG",
        "position":54,
    }

    molecule_name="fgf2"
    chain_id="A" # from testing. automate how?

    in_file=sys.argv[1]
    print("mutating structure: ", in_file)
    if len(sys.argv)>2:
        out_file=sys.argv[2]
    else:
        out_file=sys.argv[1][:-4]+f"_mutated_{mutation['from']}{mutation['position']}{mutation['to']}.pdb"

    prepare(in_file,out_file,molecule_name,chain_id,mutation)
