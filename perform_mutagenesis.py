#!/bin/env python3

mutation={
    "from":"Q",
    "from_long":"GLN",
    "to":"R",
    "to_long":"ARG",
    "position":54,
}
out_directory_base="mutation_performed"
folder_prefix=""

molecule_name="fgf2"
chain_id="A" # from testing. automate how?

#
# File: PyMOLMutateAminoAcids.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2021 Manish Sud. All rights reserved.
#
# The functionality available in this script is implemented using PyMOL, a
# molecular visualization system on an open source foundation originally
# developed by Warren DeLano.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.

#https://pymolwiki.org/index.php/Launching_From_a_Script # how to use pymol as cli
#https://pymolwiki.org/index.php/Iterate # how to get information from selected atoms
#https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/wizard/mutagenesis.py # pymol mutagenesis internals, including bump_scores (strain value for rotamers), and an example on pymol.storage for pymol interal->external communication
#http://www.mayachemtools.org/docs/scripts/html/code/PyMOLMutateAminoAcids.html # the basis for this script, though not much is left
#http://www.mayachemtools.org/docs/modules/html/PyMOLUtil.py.html#GetChains # reference for pymol chain retrieval, though not much is left
#http://www.endmemo.com/bio/codon.php # codon table because i cannot remember them myself

import sys

try:
    import pymol
    # wait for pymol to finish launching cli (-c), in quiet mode (-q), and disable loading of pymolrc and plugins (-k)
    pymol.finish_launching(['pymol', '-ckq'])
except ImportError as ErrMsg:
    sys.stderr.write(f"Failed to import PyMOL package, most likely because pymol is not installed: {ErrMsg}\n")
    sys.exit(1)

in_file=sys.argv[1]
print("mutating structure: ", in_file)
if len(sys.argv)>2:
    out_file=sys.argv[2]
else:
    out_file=sys.argv[1][:-4]+f"_mutated_{mutation['from']}{mutation['position']}{mutation['to']}.pdb"

def main():
    """Mutate specified residues across chains and generate an output file."""

    pymol.cmd.reinitialize() # necessary?
    pymol.cmd.load(in_file, molecule_name) # LoadMolecule

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
        print(f"WT residue or position wrong: mutation={{{mutation}}}, but actually residue_index: {pymol.stored.resi} and residue name: {pymol.stored.resn}")
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
        print(f"WT residue or position wrong: mutation={{{mutation}}}, but actually residue_index: {pymol.stored.resi} and residue name: {pymol.stored.resn}")
        sys.exit(1)

    #  save to output file
    pymol.cmd.save(out_file, molecule_name)
    
    # pymol.cmd.delete(MolName)

main()