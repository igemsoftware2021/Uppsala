#!/bin/env python3
"""
replace the atom type with a properly aligned atom type according to lucien/hugos instructions to work with protprep/qresfep
"""
import sys
infilename=sys.argv[1]
outfilename=infilename[:-4]+"_pymol_naming_shifted.pdb" if len(sys.argv)==2 else sys.argv[2]
with open(infilename) as infile:
        with open(outfilename,"w+") as outfile:
                for line_number,line in enumerate(infile.readlines()):
                        line_name=line[:4]
                        if line_name!="ATOM":
                                print("line type is not atom:",line_number+1,line)
                                outfile.write(line)
                                continue
                        assert line[4]==" "
                        atom_number=line[5:11]
                        atom_number=int(atom_number)
                        assert line[11]==" "
                        atom_type=line[12:17]
                        if (atom_type[0]!=" "):
                                print("first character of atom type is not space:",line_number,atom_type)
                                outfile.write(line.replace(atom_type," "+atom_type[:-1]))
                                atom_type=line[12:17]
                                print("fixed line:", line_number+1,atom_type)
                        else:
                                outfile.write(line)
