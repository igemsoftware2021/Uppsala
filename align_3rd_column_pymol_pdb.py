#!/bin/env python3
"""
replace the atom type with a properly aligned atom type according to lucien/hugos instructions to work with protprep/qresfep
"""
import sys
infilename=sys.argv[1]
outfilename=infilename[:-4]+"_pymol_naming_fixed.pdb" if len(sys.argv)==2 else sys.argv[2]
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
                        atom_type=line[12:16]
                        if (atom_type[0]!=" " or (atom_type[1]!=" " and atom_type[0]==" ") or (atom_type[2]!=" " and atom_type[0]==" " and atom_type[1]==" ")) and atom_type[0].isnumeric():
                                print("first character of atom type is not space:",line_number,atom_type)

				# move number to back of name
                                if line[15]==" ":
                                        #line[15]=line[12]
                                        line=line.replace(atom_type," "+atom_type[1:3]+atom_type[0])
                                elif line[16]==" ":
                                        #line[16]==line[12]
                                        line=line.replace(line[12:17]," "+atom_type[1:4]+atom_type[0])

                                #line[12]=" "
                                atom_type=line[12:17]
                                print("fixed line:", line_number+1,atom_type)
                        outfile.write(line)
