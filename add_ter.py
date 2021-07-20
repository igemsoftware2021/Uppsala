#!/bin/env python3

import sys

infilename=sys.argv[1]
outfilename=infilename[:-4]+"_with_term.pdb"

with open(infilename) as infile:
			with open(outfilename,"w+") as outfile:
				for line_number,line in enumerate(infile.readlines()):
			#print(line[13:])
			prepend_line_content="N   ASN   131      63.387  18.498 -44.664  0.00  0.00"
			if line[13:(13+len(prepend_line_content))]==prepend_line_content:
				outfile.write("TER\n")
			outfile.write(line)
