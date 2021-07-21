#the structure files need to be prepared for use with Q first! (automate that as well somehow?)
complex_pdb="FGF+R2_Q.pdb"
single_pdb="FGF2_Q.pdb"
#mutations="Q54R Q54K Q56M F93Y V88I F93L F93W F93Y L98M"
mutations="E58L"

import setup as setupfile
import sys

setup=setupfile.Wrapper(directory="test_out",complex_pdb=complex_pdb,single_pdb=single_pdb,mutations=mutations.split(" "))
if len(sys.argv)==1: #cannot be 0, can only be 1 or higher
    print("preparing")
    setup.build()
    setup.prepare()
elif sys.argv[1]=="run":
    print("running")
    setup.run()
elif sys.argv[1]=="analyze":
    print("analyzing")
    setup.analyze()
