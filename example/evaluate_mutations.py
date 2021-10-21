#the structure files need to be prepared for use with Q first!
complex_pdb="FGF+R2_Q.pdb"
single_pdb="FGF2_Q.pdb"
mutations="Q54K E58L"

import mutapipe
import sys

# create a configuration that contains all initially-relevant information for mutapipe, on a local computer
setup=mutapipe.Wrapper(directory="test_out",complex_pdb=complex_pdb,single_pdb=single_pdb,mutations=mutations.split(" "))
# then run initial preparation steps: mutapipe.build() and mutapipe.prepare()
if len(sys.argv)==1:
    print("preparing")
    setup.build()
    setup.prepare()
# then copy the whole directory to UPPMAX/the supercomputer
elif sys.argv[1]=="run":
    print("running")
    # and run mutapipe.run() there
    setup.run()
# manually wait until all jobs on the supercomputer are finished
# copy the whole directory from the supercomputer to the local computer
elif sys.argv[1]=="analyze":
    print("analyzing")
    # and run mutapipe.analyze()
    # this will print the relevant results to your command line (in less than a minute, so do not worry if you forget to write them down the first time)
    setup.analyze()