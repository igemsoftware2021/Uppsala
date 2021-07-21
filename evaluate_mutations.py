complex_pdb="FGF+R2_Q.pdb"
single_pdb="FGF2_Q.pdb"
mutations="Q54R Q54K Q56M F93Y V88I F93L F93W F93Y L98M".split(" ")

import setup as setupfile

setup=setupfile.Wrapper(directory="test_out",complex_pdb=complex_pdb,single_pdb=single_pdb,mutations=mutations)
setup.build()
setup.prepare()
setup.run()
setup.analyze()
