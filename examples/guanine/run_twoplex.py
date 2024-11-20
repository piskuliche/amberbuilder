from amberbuilder import Builder
from amberbuilder.rbfe_tools import rbfe_prep
import glob


leaprc=['source leaprc.RNA.OL3',
        'source leaprc.gaff2',
        'source leaprc.water.tip4pew',]

setup_rbfe = rbfe_prep(leaprc=leaprc)
setup_rbfe.add_edge("DOG", "NNG")
setup_rbfe.add_edge("DOG", "NOG")
setup_rbfe.prep_edges()