from amberbuilder.builder import Builder
import glob

leaprc=['source leaprc.RNA.OL3',
        'source leaprc.gaff2']

NewBuild = Builder(boxshape="orthorhombic",
                    com_buffer=10, 
                    aq_buffer=16,
                    neutralize=True, 
                    ion_concentration=0.14,
                    add_na=67, add_cl=0,
                    solvent="tip4pew", 
                    leaprc=leaprc, 
                    nucleic=True)

targets = glob.glob("initial/*.pdb")
print(targets)
for target in targets:
    NewBuild.add_target(target)

NewBuild.build(aqueous=True)
NewBuild.build(aqueous=False)
#

