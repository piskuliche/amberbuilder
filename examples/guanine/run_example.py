from amberbuilder import Builder
import glob

leaprc=['source leaprc.RNA.OL3',
        'source leaprc.gaff2']

NewBuild = Builder(boxshape="octahedral",
                    box_buffer=10, 
                    neutralize=True, 
                    ion_concentration=0.14,
                    solvent="tip4pew", leaprc=leaprc)

targets = glob.glob("targets/*.pdb")
for target in targets:
    NewBuild.add_target(target)

NewBuild.build()