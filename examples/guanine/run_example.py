from amberbuilder import Builder
import glob

NewBuild = Builder(boxshape="orthorhombic",
                    box_buffer=10, 
                    neutralize=True, 
                    ion_concentration=0.14,
                    solvent="tip4pew")

targets = glob.glob("targets/*.pdb")
for target in targets:
    NewBuild.add_target(target)

NewBuild.build()