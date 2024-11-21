Example 01: Guanine Riboswitch
==============================

This system uses a simple Guanine riboswitch system to demonstrate the basic functionality of the package. This riboswitch
is a nucleic acid riboswitch (RNA) that binds to guanine. There are a set of six ligands which have been generated to bind to 
this riboswitch, for which we will be generating the topology and the coordinate files needed to run an amber simulation of these
systems.

Learning Outcomes:
------------------

1) Learn how to build consistent boxes for a nucleic acid system.
2) Learn how to build aqueous phase systems for ligands.
3) Learn how to build a system in both octahedral and orthorhombic boxes.


Files 
-----

The files for this example can be found in `amberbuilder/amberbuilder/examples/01_GuanineRiboswitch/` directory of the source code.


Tutorial 
--------

To get started with any amberbuilder run, you need to ensure that you have your ligands fully parametrized and in the correct format. Ideally, your ligands are also all docked into the binding site of your target protein/nucleic acid.
In the current example, we have already done these, and located these files in the initial directory. The necessary files are ligand.pdb ligand.frcmod and ligand.lib where ligand is replaced with the name of your ligand. Within these files
the resname for the ligand needs to be "LIG" across the board.

To build the system, we will be using the Builder class from the amberbuilder package. This class is the main class for building what we refer to as "Nodes", individual ligand + target systems. To start out, we 
will need to import the necessary modules.

.. code-block:: python

    from amberbuilder.builder import Builder
    import glob

Here, we imported both the Builder class that does the hard work, and the glob module which will 
be used to find all the ligands in the initial directory.

Next, we need to set up a leaprc for the system - this is what will be called by tleap. Include everything you 
need outside of the ligands and the solvent in this file. 

.. code-block:: python


    leaprc=['source leaprc.RNA.OL3',
        'source leaprc.gaff2']

Now, we will initialize the builder class with initial parameters.

.. code-block:: python

    NewBuild = Builder(boxshape="orthorhombic",
                        com_buffer=10, 
                        aq_buffer=16,
                        neutralize=True, 
                        ion_concentration=0.14,
                        add_na=67, add_cl=0,
                        solvent="tip4pew",
                        leaprc=leaprc, 
                        nucleic=True)

Here there are a number of command line aguments that we use which influence the calculation. We describe those briefly here:

- boxshape: The shape of the box to be built. This can be either "octahedral" or "orthorhombic".
- com_buffer: The buffer around the target complex to build solvent.
- aq_buffer: The buffer around the solute in the aqueous phase.
- neutralize: Whether to neutralize the system.
- ion_concentration: The concentration of ions in the system.
- add_na: The number of sodium ions to add. We have set this to 67 to neutralize the system.
- add_cl: The number of chloride ions to add.
- solvent: The water model to use.
- leaprc: List of leaprc files to source.
- nucleic: Whether the system is a nucleic acid system.

Note - neutralization behavior is influenced by the choice of add_na, add_cl, and nucleic options. If nucleic is 
selected, then the neutralization happens PRIOR to adding the ions. If not, its added after and replaces some number
of solvent molecules.

Next we get the ligands and add the targets (Nodes) to the builder.

.. code-block:: python

    targets = glob.glob("initial/*.pdb")
    for target in targets:
        NewBuild.add_target(target)

This doesn't do any actual building at this poitn, it just adds the targets and prepares
the builder for the build step.

The build actually happens by running the build method. Since we want both an aqueous system around the ligand
and a complex system, we call this twice. Once to build the smaller aqueous system and once to build the 
larger complex system. This behavior is controlled via the aqueous flag.

.. code-block:: python

    NewBuild.build(aqueous=True)
    NewBuild.build(aqueous=False)

Outputs from the build end up in outputs/ligand/com and outputs/ligand/aq directories. The entire script that does this
is called run_example.py.


To set up RBFE edges with a dual topology, we can use the run_twoplex.py code. This code takes the Nodes
generated above, and then sets up the dual topology system.

.. code-block:: python

    from amberbuilder.builder import Builder
    from amberbuilder.rbfe_tools import rbfe_prep
    import glob


    leaprc=['source leaprc.RNA.OL3',
            'source leaprc.gaff2',
            'source leaprc.water.tip4pew',]

    setup_rbfe = rbfe_prep(leaprc=leaprc)
    setup_rbfe.add_edge("DOG", "NNG")
    setup_rbfe.add_edge("DOG", "NOG")
    setup_rbfe.prep_edges()




Full code
---------

.. code-block:: python

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

.. code-block:: python

    from amberbuilder.builder import Builder
    from amberbuilder.rbfe_tools import rbfe_prep
    import glob


    leaprc=['source leaprc.RNA.OL3',
            'source leaprc.gaff2',
            'source leaprc.water.tip4pew',]

    setup_rbfe = rbfe_prep(leaprc=leaprc)
    setup_rbfe.add_edge("DOG", "NNG")
    setup_rbfe.add_edge("DOG", "NOG")
    setup_rbfe.prep_edges()