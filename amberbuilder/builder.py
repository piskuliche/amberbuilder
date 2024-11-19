import warnings
import MDAnalysis as mda
import numpy as np


from pathlib import Path
from scipy.linalg import eigh
from MDAnalysis.lib import transformations, mdamath
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

from amberbuilder.interfaces import AddToBox, Leap
from amberbuilder import mdatools


class Builder:
    def __init__(self, boxshape="orthorhombic", box_buffer=10, neutralize=True, ion_concentration=0.14, add_na=0, add_cl=0, solvent='tip3p', density=0.997, leaprc=[]):
        """ This class builds a consistent set of boxes for amber simulations across a wide range of targets.
        
        Parameters
        ----------
        boxshape : str
            The shape of the box. Either 'orthorhombic' or 'octahedral'.
        box_buffer : float
            The buffer between the solute and the edge of the box.
        neutralize : bool
            Whether to neutralize the system.
        ion_concentration : float
            The concentration of ions in the box.
        add_na : int
            The number of additional NA ions to add.
        add_cl : int
            The number of additional CL ions to add.
        solvent : str
            The solvent to use. Either 'tip3p' or 'tip4pew'.
        density : float
            The density of the solvent.
        leaprc : list
            A list of leaprc files to use.
        
        """
        self.boxshape = str(boxshape)
        self.box_buffer = float(box_buffer)
        self.neutralize = bool(neutralize)
        self.ion_concentration = float(ion_concentration)
        self.add_na = int(add_na)
        self.add_cl = int(add_cl)
        self.solvent = str(solvent)
        self.density = float(density)
        self.leaprc = leaprc

        # Things that are set later
        self.target_files = []
        self.target_universes = []
        self.rotations = None
        self.aligned_superuniverse = None

        # Test the inputs
        self._test_input()
        pass

    def build(self):
        """ Build a set of Amber simulation files.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None

        Raises
        ------
        ValueError
            If no targets have been added.
        
        """

        print("Building a set of Amber simulation files")
        try:
            self._shell_build()
        except Exception as e:
            raise RuntimeError(f"An error occurred while building the Amber simulation files: \n {e}")
        return


    def add_target(self, target):
        """ Add a target to the list of targets.
        
        Parameters
        ----------
        target : str
            The path to the target file.
        
        Returns
        -------
        None
        
        Raises
        ------
        FileNotFoundError
            If the target file is not found.
        
        """
        print(f"Adding target {target}")
        if Path(target).exists():
            self.target_files.append(target)
        else:
            raise FileNotFoundError(f"Target {target} not found")
        return
    
    def _shell_build(self):
        """ Build amber simulation files using a shell thickness approach.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If no targets have been added.
        
        """ 
        print("Building the Amber simulation files")
        if len(self.target_files) == 0:
            raise ValueError("_build: No targets have been added")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            print("*********************************************")
            # Read the targets, and creates a superuniverse, and removes the COG
            self.superuniverse = self._altreader()

            self.super_cog = self.superuniverse.atoms.center_of_geometry()
            self.superuniverse.atoms.positions -= self.super_cog
            mdatools.WritePDB(self.superuniverse, "superuniverse.pdb")
            print("*********************************************")
            # Rotates the superuniverse to the principal axes
            print("*********************************************")
            # Pack the box using TLEAP and add ions
            self._pack_box("superuniverse.pdb")
            print("*********************************************")
            # Split the targets and apply the rotation to each
            self._remove_and_align()
            self._write_empty_target()
            self._combine_into_empty()
            print("*********************************************")
            print("MDAnalysis warnings:")
            for warning in w:
                print(f"Warnings: {warning.message}")
            print("*********************************************")
            print("Finished building the Amber simulation files")

    def _targetreader(self):
        """ This function reads the target files individually as universes and uses them to build a superuniverse.

        This function reads the target files individually as universes and uses them to construct a superuniverse.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        superuniverse : mda.Universe
            The superuniverse containing all of the target files.
        
        Raises
        ------
        AssertionError
            If the number of atoms in the superuniverse does not match the sum of the atoms in the target files.
        
        """
        print("Reading the target files into a superuniverse.")
        group_list = []
        total_atoms = 0
        for target in self.target_files:
            self.target_universes.append(mda.Universe(target))
            group_list.append(self.target_universes[-1].select_atoms("all"))
            total_atoms += self.target_universes[-1].atoms.n_atoms
            
        superuniverse = mda.Merge(*group_list)
        super_natoms = superuniverse.atoms.n_atoms

        # Check assertions
        assert super_natoms == total_atoms, f"Error: {super_natoms} atoms in superuniverse, {total_atoms} atoms in target files"
        print(f"Read {len(self.target_files)} target files into a superuniverse with {super_natoms} atoms.")
        return superuniverse
    
    def _altreader(self):
        """ This function reads the target files individually as universes and uses them to build a superuniverse.
        
        This function reads the target files individually as universes and uses them to construct a superuniverse.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        superuniverse : mda.Universe
            The superuniverse containing all of the target files.
            
        Raises
        ------
        AssertionError
            If the number of atoms in the superuniverse does not match the sum of the atoms in the target files.
        
        """
        print("Reading the target files into a superuniverse.")
        group_list = []
        for i, target in enumerate(self.target_files):
            u = mda.Universe(target)
            self.target_universes.append(u)
            if i == 0:
                group_list.append(u.select_atoms("all"))
            else:
                if "LIG" in u.residues.resnames:
                    ag = u.select_atoms("resname LIG")
                    ag.residues.resnames = [f"REM"]
                    group_list.append(ag)
                else:
                    raise ValueError(f"Error: Ligand not found in target file {target}")
                
        superuniverse = mda.Merge(*group_list)
        print(superuniverse.atoms)
        print(f"Read {len(self.target_files)} target files into a superuniverse with {superuniverse.atoms.n_atoms} atoms.")
        return superuniverse

        


    def _pack_box(self, pdbfilename):
        """ This function packs the box with water and ions. 
        
        Parameters
        ----------
        pdbfilename : str
            The name of the pdb file to pack.
        
        Returns
        -------
        None
        
        Raises
        ------
        None 

        """
        newleap = []
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"mol = loadpdb {pdbfilename}")
        if self.boxshape == "orthorhombic":
            newleap.append(f"solvatebox mol {self.solvent.upper()}BOX {self.box_buffer}")
        elif self.boxshape == "octahedral":
            newleap.append(f"solvateoct mol {self.solvent.upper()}BOX {self.box_buffer}")
        newleap.append("savepdb mol box.pdb")
        newleap.append("quit")

        # Write and run leap
        self.write_file("tleap.box.in", newleap)
        leap = Leap()
        leap.call(f="tleap.box.in")

        # Get the number of ions
        num_NA, num_CL = mdatools.GetNumIons("box.pdb", self.ion_concentration)
        newleap = []
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"mol = loadpdb box.pdb")
        newleap.append(f"addionsrand mol Na+ {int(num_NA)} Cl- {int(num_CL)} 6.0")
        newleap.append("savepdb mol packed.pdb")
        newleap.append("quit")
        with open("tleap.ions.in", "w") as W:
            W.write("\n".join(newleap))

        leap = Leap()
        leap.call(f="tleap.ions.in")
        return

    def _split_targets(self):
        """ This function splits the targets to separate pdb files.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Raises
        ------
        None
        
        """
        print("Splitting the targets")
        u = mda.Universe("packed.pdb")
        atoms = u.select_atoms("resname WAT NA+ CL-")
        print(atoms)
        for i, target in enumerate(self.target_universes):
            target.atoms.positions -= self.super_cog
            newu = mda.Merge(atoms, target.atoms)
            mdatools.WritePDB(newu, f"target_{i}.pdb")
        return
    
    def _remove_and_align(self):
        """ This function removes the ligand """
        u = mda.Universe("packed.pdb")
        just_lig = u.select_atoms("resname LIG")
        mdatools.WritePDB(just_lig, "ligand.pdb")
        goal = Chem.MolFromPDBFile("ligand.pdb", removeHs=True)

        for i, target in enumerate(self.target_universes):
            print(np.unique(target.residues.resnames))
            target_lig = target.select_atoms("resname LIG REM")
            mdatools.WritePDB(target_lig, f"ligand_{i}.pdb")
            tgt = Chem.MolFromPDBFile(f"ligand_{i}.pdb", removeHs=True)
            
            # Get the MCS
            mcs = rdFMCS.FindMCS([goal, tgt])
            common_smarts = mcs.smartsString
            common_mol = Chem.MolFromSmarts(common_smarts)
            
            # Get the atom mapping
            match1 = goal.GetSubstructMatch(common_mol)
            match2 = tgt.GetSubstructMatch(common_mol)

            # Align the molecules based on the common scaffold.
            AllChem.AlignMol(tgt, goal, atomMap=list(zip(match2, match1)))

            Chem.MolToPDBFile(tgt, f"aligned_{i}.pdb")
        
    def _write_empty_target(self):
        """ Write an empty target file to the directory.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Raises
        ------
        None
        
        """
        u = mda.Universe("packed.pdb")
        atoms = u.select_atoms("all and not resname LIG REM")
        mdatools.WritePDB(atoms, "empty_target.pdb")
            
    def _combine_into_empty(self):
        """ Combine the aligned ligands into an empty target file.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Raises
        ------
        None
        
        """
        base_target = mda.Universe("empty_target.pdb")
        for i in range(len(self.target_universes)):
            u = mda.Universe(f"aligned_{i}.pdb")
            u.residues.resnames = [f"LIG"]
            combo_target = mda.Merge(base_target.atoms, u.atoms)
            mdatools.WritePDB(combo_target, f"complex_{i}.pdb")
            self._write_final_tleap(i)
        return
    
    def _write_final_tleap(self, i):
        """ Write the final tleap file for the system.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Raises
        ------
        None
        
        """
        lib_file = self.target_files[i].replace('.pdb', '.lib')
        frcmod_file = self.target_files[i].replace('.pdb', '.frcmod')
        if not Path(lib_file).exists():
            raise FileNotFoundError(f"Could not find {lib_file}")
        if not Path(frcmod_file).exists():
            raise FileNotFoundError(f"Could not find {frcmod_file}")    
        
        newleap = []
        newleap.extend(self.leaprc)
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"loadoff {lib_file}")
        newleap.append(f"loadamberparams {frcmod_file}")
        newleap.append(f"mol = loadpdb complex_{i}.pdb")

        if self.add_na > 0:
            newleap.append(f"addions mol Na+ {self.add_na}")
        if self.add_cl > 0:
            newleap.append(f"addions Cl- {self.add_cl}")
        newleap.append(f"saveamberparm mol complex_{i}.parm7 complex_{i}.rst7")
        newleap.append("quit")
        with open(f"tleap.final_{i}.in", "w") as W:
            W.write("\n".join(newleap))

        leap = Leap()
        leap.call(f=f"tleap.final_{i}.in")
        return

    

    @staticmethod
    def write_file(filename, lines):
        """ Write a list of lines to a file.
        
        Parameters
        ----------
        filename : str
            The name of the file to write.
        lines : list
            The lines to write to the file.
        
        Returns
        -------
        None
        
        Raises
        ------
        None
        
        """
        with open(filename, "w") as W:
            W.write("\n".join(lines))
        return
    
    

    def _test_input(self):
        assert self.boxshape in ["orthorhombic", "octahedral"], "Box shape must be either 'orthorhombic' or 'octahedral'"
        assert self.solvent in ["tip3p", "tip4pew"], "Solvent must be either 'tip3p' or 'tip4pew'"
        assert self.density > 0.0, "Density must be greater than 0.5 g/cm^3"
        assert self.density < 2.0, "Density must be less than 2.0 g/cm^3"
        assert self.box_buffer >= 0, "Box buffer must be a positive number"
        assert self.ion_concentration >= 0, "Ion concentration must be a positive number"
        return
    

if __name__ == "__main__":
    builder = Builder()
    builder.build()