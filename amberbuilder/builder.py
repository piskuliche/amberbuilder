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
    def __init__(self, boxshape="orthorhombic", com_buffer=10, aq_buffer=16, neutralize=True, ion_concentration=0.14, 
                 add_na=0, add_cl=0, solvent='tip3p', density=0.997, leaprc=[], nucleic=False, verbose=False):
        """ This class builds a consistent set of boxes for amber simulations across a wide range of targets.
        
        Parameters
        ----------
        boxshape : str
            The shape of the box. Either 'orthorhombic' or 'octahedral'. [default='orthorhombic']
        com_buffer : float
            The buffer between the solute and the edge of the box. [default=10]
        aq_buffer : float
            The buffer between the solute and the edge of the aqueous box. [default=16]
        neutralize : bool
            Whether to neutralize the system. [default=True]
        ion_concentration : float
            The concentration of ions in the box. [default=0.14]
        add_na : int
            The number of additional NA ions to add. [default=0]
        add_cl : int
            The number of additional CL ions to add. [default=0]
        solvent : str
            The solvent to use. Either 'tip3p' or 'tip4pew'. [default='tip3p']
        density : float
            The density of the solvent. [default=0.997]
        leaprc : list
            A list of leaprc files to use. [default=[]]
        nucleic : bool
            Whether the system is a nucleic acid. This impacts neutralization. [default=False]
        verbose : bool
            Whether to print verbosely. [default=False]
        
        """
        self.boxshape = str(boxshape)
        self.com_buffer = float(com_buffer)
        self.aq_buffer = float(aq_buffer)
        self.neutralize = bool(neutralize)
        self.ion_concentration = float(ion_concentration)
        self.add_na = int(add_na)
        self.add_cl = int(add_cl)
        self.solvent = str(solvent)
        self.density = float(density)
        self.leaprc = leaprc
        self.nucleic = nucleic
        self.verbosity = verbose

        # Things that are set later
        self.target_files = []
        self.target_universes = []
        self.rotations = None
        self.aligned_superuniverse = None

        # Test the inputs
        self._test_input()
        pass

    def build(self, aqueous=False):
        """ This builds a set of amber simulation files using a shell thickness approach.

        This build calls the shell_build function to build a set of amber simulation files using a shell thickness approach.
        This also calls the aqueous build function to build a set of amber simulation files in the aqueous phase.
        
        Parameters
        ----------
        aqueous : bool
            Whether to build the system in the aqueous phase. [default=False]
        
        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            If an error occurs while building the Amber simulation files.
        
        """
        print("Building a set of Amber simulation files")
        try:
            self._shell_build(aqueous=aqueous)
        except Exception as e:
            raise RuntimeError(f"An error occurred while building the Amber simulation files: \n {e}")
        return


    def add_target(self, target):
        """ Add a target to the list of targets.

        This function takes a target (in some initial directory) and then it adds it to the list of targets to build. It also 
        copies the target, and its associated lib and frcmod files to the targets directory. If the target is a nucleic acid, 
        it will add the neutralizing ions to the system prior to moving the target to the targets directory.
        
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
            If one of the target files is not present.
        
        """
        import shutil
        print(f"Adding target {target}")
        # Check that necessary files exist.
        pdb = Path(target)
        lib_file = Path(target.replace(".pdb", ".lib"))
        frcmod_file = Path(target.replace(".pdb", ".frcmod"))
        for file in [pdb, lib_file, frcmod_file]:
            if not file.exists():
                raise FileNotFoundError(f"Could not find {file}")
        
        # Create the target directory - this is where files are copied TO
        target_dir = Path("targets/")
        if not target_dir.exists():
            target_dir.mkdir()
        
        # Copy the files to the target directory
        if not self.nucleic: # If it's not a nucleic acid, just copy the files
            shutil.copy(pdb, target_dir / pdb.name)
        else: # If it is a nucleic acid, neutralize it first - this uses the self.add_na and self.add_cl values
            self.write_nucleic_neutralize(str(pdb))

        # Add the target to the list of targets
        self.target_files.append(str(target_dir / pdb.name))

        # Copy the lib and frcmod files to the target directory
        shutil.copy(lib_file, target_dir / lib_file.name)
        shutil.copy(frcmod_file, target_dir / frcmod_file.name)

        return
    
    def _shell_build(self, aqueous=False):
        """ Build amber simulation files using a shell thickness approach.

        This function builds the amber simulation files around a common target using a shell thickness approach. It reads the
        target files into a superuniverse, translates the superuniverse to the center of geometry, and then  uses TLEAP to pack
        the box with water and ions. The ligands are then aligned to the common scaffold, and written out as separate pdb, parm, and rst7 files.
        
        Parameters
        ----------
        aqueous : bool
            Whether to build the system in the aqueous phase. [default=False]
        
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
        
        if aqueous:
            self.box_buffer = self.aq_buffer
        else:
            self.box_buffer = self.com_buffer
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            print("*********************************************")

            # Read the targets, and creates a superuniverse, and removes the COG
            self.superuniverse = self._target_reader(aqueous=aqueous)
            mdatools.WritePDB(self.superuniverse, "superuniverse.pdb")
            print("*********************************************")
            # Pack the box using TLEAP and add ions
            self._pack_box("superuniverse.pdb")
            print("*********************************************")
            # Split the targets and apply the rotation to each
            self._remove_and_align()
            self._write_empty_target()
            self._combine_into_empty(aqueous=aqueous)
            print("*********************************************")
            # Clean up the intermediate files

            self._clean(aqueous=aqueous)

            if self.verbosity:
                print("*********************************************")
                print("MDAnalysis warnings:")
                for warning in w:
                    print(f"Warnings: {warning.message}")
            print("*********************************************")
            print("Finished building the Amber simulation files")
    
        
    def _aqueous_build(self):
        """ Build amber aqueous simulation files using a shell thickness approach.

        This function builds the amber simulation files for the ligands in the aqueous phase using a shell thickness approach. It reads the
        target files into a superuniverse, translates the superuniverse to the center of geometry, and then uses TLEAP to pack the box with water
        and ions. The ligands are then aligned to the common scaffold, and written out as separate pdb, parm, and rst7 files.
        
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
        raise NotImplementedError("Aqueous build not yet implemented")
       

    
    def _target_reader(self, aqueous=False):
        """ This function reads the target files individually as universes and uses them to build a superuniverse.
        
        This function reads the target files into a superuniverse. These assume the first target file as the reference, and then 
        adds the ligands from the other target files to the superuniverse. The additional ligands are renamed to REM, and the superuniverse is
        returned. THUS, any neutralizing ions really only need be included in the FIRST target file, as those are the only ones that will be added to 
        the superuniverse.
        
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
        self.target_universes = []
        for i, target in enumerate(self.target_files):
            u = mda.Universe(target)
            self.target_universes.append(u)
            if i == 0 and not aqueous:
                group_list.append(u.select_atoms("all"))
            elif i == 0 and aqueous:
                group_list.append(u.select_atoms("resname LIG"))
            else:
                if "LIG" in u.residues.resnames:
                    ag = u.select_atoms("resname LIG")
                    ag.residues.resnames = [f"REM"]
                    group_list.append(ag)
                else:
                    raise ValueError(f"Error: Ligand not found in target file {target}")
                
        superuniverse = mda.Merge(*group_list)
        print(f"Read {len(self.target_files)} target files into a superuniverse with {superuniverse.atoms.n_atoms} atoms.")
        return superuniverse


    def _pack_box(self, pdbfilename):
        """ This function packs the box with water and ions using the TLEAP program.

        This function packs the box with water and ions using the TLEAP program. It first solvates the box with water, and then adds ions
        to the system based on the concentration. This happens in two subsequent steps, first it solvates the box with tleap, then it calculates the number of
        water molecules to add ions to the system, and then it adds the ions to the system.
        
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
        print(f"Adding {num_NA} NA ions and {num_CL} CL ions to the system.")

        # Use tleap to add the ions.
        newleap = []
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"mol = loadpdb box.pdb")
        newleap.append(f"addionsrand mol Na+ {int(num_NA)} Cl- {int(num_CL)} 6.0")
        newleap.append("savepdb mol packed.pdb")
        newleap.append("quit")
        self.write_file("tleap.ions.in", newleap)

        leap = Leap()
        leap.call(f="tleap.ions.in")
        return 
    
    def _remove_and_align(self):
        """ This function removes the ligands from the superuniverse and aligns them to the common scaffold.

        This function removes the ligands from the superuniverse and aligns them to the common scaffold. It then writes the aligned ligands
        to pdb files.

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
        # Create a universe from the packed box
        u = mda.Universe("packed.pdb")

        # Get JUST the ligand to ALIGN TO and read it into rdkit
        just_lig = u.select_atoms("resname LIG")
        mdatools.WritePDB(just_lig, "ligand.pdb")
        goal = Chem.MolFromPDBFile("ligand.pdb", removeHs=True)

        # Loop over targets
        for i, target in enumerate(self.target_universes):
            # Get the ligand TO ALIGN from the universe and read it into rdkit
            target_lig = target.select_atoms("resname LIG REM")
            mdatools.WritePDB(target_lig, f"ligand_{i}.pdb")
            tgt = Chem.MolFromPDBFile(f"ligand_{i}.pdb", removeHs=True)
            
            # Find the maximum common substructure
            mcs = rdFMCS.FindMCS([goal, tgt])
            common_smarts = mcs.smartsString
            common_mol = Chem.MolFromSmarts(common_smarts)
            
            # Map the atoms to find the common atoms
            match1 = goal.GetSubstructMatch(common_mol)
            match2 = tgt.GetSubstructMatch(common_mol)

            # Align the molecules based on the common scaffold.
            AllChem.AlignMol(tgt, goal, atomMap=list(zip(match2, match1)))

            # Write the aligned ligand to a pdb file
            Chem.MolToPDBFile(tgt, f"aligned_{i}.pdb")
        return
        
    def _write_empty_target(self):
        """ Writes a version of the target file with the ligands removed.

        This function writes a version of the target file with the ligands removed. This is used to combine the aligned ligands into a single
        target file.
        
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
            
    def _combine_into_empty(self, aqueous=False):
        """ Combine the aligned ligands into an empty target file.

        This function combines the aligned ligands into an empty target file. This is used to create the final complex.
        
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
            self._write_final_tleap(i, aqueous=aqueous)
        return
    
    def _write_final_tleap(self, i, aqueous=False):
        """ Write the final tleap file for the system.

        This function writes the final tleap file for the system. This file is used to create the final amber simulation files, including
        the parm and rst7 files.
        
        Parameters
        ----------
        i : int
            The index of the target file.
        aqueous : bool
            Whether the system is in the aqueous phase. [default=False]
        
        Returns
        -------
        None
        
        Raises
        ------
        FileNotFoundError
            If the lib or frcmod files are not found.
        
        """
        subfolder = "com" if not aqueous else "aq"
        molname = self.target_files[i].replace('.pdb', '').replace('targets/', '')
        name = self.target_files[i].replace('.pdb', '').replace('targets/', f'outputs/{molname}/{subfolder}/')
        lib_file = self.target_files[i].replace('.pdb', '.lib')
        frcmod_file = self.target_files[i].replace('.pdb', '.frcmod')

        # Check that the lib and frcmod files exist
        for file in [lib_file, frcmod_file]:
            if not Path(file).exists():
                raise FileNotFoundError(f"Could not find {file}")  
        
        outputpath = Path(f"outputs/{molname}/{subfolder}/")
        if not outputpath.exists():
            outputpath.mkdir(parents=True)
        
        # Write the final tleap file
        newleap = []
        newleap.extend(self.leaprc)
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"loadoff {lib_file}")
        newleap.append(f"loadamberparams {frcmod_file}")
        newleap.append(f"mol = loadpdb complex_{i}.pdb")
        if not self.nucleic:
            if self.add_na > 0:
                newleap.append(f"addions mol Na+ {self.add_na}")
            if self.add_cl > 0:
                newleap.append(f"addions Cl- {self.add_cl}")
        newleap.append(f"saveamberparm mol {name}.parm7 {name}.rst7")
        newleap.append("quit")
        self.write_file(f"tleap.final_{i}.in", newleap)

        leap = Leap()
        leap.call(f=f"tleap.final_{i}.in")
        return
    
    def write_nucleic_neutralize(self, pdbfilename):
        """ Write the tleap file for a nucleic acid system.

        This script writes a tleap file to neutralize a nucleic acid system. It first loads the pdb file, and then adds the ions to the system.
        
        Parameters
        ----------
        i : int
            The index of the target file.
        
        Returns
        -------
        None
        
        Raises
        ------
        None
        
        """
        pdbfilename = str(pdbfilename)
        lib_file = pdbfilename.replace('.pdb', '.lib')
        frcmod_file = pdbfilename.replace('.pdb', '.frcmod')
        outputpdb = pdbfilename.replace('initial/', 'targets/')
        for file in [lib_file, frcmod_file]:
            if not Path(file).exists():
                raise FileNotFoundError(f"Could not find {file}")
        
        newleap = []
        newleap.extend(self.leaprc)
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"loadoff {lib_file}")
        newleap.append(f"loadamberparams {frcmod_file}")
        newleap.append(f"mol = loadpdb {pdbfilename}")
        newleap.append(f"addions mol Na+ {self.add_na}")
        newleap.append(f"savepdb mol {outputpdb}")
        newleap.append("quit")

        tleap_file = f"tleap.nuc_neut.in"
        with open(tleap_file, "w") as W:
            W.write("\n".join(newleap))

        leap = Leap()
        leap.call(f=tleap_file)

        return
    
    def _clean(self, aqueous=False):
        """ Cleans up all of the intermediate files created during the build process.

        This function cleans up all of the intermediate build files created during the build process and sorts them
        into the appropriate directories.
        
        Parameters
        ----------
        aqueous : bool
            Whether the system is in the aqueous phase. [default=False]
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If the subfolder is not 'com' or 'aq'.
        
        """
        print("Cleaning up intermediate files")

        subfolder = "com" if not aqueous else "aq"

        temporary = Path("temporary/")
        key = Path(f"temporary/keyfiles/{subfolder}")
        aligned = Path(f"temporary/aligned/{subfolder}")
        ligands = Path(f"temporary/ligands/{subfolder}")
        tleap = Path(f"temporary/tleap/{subfolder}")
        for folder in [key, aligned, ligands, tleap]:
            if not folder.exists():
                folder.mkdir(parents=True)
        
        # Move the key files
        key_files = ["box.pdb", "packed.pdb", "superuniverse.pdb", "ligand.pdb", "empty_target.pdb"]
        for key_file in key_files:
            if Path(key_file).exists():
                Path(key_file).replace(key / key_file)
        
        # Move the aligned ligands
        aligned_files = Path().glob("aligned_*.pdb")
        for aligned_file in aligned_files:
            aligned_file.replace(aligned /  aligned_file)
        
        # Move the ligands
        ligand_files = Path().glob("ligand_*.pdb")
        for ligand_file in ligand_files:
            ligand_file.replace(ligands / ligand_file)
        
        # Move the tleap files
        tleap_files = Path().glob("tleap.*.in")
        for tleap_file in tleap_files:
            tleap_file.replace(tleap / tleap_file)

        complex_files = Path().glob("complex_*.pdb")
        for complex_file in complex_files:
            complex_file.unlink()
        
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
        assert self.com_buffer >= 0, "Box buffer must be a positive number"
        assert self.aq_buffer >= 0, "Aqueous buffer must be a positive number"
        assert self.ion_concentration >= 0, "Ion concentration must be a positive number"
        return
    

if __name__ == "__main__":
    builder = Builder()
    builder.build()