import warnings
import MDAnalysis as mda
import numpy as np


from pathlib import Path
from scipy.linalg import eigh

from amberbuilder.interfaces import AddToBox
from MDAnalysis.lib import transformations, mdamath

class Builder:
    def __init__(self, boxshape="orthorhombic", box_buffer=10, neutralize=True, ion_concentration=0.14, add_na=0, add_cl=0, solvent='tip3p', density=0.997):
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
        
        """
        self.boxshape = str(boxshape)
        self.box_buffer = float(box_buffer)
        self.neutralize = bool(neutralize)
        self.ion_concentration = float(ion_concentration)
        self.add_na = int(add_na)
        self.add_cl = int(add_cl)
        self.solvent = str(solvent)
        self.density = float(density)

        # Things that are set later
        self.target_files = []
        self.target_universes = []
        self.rotations = None

        # Test the inputs
        self._test_input()
        self.gen_solvent_pdb(self.solvent)
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
            self._build()
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
    
    def _build(self):
        print("Building the Amber simulation files")
        if len(self.target_files) == 0:
            raise ValueError("_build: No targets have been added")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            print("*********************************************")
            self._targetreader()
            print("*********************************************")
            self._findbox()
            print("*********************************************")
            if self.boxshape == "orthorhombic":
                self._pack_ortho_box()
            elif self.boxshape == "octahedral":
                self._pack_octa_box()
            if self.add_na > 0 or self.add_cl > 0:
                print("*********************************************")
                self._add_more_ions()
            print("*********************************************")
            self._split_targets()
            print("*********************************************")
            print("MDAnalysis warnings:")
            for warning in w:
                print(f"Warnings: {warning.message}")
            print("*********************************************")
            print("Finished building the Amber simulation files")

    def _targetreader(self):
        """ This function reads the target files individually as universes and uses them to build a superuniverse.
        
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
        print("Reading the target files into a superuniverse.")
        group_list = []
        total_atoms = 0
        for target in self.target_files:
            self.target_universes.append(mda.Universe(target))
            group_list.append(self.target_universes[-1].select_atoms("all"))
            total_atoms += self.target_universes[-1].atoms.n_atoms
        self.superuniverse = mda.Merge(*group_list)
        super_natoms = self.superuniverse.atoms.n_atoms
        assert super_natoms == total_atoms, f"Error: {super_natoms} atoms in superuniverse, {total_atoms} atoms in target files"
        with mda.Writer("superuniverse.pdb") as W:
            W.write(self.superuniverse.atoms)
        print(f"Read {len(self.target_files)} target files into a superuniverse with {super_natoms} atoms.")

    def _findbox(self):
        """ Orients the superunvierse along the z-axis and finds the box size.
        
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
        print("Finding the box size")
        self._orient()
        min_coords = self.superuniverse.atoms.positions.min(axis=0)
        max_coords = self.superuniverse.atoms.positions.max(axis=0)
        lengths = max_coords - min_coords + 2 * self.box_buffer
        self.superuniverse.dimensions = [lengths[0], lengths[1], lengths[2], 
                                         90.0, 90.0, 90.0]
        with mda.Writer("oriented.pdb") as W:
            W.write(self.superuniverse.atoms)
        print("Selected box dimensions: ", self.superuniverse.dimensions)
        
    def _orient(self):
        """ This function orients the system along the z-axis.
        
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
        print("Orienting the system")
        #aligned = self.superuniverse.atoms.align_principal_axis(2,[0,0,1])
        p = self.superuniverse.atoms.principal_axes()[2]
        angle = np.degrees(mdamath.angle(p, [0,0,1]))
        ax = transformations.rotaxis(p, [0,0,1])
        aligned = self.superuniverse.atoms.rotateby(angle, ax)

        self.superuniverse = aligned
        self.rotation = [angle, ax]
        return
    
    def _add_more_ions(self):
        print("Neutralizing the system")
        new_box = AddToBox()
        new_box.call(c="packed.pdb", a="NA.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=self.add_na)
        new_box.call(c="packed.pdb", a="CL.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=self.add_cl)
    
    def _pack_ortho_box(self):
        """ This function packs the box with water and ions."""
        print("Packing the box with water and ions")
        target_vol=self._supervolume()
        super_vol = self.ortho_volume(self.superuniverse.dimensions)
        water_vol = super_vol - target_vol
        nwaters = self.water_in_volume(water_vol, density=self.density)

        print(f"Adding {nwaters} water molecules to the box")
        # Add the solvent
        new_box = AddToBox()
        new_box.call(c="oriented.pdb", a=f"{self.solvent}.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=nwaters)
        test = mda.Universe("packed.pdb")
        wats = test.select_atoms("resname WAT")
        print(f"Added {len(wats.residues)} water molecules to the box")

        # Add the ions
        num_NA, num_CL = self.Get_Num_Ions("packed.pdb", self.ion_concentration)
        new_box.call(c="packed.pdb", a="NA.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=num_NA)
        new_box.call(c="packed.pdb", a="CL.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=num_CL)
        print(f"Added {num_NA} NA ions and {num_CL} CL ions to the box")
        return
    
    def _pack_octa_box(self):
        return
    
    def _supervolume(self):
        """ This function calculates the volume of the superuniverse using a grid-based approach.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        float
            The volume of the superuniverse.
        
        Raises
        ------
        None
        
        """
        print("Calculating the superuniverse volume")
        positions = self.superuniverse.atoms.positions
        min_coords = positions.min(axis=0)
        max_coords = positions.max(axis=0)
        grid_spacing = 1.0 # Angstroms

        x = np.arange(min_coords[0], max_coords[0], grid_spacing)
        y = np.arange(min_coords[1], max_coords[1], grid_spacing)   
        z = np.arange(min_coords[2], max_coords[2], grid_spacing)

        grid = np.meshgrid(x, y, z, indexing="ij")
        grid_points = np.vstack([np.ravel(g) for g in grid]).T

        distances = np.linalg.norm(positions[:, np.newaxis, :] - grid_points[np.newaxis, :, :], axis=2)
        within_cutoff = np.any(distances < 1.4, axis=0)  # 1.4 Å is a typical van der Waals radius
        
        volume = np.sum(within_cutoff) * (grid_spacing ** 3)
        print(f"Superuniverse volume: {volume} Å^3")
        return volume

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
            newu = mda.Merge(atoms, target.atoms.rotateby(self.rotation[0], self.rotation[1]))
            with mda.Writer(f"output_target_{i}.pdb") as W:
                W.write(newu.atoms)
        return
    
    @staticmethod
    def gen_solvent_pdb(solvent):
        assert solvent in ["tip3p", "tip4pew"], "Solvent must be either 'tip3p' or 'tip4pew'"
        if solvent == "tip3p":
            lines = """CRYST1  179.034  179.034  179.034 109.47 109.47 109.47 P 1           1
ATOM  1      O   WAT       1    -7.336   0.430   6.058  1.00  0.00
ATOM  2      H1  WAT       1    -7.789  -0.296   5.586  1.00  0.00
ATOM  3      H2  WAT       1    -6.850  -0.048   6.757  1.00  0.00"""
        elif solvent == "tip4pew":
            lines = """CRYST1  179.034  179.034  179.034 109.47 109.47 109.47 P 1           1
ATOM  18705  O   WAT  1868      18.351 -37.529  84.385  1.00  0.00
ATOM  18706  H1  WAT  1868      19.241 -37.676  84.701  1.00  0.00
ATOM  18707  H2  WAT  1868      18.051 -38.393  84.102  1.00  0.00
ATOM  18708 EPW  WAT  1868      18.426 -37.659  84.388  1.00  0.00"""
        with open(f"{solvent}.pdb", "w") as W:
            W.write(lines)
        
        lines = """CRYST1  179.034  179.034  179.034 109.47 109.47 109.47 P 1           1
ATOM  18705  O   NA+  1868      18.351 -37.529  84.385  1.00  0.00"""
        with open(f"NA.pdb", "w") as W:
            W.write(lines)
        lines = """CRYST1  179.034  179.034  179.034 109.47 109.47 109.47 P 1           1
ATOM  18705  O   CL-  1868      18.351 -37.529  84.385  1.00  0.00"""
        with open(f"CL.pdb", "w") as W:
            W.write(lines)
        return
    @staticmethod
    def ortho_volume(dimensions):
        """ Calculate the volume of a orthorhombic box.
        
        Parameters
        ----------
        dimensions : list
            The dimensions of the box.
        
        Returns
        -------
        float
            The volume of the box.
        
        """

        return dimensions[0] * dimensions[1] * dimensions[2]
    
    @staticmethod
    def octa_volume(dimensions):
        """ Calculate the volume of an octahedral box.
        
        Parameters
        ----------
        dimensions : list
            The dimensions of the box.
        
        Returns
        -------
        float
            The volume of the box.
        
        """
        raise NotImplementedError("octa_volume has not been implemented yet")

    @staticmethod
    def water_in_volume(volume, density=0.997):
        """ Calculate the number of water molecules in a volume.
        
        Parameters
        ----------
        volume : float
            The volume of the box.
        density : float
            The density of water.
        
        Returns
        -------
        int
            The number of water molecules in the box.
        
        """
        Na = 6.022e23
        mw = 18.01468
        angpercm = 1e8
        nwaters = int(volume * density * Na / mw * angpercm ** -3)
        return nwaters
    
    @staticmethod
    def Get_Num_Ions(pdb_file, concentration):
        """ Get the number of ions needed for the system. """
        water_concentration = 55.
        u = mda.Universe(pdb_file)
        num_waters = len(u.select_atoms("resname WAT").residues)

        conc_na = ((num_waters ) * concentration / water_concentration) 
        conc_cl = ((num_waters ) * concentration / water_concentration) 

        if conc_na > 0:
            num_NA = int(conc_na)
        else:
            raise ValueError("ERROR: Negative concentration of NA ions")
        if conc_cl > 0:
            num_CL = int(conc_cl)
        else:
            raise ValueError("ERROR: Negative concentration of CL ions")
        return num_NA, num_CL

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