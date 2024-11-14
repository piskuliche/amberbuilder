import warnings
import MDAnalysis as mda


from pathlib import Path
from scipy.linalg import eigh

class Builder:
    def __init__(self, boxshape="cubic", box_buffer=10, neutralize=True, ion_concentration=0.14):
        """ This class builds a consistent set of boxes for amber simulations across a wide range of targets.
        
        Parameters
        ----------
        boxshape : str
            The shape of the box. Either 'cubic' or 'octahedral'.
        box_buffer : float
            The buffer between the solute and the edge of the box.
        neutralize : bool
            Whether to neutralize the system.
        ion_concentration : float
            The concentration of ions in the box.
        
        """
        self.boxshape = str(boxshape)
        self.box_buffer = float(box_buffer)
        self.neutralize = bool(neutralize)
        self.ion_concentration = float(ion_concentration)
        self.target_files = []
        self.target_universes = []
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
            self._packbox()
            print("*********************************************")
            self._neutralize()
            print("*********************************************")
            for warning in w:
                print(f"Warnings: {warning.message}")

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
        super_natoms = self.superuniverse.atoms.n_atomsprint("*********************************************")
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
        aligned = self.superuniverse.atoms.align_principal_axis(2,[0,0,1])
        with mda.Writer("oriented.pdb") as W:
            W.write(aligned.atoms)
        self.superuniverse = aligned
        return
    
    def _neutralize(self):
        print("Neutralizing the system")
        raise NotImplementedError("_neutralize has not been implemented yet")
    
    def _packbox(self):
        print("Packing the box with water and ions")
        
        raise NotImplementedError("_packbox has not been implemented yet")

    def _test_input(self):
        assert self.boxshape in ["cubic", "octahedral"], "Box shape must be either 'cubic' or 'octahedral'"
        assert self.box_buffer >= 0, "Box buffer must be a positive number"
        assert self.ion_concentration >= 0, "Ion concentration must be a positive number"
        return
    

if __name__ == "__main__":
    builder = Builder()
    builder.build()