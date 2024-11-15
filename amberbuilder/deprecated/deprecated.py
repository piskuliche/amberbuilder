def _pack_ortho_box(self):
        """ This function packs the box with water and ions for an orthorhombic box. """
        print("Packing the box with water and ions")

        # Find the volume fo the superuniverse and calculate the number of water molecules needed
        # to fill the box to the target density.
        target_vol=self._supervolume()
        super_vol = self.ortho_volume(self.superuniverse.dimensions)
        water_vol = super_vol - target_vol
        nwaters = self.water_in_volume(water_vol, density=self.density)

        # Add the water molecules to the box using AddToBox
        print(f"Adding {nwaters} water molecules to the box")
        # Add the solvent
        new_box = AddToBox()
        new_box.call(c="oriented.pdb", a=f"{self.solvent}.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=nwaters)
        test = mda.Universe("packed.pdb")
        wats = test.select_atoms("resname WAT")
        print(f"Added {len(wats.residues)} water molecules to the box")

        # Add the ions to the box using AddToBox
        num_NA, num_CL = self.Get_Num_Ions("packed.pdb", self.ion_concentration)
        new_box.call(c="packed.pdb", a="NA.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=num_NA)
        new_box.call(c="packed.pdb", a="CL.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=num_CL)
        print(f"Added {num_NA} NA ions and {num_CL} CL ions to the box")
        return
    
    def _pack_octa_box(self):
        """ This function packs the box with water and ions for a truncated octahedron using tleap.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Raises
        ------
        NotImplementedError
            If the box shape is octahedral.
        
        """
        """
        # Find the max dimension - this will be the diameter of the truncated octahedron
        max_dim = self.superuniverse.dimensions[:3].max()
        rad = np.round(max_dim / 2)
        # Build an octahedron with the same diameter around a solvent molecule using tleap.
        newleap = []
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"mol = loadpdb {self.solvent}.pdb")
        newleap.append(f"solvateoct mol {self.solvent.upper()}BOX {rad}")
        newleap.append("savepdb mol octbox.pdb")
        newleap.append("quit")
        with open("tleap.in", "w") as W:
            W.write("\n".join(newleap))

        leap = Leap()
        leap.call(f="tleap.in")
        """
        newleap = []
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"mol = loadpdb oriented.pdb")
        newleap.append(f"solvateoct mol {self.solvent.upper()}BOX {self.box_buffer}")
        newleap.append("savepdb mol octbox.pdb")
        newleap.append("quit")
        with open("tleap.octbox.in", "w") as W:
            W.write("\n".join(newleap))
        leap = Leap()
        leap.call(f="tleap.octbox.in")

        # Add ions to the box using tleap.
        num_NA, num_CL = self.Get_Num_Ions("octbox.pdb", self.ion_concentration)
        print(f"Adding {num_NA} NA ions and {num_CL} CL ions to the box")
        newleap = []
        for leap_line in self.leaprc:
            newleap.append(leap_line)
        newleap.append(f"source leaprc.water.{self.solvent}")
        newleap.append(f"mol = loadpdb octbox.pdb")
        newleap.append(f"addionsrand mol Na+ {int(num_NA)} Cl- {int(num_CL)} 6.0")
        newleap.append("savepdb mol packed.pdb")
        newleap.append("quit")
        with open("tleap.ions.in", "w") as W:
            W.write("\n".join(newleap))

        leap = Leap()
        leap.call(f="tleap.ions.in")

        return


    
    def _supervolume(self):
        """ This function calculates the volume of the superuniverse using a grid-based approach.

        This function calculates the volume of the superuniverse using a grid-based approach. The superuniverse
        is gridded with a spacing of 1 Å, and the volume is calculated by counting the number of grid points
        within 1.4 Å of any atom in the superuniverse.
        
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
        grid_spacing = 1.0 # Angstromsself.superuniverse

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


@staticmethod
    def gen_solvent_pdb(solvent):
        """ Generate a pdb file for the solvent and ions."""
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
def _add_more_ions(self):
    """ This function adds additional ions to the box. Currently only NA and CL ions are supported.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    
    Raises
    ------
    NotImplementedError
        If the box shape is octahedral.
        
    """
    print("Neutralizing the system")
    if self.boxshape == "octahedral":
        raise NotImplementedError("Adding more ions to an octahedral box has not been implemented yet")
    new_box = AddToBox()
    new_box.call(c="packed.pdb", a="NA.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=self.add_na)
    new_box.call(c="packed.pdb", a="CL.pdb", o="packed.pdb", RP=2.0, RW=2.0, na=self.add_cl)
