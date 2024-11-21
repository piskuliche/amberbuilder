import MDAnalysis as mda
import numpy as np

from MDAnalysis.lib import transformations, mdamath

def OrientUniverse(atomgroup):
    """ This function orients an atomgroup along the z-axis.
    The orientation is determined by the principal axis of the superuniverse, and the superuniverse is rotated
    to align the principal axis with the z-axis.
    
    Parameters
    ----------
    atomgroup : mda.AtomGroup
        The atomgroup to orient.
    
    Returns
    -------
    aligned : mda.Universe
        The aligned universe
    angle : float
        The angle of rotation.
    ax : list
        The axis of rotation.
    
    Raises
    ------
    None
    
    """
    print("Orienting the system")

    # Finds the principal axis of the superuniverse and rotates it to align with the z-axis
    p = atomgroup.principal_axes()[2]
    angle = np.degrees(mdamath.angle(p, [0,0,1]))
    ax = transformations.rotaxis(p, [0,0,1])
    aligned = ApplyRotation(atomgroup, angle, ax)

    return aligned, angle, ax


def ApplyRotation(atomgroup, angle, ax):
    """ This function applies a rotation to an atomgroup.
    
    Parameters
    ----------
    atomgroup : mda.AtomGroup
        The atomgroup to rotate.
    angle : float
        The angle of rotation.
    ax : list
        The axis of rotation.
    
    Returns
    -------
    rotated : mda.Universe
        The rotated universe.
    
    Raises
    ------
    None
    
    """
    print("Applying the rotation")
    rotated = atomgroup.rotateby(angle, ax)
    return rotated

def FindBox(atomgroup, box_buffer):
    """ Finds the size of the box.

    This function finds the box size. The box size is determined 
    by the maximum and minimum coordinates of the superuniverse, with a buffer added to each side as specified by
    the box_buffer attribute.
    
    Parameters
    ----------
    atomgroup : mda.AtomGroup
        The atomgroup to find the box size for.
    box_buffer : float
        The buffer to add to each side of the box.
    
    Returns
    -------
    dimensions : list
        The dimensions of the box. (lengths[0], lengths[1], lengths[2], 90.0, 90.0, 90.0)
    
    Raises
    ------
    None
    
    """
    print("Finding the box size")
    min_coords = atomgroup.positions.min(axis=0)
    max_coords = atomgroup.positions.max(axis=0)
    lengths = max_coords - min_coords + 2 * box_buffer
    dimensions = [lengths[0], lengths[1], lengths[2], 
                    90.0,       90.0,       90.0]
    print("Selected box dimensions: ", dimensions)
    return dimensions

def WritePDB(universe, filename):
        """ Write a pdb file from a universe.
        
        Parameters
        ----------
        universe : mda.Universe
            The universe to write.
        filename : str
            The name of the file to write.
        
        Returns
        -------
        None
        
        Raises
        ------
        None
        
        """
        with mda.Writer(filename) as W:
            W.write(universe.atoms)

def GetNumIons(pdb_file, concentration):
    """ Get the number of ions needed for the system. 
    
    Parameters
    ----------
    pdb_file : str
        The path to the pdb file.
    concentration : float
        The concentration of the ions.
    
    Returns
    -------
    int
        The number of NA ions.
    int
        The number of CL ions.
    
    Raises
    ------ 
    ValueError
        If the concentration of ions is negative.
    
    """
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

def CheckComposition(pdb_file):
    """ Checks the composition of the system.
    
    Parameters
    ----------
    pdb_file : str
        The path to the pdb file.
    
    Returns
    -------
    None
    
    Raises
    ------
    ValueError
        If the system is not neutral.

    """
    protein_resnames = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    nucleic_resnames = ["DA", "DC", "DG", "DT", "A", "C", "G", "U", "G5", "C5", "A5", "U5", "G3", "C3", "A3", "U3"]
    water_resnames = ["WAT", "HOH"]
    ion_resnames = ["NA", "CL", "Na+", "Cl-","MG"]
    u = mda.Universe(pdb_file)
    resnames_found = u.residues.resnames
    count = {}
    for resname in resnames_found:
        if resname not in count:
            count[resname] = 1
        else:
            count[resname] += 1
    
    for key in [protein_resnames, nucleic_resnames, water_resnames, ion_resnames]:
        total = 0
        for resname in key:
            if resname in count:
                total += count[resname]
        if total == 0:
            continue
        else:
            print(f"Found {total} residues of type {key}")