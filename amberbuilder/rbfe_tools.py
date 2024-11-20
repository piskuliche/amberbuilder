import MDAnalysis as mda

from pathlib import Path


from amberbuilder import mdatools
from amberbuilder.interfaces import Leap

class rbfe_prep:
    def __init__(self, leaprc=[]):
        self.edges = []
        self.leaprc = leaprc
        return
    def add_edge(self, node1, node2):
        self.edges.append([node1, node2])
        return
    
    def prep_edge(self, edge):
        node0, node1 = self._read_edge_complex(edge)
        twoplex = self._build_twostates(node0, node1)
        twoplex = self._reresidue(twoplex)
        self._rename_node_lib_resname(edge[0], "L00")
        self._rename_node_lib_resname(edge[1], "L01")
        mdatools.WritePDB(twoplex, f"twoplex_{edge[0]}_{edge[1]}.pdb")
        self._releap(edge)
        self._clean_edge(edge)
        return 
    
    def prep_edges(self):
        for edge in self.edges:
            try:
                self.prep_edge(edge)
            except Exception as e:
                raise RuntimeError(f"Error in preparing edge {edge[0]}_{edge[1]}") from e
        return
    
    def _read_edge_complex(self, edge):
        node0 = self.amber_universe(edge[0])
        node1 = self.amber_universe(edge[1])
        node0 = self.replace_resname(node0, "LIG", "L00")
        node1 = self.replace_resname(node1, "LIG", "L01")
        return node0, node1
    
    def _build_twostates(self, node0, node1):
        """ Build a two-state complex from two nodes

        Parameters
        ----------
        node0 : mda.Universe
            The first node
        node1 : mda.Universe
            The second node
        
        Returns
        -------
        complex : mda.Universe
            The two-state complex universe

        """
        wat = node0.select_atoms("resname WAT")
        na = node0.select_atoms("resname Na+")
        cl = node0.select_atoms("resname Cl-")
        l00 = node0.select_atoms("resname L00")
        l01 = node1.select_atoms("resname L01")
        other = node0.select_atoms("not (resname WAT or resname Na+ or resname Cl- or resname LOO)")
        return mda.Merge(l00, l01, other, na, cl, wat)
    
    def _rename_node_lib_resname(self, node, newresname):
        """ Copy a node lib file and replace the resname within it.
        
        Parameters
        ----------
        node : str
            The node name
        newresname : str
            The new resname
        
        Returns
        -------
        newlib : str
            The new lib file name
        """
        edge_lib = Path(f"edge_lib/")
        if not edge_lib.exists():
            edge_lib.mkdir()
        node_lib = Path(f"targets/{node}.lib")
        # Copy the file to edge_lib
        newlib = edge_lib / f"{node}_{newresname}.lib"
        with open(node_lib, "r") as F:
            with open(newlib, "w") as G:
                for line in F:
                    if "LIG" in line:
                        line = line.replace("LIG", newresname)
                    G.write(line)
        return


    def _reresidue(self, universe):
        """ Renumber residues in a universe."""
        resnum = 1
        for residue in universe.residues:
            residue.resid = resnum
            resnum += 1
        return universe
    
    def _releap(self, edge):
        leap_lines = []
        for leap_line in self.leaprc:
            leap_lines.append(leap_line)
        leap_lines.append(f"loadamberparams targets/{edge[0]}.frcmod")
        leap_lines.append(f"loadamberparams targets/{edge[1]}.frcmod")
        leap_lines.append(f"loadoff edge_lib/{edge[0]}_L00.lib")
        leap_lines.append(f"loadoff edge_lib/{edge[1]}_L01.lib")
        leap_lines.append(f"complex = loadpdb twoplex_{edge[0]}_{edge[1]}.pdb")
        leap_lines.append("saveamberparm complex unisc.parm7 unisc.rst7")
        leap_lines.append("quit")
        with open("tleap.releap.in", "w") as F:
            for line in leap_lines:
                F.write(line + "\n")
        
        leap = Leap()
        leap.call(f="tleap.releap.in")
        return
    
    
    def _clean_edge(self, edge):
        tleap_dir = Path("tleap/")
        tleap = Path("tleap.releap.in")
        if not tleap_dir.exists():
            tleap_dir.mkdir()
        if tleap.exists():
            tleap.rename(tleap_dir / tleap)

        for file in Path(".").glob(f"twoplex_{edge[0]}_{edge[1]}*"):
            file.unlink()
        
        outputs = Path(f"outputs/{edge[0]}_{edge[1]}")
        if not outputs.exists():
            outputs.mkdir()
        unisc = Path("unisc.parm7")
        unisc_rst7 = Path("unisc.rst7")
        if unisc.exists():
            unisc.rename(outputs / unisc)
        if unisc_rst7.exists():
            unisc_rst7.rename(outputs / unisc_rst7)

        return
    
    @staticmethod
    def amber_universe(shared_name):
        return mda.Universe(f"outputs/{shared_name}.parm7", f"outputs/{shared_name}.rst7", topology_format = 'PARM7', format="INPCRD")
    
    @staticmethod
    def replace_resname(universe, oldresname, newresname):
        """ Replace a resname """
        for residue in universe.residues:
            if residue.resname == oldresname:
                residue.resname = newresname
        return universe
    
