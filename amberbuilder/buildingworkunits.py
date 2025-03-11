import warnings
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.analysis import align


class BoxBuildABFE:

    def __init__(self, main_input, wu_name:str = "BoxBuildABFE"):
        pass

    def align_complexes(reference_file: Path, mobile_files: list[Path], cwd: Path,
                        select_ref="protein", select_mobile="protein",
                        out_prefix=""):
        # Sort to keep consistency
        mobile_files = sorted(mobile_files)

        Path(cwd).mkdir(exist_ok=True)

        ref_universe = mda.Universe(reference_file)
        ref_selection = ref_universe.select_atoms(select_ref)
        if len(ref_selection) == 0:
            raise ValueError(f"Reference selection '{select_ref}' returned no atoms.")

        # Align and save
        for mobile_file in mobile_files:
            if not mobile_file.exists():
                raise ValueError(f"Mobile file not found: {mobile_file}")

            mobile_universe = mda.Universe(mobile_file)
            mobile_selection = mobile_universe.select_atoms(select_mobile)
            if len(mobile_selection) == 0:
                print(f"Warning: Mobile selection '{select_mobile}' in {mobile_file} returned no atoms. Skipping.")
                continue

            alignment = align.AlignTraj(mobile_universe, ref_universe,
                                        select=select_mobile,
                                        # TODO: watch out with this. We may run into issues when building thousands of systems
                                        in_memory=True)
            alignment.run()

            mobile_universe.atoms.write(cwd / (out_prefix + mobile_file.name))

    def aggregate(self, filepaths: list[Path], biomol="protein") -> mda.Universe:
        print("Reading the target files into a superuniverse.")
        all_systems = []
        self.universes = []

        for cpx in filepaths:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(cpx)
            self.universes.append(u)
            all_systems.append(u.atoms)

            lig_atm_group = u.select_atoms("resname LIG")
            if len(lig_atm_group) == 0:
                err_msg = f"Error: Ligand not found in target file {cpx}"
                raise ValueError(err_msg)
            # lig_atm_group.residues.resnames = ["REM"]
            all_systems.append(lig_atm_group)

        superuniverse = mda.Merge(*all_systems)
        print(
            f"Read {len(self.target_files)} target files into a superuniverse with {superuniverse.atoms.n_atoms} atoms.")
        superuniverse.atoms.write("superuniverse.pdb")
        return superuniverse
