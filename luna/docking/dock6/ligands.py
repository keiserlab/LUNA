import os
from pathlib import Path

from luna.util.exceptions import IllegalArgumentError
from luna.util.file import (is_file_valid, get_file_format,
                            create_directory, new_unique_filename,
                            remove_directory)
from luna.wrappers.obabel import convert_molecule
from luna.wrappers.base import MolWrapper
from luna.wrappers.chimera import add_charges

import logging
logger = logging.getLogger()


class PrepareLigand:

    def __init__(self, add_h=False, ph=7.4,
                 use_canonical_smiles=False,
                 add_charges=False, charge_method="am1-bcc",
                 minimize=False, steps=150, forcefield="MMFF94",
                 use_sd=True, min_opts=None):

        self.add_h = add_h
        self.ph = ph

        # If Smiles are provided, transform it to their canonical form.
        self.use_canonical_smiles = use_canonical_smiles

        # Arguments to add charges with Antechamber.
        self.add_charges = add_charges
        self.charge_method = charge_method

        # Arguments to minimize ligand with Open Babel.
        self.minimize = minimize
        self.steps = steps
        self.forcefield = forcefield
        self.use_sd = use_sd
        self.min_opts = min_opts or {}

    def prep(self, ligand_input, output_file, is_smiles=False,
             ligand_name=None):

        if not is_smiles:
            ligand_input = str(Path(ligand_input).absolute())

        output_path = Path(output_file).parent.absolute()
        filename = Path(output_file).stem

        cur_path = os.path.abspath(os.getcwd())

        tmp_path = new_unique_filename(output_path)
        create_directory(tmp_path)

        os.chdir(tmp_path)

        ligand_file = None
        try:
            if is_smiles:
                if not isinstance(ligand_input, str):
                    raise TypeError("Invalid type for the ligand SMILES. "
                                    "A string was expected instead.")

                # Transform Smiles to its canonical form if required.
                if self.use_canonical_smiles:
                    ligand_input = self.to_canonical_form(ligand_input)

                # Add hydrogen if required.
                if self.add_h:
                    ligand_input = self.add_h_to_ligand(ligand_input)

                ligand_file = f"{tmp_path}/{filename}_coords.mol2"
                self.generate_3d_coords(ligand_input, ligand_file)

            elif isinstance(ligand_input, str):
                if not is_file_valid(ligand_input):
                    raise IllegalArgumentError("'%s' does not exist or "
                                               "is not a file." % ligand_input)

                ligand_file = ligand_input
            else:
                raise IllegalArgumentError("'ligand_input' must be a "
                                           "string representing a molecular "
                                           "file or SMILES.")

            if self.add_charges:
                tmp_file = f"{tmp_path}/{filename}_charged.mol2"
                self.add_charges_to_ligand(ligand_file, tmp_file)
                ligand_file = tmp_file

            if self.minimize:
                tmp_file = f"{tmp_path}/{filename}_min.mol2"
                self.minimize_ligand(ligand_file, tmp_file)
                ligand_file = tmp_file

            # Molecular file format.
            ext = get_file_format(ligand_file, 1)

            if ligand_name:
                tmp_file = f"{tmp_path}/{filename}_renamed.{ext}"
                self.rename_ligand_at_file(ligand_name, ligand_file, tmp_file)
                ligand_file = tmp_file

            # Update the charge method in the resulting file.
            if ext == "mol2" and self.add_charges:
                self.set_charge_method(ligand_file)

        except Exception as e:
            ligand_file = None
            logger.exception(e)
            raise

        # Safe finish: set back the working directory to its initial value.
        finally:
            os.chdir(cur_path)
            if ligand_file and os.path.exists(ligand_file):
                os.rename(ligand_file, Path(output_file).absolute())
            remove_directory(tmp_path)

    def to_canonical_form(self, ligand_input):
        return convert_molecule(ligand_input, input_format="smi",
                                output_format="can")

    def add_h_to_ligand(self, ligand_input):
        # If it is a string, but not a file path.
        if not is_file_valid(ligand_input):
            if self.ph is None:
                opts = {"h": ""}
            else:
                opts = {"p": self.ph}

            return convert_molecule(ligand_input, input_format="smi",
                                    output_format="smi", opts=opts)
        else:
            # Not implemented yet.
            pass

    def generate_3d_coords(self, input_file, output_file=None):
        if not isinstance(input_file, str):
            raise IllegalArgumentError("This function accepts only "
                                       "SMILES string.")

        return convert_molecule(input_file, input_format="smi",
                                output_file=output_file, output_format="mol2",
                                opts={"gen3d": ""})

    def rename_ligand_at_file(self, new_name, input_file, output_file=None):
        output_file = output_file or input_file

        output_format = get_file_format(output_file, 1)
        target_row = 0
        if output_format == "mol2":
            target_row = 1

        lines = []
        with open(input_file, "r") as IN:
            c = 0
            for line in IN:
                if c == target_row:
                    lines.append(new_name + "\n")
                else:
                    lines.append(line)

                c += 1

        with open(output_file, "w") as OUT:
            OUT.write("".join(lines))

    def remove_unity_records(self, input_file, output_file=None):
        output_file = output_file or input_file

        lines = []
        with open(input_file, "r") as IN:
            ignore = False
            for line in IN:
                if line.strip() == "@<TRIPOS>UNITY_ATOM_ATTR":
                    ignore = True
                elif line.strip() == "@<TRIPOS>BOND":
                    ignore = False

                if not ignore:
                    lines.append(line)

        with open(output_file, "w") as OUT:
            OUT.write("".join(lines))

    def get_unity_records(self, input_file):
        lines = []
        with open(input_file, "r") as IN:
            ignore = True
            for line in IN:
                if line.strip() == "@<TRIPOS>UNITY_ATOM_ATTR":
                    ignore = False
                elif line.strip() == "@<TRIPOS>BOND":
                    break

                if not ignore:
                    lines.append(line.rstrip("\n"))
        return lines

    def add_charges_to_ligand(self, input_file, output_file):
        input_format = get_file_format(input_file, 1)

        try:
            mol_obj = MolWrapper.from_mol_file(input_file, input_format,
                                               mol_obj_type="rdkit")
            total_charge = mol_obj.get_total_charge()
        except Exception:
            mol_obj = MolWrapper.from_mol_file(input_file, input_format,
                                               mol_obj_type="openbabel")
            total_charge = mol_obj.get_total_charge()

        add_charges(input_file, output_file, total_charge, self.charge_method)

        if input_format == "mol2":
            # Add charge information (Unity records) to the output file
            # so that the molecule can be parsed by Open Babel without
            # raising kekulization issues.
            unity_records = self.get_unity_records(input_file)

            lines = []
            with open(output_file) as IN:
                for line in IN:
                    line = line.rstrip("\n")

                    if line == "@<TRIPOS>BOND":
                        lines.extend(unity_records)

                    lines.append(line)

            with open(output_file, "w") as OUT:
                OUT.write("\n".join(lines))

        self._compare_mols(input_file, output_file)

    def _compare_mols(self, mol_file1, mol_file2):
        mol_format1 = get_file_format(mol_file1, 1)
        mol_format2 = get_file_format(mol_file2, 1)

        # Molecule 1
        mol_obj1 = MolWrapper.from_mol_file(mol_file1, mol_format1,
                                            mol_obj_type="openbabel")
        #
        # Molecule 2
        mol_obj2 = MolWrapper.from_mol_file(mol_file2, mol_format2,
                                            mol_obj_type="openbabel")

        try:
            mol_obj1 = MolWrapper(mol_obj1.as_rdkit())
            mol_obj2 = MolWrapper(mol_obj2.as_rdkit())

            #
            # Maximum Common Substructure.
            mcs = mol_obj1.find_mcs(mol_obj2)

            assert (mcs.numAtoms == mol_obj1.get_num_heavy_atoms()
                    and mcs.numBonds == len(mol_obj1.get_bonds()))

        except AssertionError:
            raise AssertionError("Molecules in the files '%s' and '%s' does "
                                 "not have a structural match."
                                 % (mol_file1, mol_file2))
        except Exception:
            # Compare molecules by InChI identifiers.
            inchi1 = convert_molecule(mol_file1, input_format="mol2",
                                      output_format="inchi")
            inchi2 = convert_molecule(mol_file2, input_format="mol2",
                                      output_format="inchi")

            if inchi1 != inchi2:
                raise AssertionError("Molecules in the files '%s' and '%s' "
                                     "does not have equal InChI identifiers."
                                     % (mol_file1, mol_file2))

    def minimize_ligand(self, input_file, output_file=None):
        opts = {"minimize": "",
                "steps": self.steps,
                "ff": self.forcefield}

        if self.use_sd:
            opts["sd"] = ""

        input_format = get_file_format(input_file, 1)

        return convert_molecule(input_file, input_format=input_format,
                                output_file=output_file, output_format="mol2",
                                opts=opts)

    def set_charge_method(self, input_file, output_file=None):
        output_file = output_file or input_file

        lines = []
        with open(input_file) as IN:
            cm = "GASTEIGER"
            if self.charge_method == "am1-bcc":
                cm = "AM1-BCC"

            c = 0
            for line in IN:
                line = line.rstrip("\n")

                if line == "@<TRIPOS>MOLECULE":
                    mol_started = True

                if c == 4:
                    line = cm

                lines.append(line)

                if mol_started:
                    c += 1

        with open(output_file, "w") as OUT:
            OUT.write("\n".join(lines))
