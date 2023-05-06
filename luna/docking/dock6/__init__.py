import os
from pathlib import Path
from subprocess import Popen, PIPE

from luna.docking.dock6.params import Params
from luna.util.default_values import DOCK6
from luna.util.exceptions import ProcessingFailed


import logging
logger = logging.getLogger()


class Dock6Manager:

    def __init__(self, dock6=DOCK6):

        self.dock6 = dock6

    def _validate_params(self, params):

        if "ligand_atom_file" not in params:
            raise KeyError("The parameter 'ligand_atom_file' has not been "
                           "defined in the input file.")
        elif not os.path.exists(params["ligand_atom_file"]):
            raise FileNotFoundError("The ligand atom file defined in "
                                    "the input file was not found.")
        elif (" " in params["ligand_atom_file"]
                or "\t" in params["ligand_atom_file"]):
            raise ValueError("Spaces or tabs were found in the ligand "
                             "atom file. However, DOCK 6 cannot properly "
                             "deal with them.")

        if "receptor_site_file" not in params:
            raise KeyError("The parameter 'receptor_site_file' has not been "
                           "defined in the input file.")
        elif not os.path.exists(params["receptor_site_file"]):
            raise FileNotFoundError("The receptor site file defined in "
                                    "the input file was not found.")
        elif (" " in params["receptor_site_file"]
                or "\t" in params["receptor_site_file"]):
            raise ValueError("Spaces or tabs were found in the receptor "
                             "site file. However, DOCK 6 cannot properly "
                             "deal with them.")

        if "grid_score_grid_prefix" not in params:
            raise KeyError("The parameter 'grid_score_grid_prefix' has not "
                           "been defined in the input file.")
        elif (" " in params["grid_score_grid_prefix"]
                or "\t" in params["grid_score_grid_prefix"]):
            raise ValueError("Spaces or tabs were found in the grid file "
                             "prefix. However, DOCK 6 cannot properly "
                             "deal with them.")
        else:
            if not os.path.exists(params["grid_score_grid_prefix"] + ".bmp"):
                raise FileNotFoundError("The bmp grid file was not found in "
                                        "the working path.")
            if not os.path.exists(params["grid_score_grid_prefix"] + ".nrg"):
                raise FileNotFoundError("The nrg grid file was not found in "
                                        "the working path.")

    def dock(self, params_file, working_path):

        cur_path = os.path.abspath(os.getcwd())

        working_path = Path(working_path).absolute()
        os.chdir(working_path)

        output_file = f"{working_path}/dock6.out"

        params = Params.load(params_file)
        self._validate_params(params)

        try:
            command = f"{self.dock6} -i '{params_file}' -o '{output_file}'"

            p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
            try:
                stdout, stderr = p.communicate()
            except Exception:
                p.kill()
                raise

            mol_file = params["ligand_outfile_prefix"] + "_scored.mol2"
            if not os.path.exists(mol_file):
                raise ProcessingFailed("The ligand output file was not "
                                       "created. Check the output file "
                                       "'%s' for more information."
                                       % output_file)
            elif os.path.getsize(mol_file) == 0:
                raise ProcessingFailed("The ligand output file is "
                                       "empty. Check the output file "
                                       "'%s' for more information."
                                       % output_file)

        except Exception as e:
            os.chdir(cur_path)
            logger.exception(e)
            raise
