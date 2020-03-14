from util import logging_ini

from mol.obabel import mol_2svg_obabel
from mol.obabel import convert_molecule

path = '../tmp/icode'

infile = "%s/1BGA.7DG.pdb" % (path)
molfile = "%s/1BGA.7DG.mol" % (path)


# CONVERT MOLECULE
opt = {"p": 7}
opt = {}
opt = {"AddPolarH": None}
opt = {"h": None}

convert_molecule(infile, molfile, opt=opt)

# LIGAND FIGURE
output = "%s/1BGA.7DG.svg" % path

opt = {"a": None, "e": None, "P": 500, "d": None}
opt = {"C": None, "e": None, "P": 500, "d": None}

mol_2svg_obabel(molfile, output, opt)
