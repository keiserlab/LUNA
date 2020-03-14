from util import logging_ini
from mol.clustering import cluster_fps_butina
from rdkit import Chem
from rdkit.Chem.AllChem import *
from rdkit.Chem import Draw
from mol.fingerprint import generate_fp_for_mols
from os.path import exists
from os import makedirs
from file.util import *
from mol.obabel import mol_2svg_obabel
from rdkit.Chem import rdFMCS
import glob


path = "../tmp/clustering"

molFiles = glob.glob('%s/ligands_mol/*.mol' % path)

mols = [Chem.MolFromMolFile(m) for m in molFiles]

# mols = mols[1:20]
mols = mols[1:10]

func = "pharm2d_fp"
# func = "morgan_fp"
# func = "topological_fp"
# func = "atom_pairs_fp"
opt = {"nBits": 12, "features": False}

fps = generate_fp_for_mols(mols, func, opt)
fpsOnly = [x["fp"] for x in fps]

clusters = cluster_fps_butina(fpsOnly, 0.07)


clustersPath = "%s/clusters" % path

clusterDict = {}
opt = {"C": None, "e": None, "P": 500, "d": None}
for clusterIdx, cluster in enumerate(clusters):
    clusterPath = "%s/%d" % (clustersPath, clusterIdx)

    if (exists(clusterPath)):
        clear_directory(clusterPath)

    makedirs(clusterPath)

    templateMol = fps[cluster[0]]["mol"]
    templateEntry = get_filename(templateMol)
    templateMol = "%s/ligands_mol/%s.mol" % (path, templateEntry)

    tempObj = Chem.MolFromMolFile(templateMol)
    Chem.rdDepictor.Compute2DCoords(tempObj)
    print(templateEntry)

    for molIdx in cluster:
        mol = fps[molIdx]["mol"]
        entry = get_filename(mol)
        mol = "%s/ligands_mol/%s.mol" % (path, entry)
        clusterDict[entry] = clusterIdx
        svg = "%s/%s.svg" % (clusterPath, entry)

        molObj = Chem.MolFromMolFile(mol)
        Chem.rdDepictor.Compute2DCoords(molObj)
        align = Chem.rdMolAlign.GetO3A(molObj, tempObj)
        TransformMol(molObj, align.Trans()[1])
        # Chem.rdDepictor.Compute2DCoords(molObj)

        # GenerateDepictionMatching3DStructure(molObj, molObj)

        # res = rdFMCS.FindMCS(mols)
        # resObj = Chem.MolFromSmarts(res.smartsString)
        # highlightAtoms = molObj.GetSubstructMatch(resObj)

        # GenerateDepictionMatching2DStructure(molObj, tempObj)
        # Draw.MolToFile(molObj, svg, highlightAtoms=highlightAtoms)
        Draw.MolToFile(molObj, svg)
