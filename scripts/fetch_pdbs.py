from MyBio.PDB.PDBList import PDBList

listToDownload = ["1H2T", "1eve", "1P5E", "3PFP", "1WBG", "1R9O"]

outPath = "../tmp/pharm"

pdbl = PDBList()
for pdbId in listToDownload:
    pdbl.retrieve_pdb_file(pdbId, pdir=outPath, file_format="pdb")
