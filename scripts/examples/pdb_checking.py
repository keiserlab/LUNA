from util import logging_ini

from bio.pdb import *

pdbId = "1IL8"
outputPath = "."

# download_pdb(pdbId, outputPath)


control = VersionControl()
print(control.is_pdb_outdated(pdbId))
print(control.is_pdb_obsolete(pdbId))
