import logging

logger = logging.getLogger()


class VersionControl():

    def __init__(self):
        from MyBio.PDB.PDBList import PDBList
        from urllib.error import URLError, HTTPError

        pdbl = PDBList()
        try:
            new, modified, obsolete = pdbl.get_recent_changes()
        except HTTPError as e:
            logger.exception(e)
            raise HTTPError("PDB server could not be reached.")
        except URLError as e:
            logger.exception(e)
            raise URLError("PDB server could not be reached.")

        self.new = set(new)
        self.modified = set(modified)
        self.obsolete = set(obsolete)

    def is_pdb_outdated(self, pdbId):
        if (pdbId in self.new or pdbId in self.modified):
            return True
        else:
            return False

    def is_pdb_obsolete(self, pdbId):
        if (pdbId in self.obsolete):
            return True
        else:
            return False
