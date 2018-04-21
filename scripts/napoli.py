# Complex validator
from input.entry_validator import PliEntryValidator
from input.complex import DBComplex

from util.default_values import *

from file.util import create_directory

from MyBio.PDB.PDBParser import PDBParser
from MyBio.util import (download_pdb, pdb_object_2block)
from MyBio.selector import ResidueSelector

from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem import MolFromPDBBlock
from rdkit.Chem.rdDepictor import Compute2DCoords

# Get nearby molecules (contacts)
from interaction.contact import get_contacts_for_entity
from interaction.calc_interactions import (calc_all_interactions,
                                           filter_interactions,
                                           InteractionConf)

from mol.groups import find_compound_groups
from mol.fingerprint import generate_fp_for_mols
from mol.chemical_feature import FeatureExtractor
from mol.obabel import convert_molecule

from file.validator import is_file_valid

from util.exceptions import (InvalidNapoliEntry)
from util.config_parser import Config

from database.loader import *
from database.napoli_model import *
from database.helpers import *
from database.util import (default_interaction_filters, prepare_complex_entries)

from data.summary import *
from data.frequency import *

from os.path import exists

import logging
logger = logging.getLogger(__name__)

class NAPOLI_PLI:

    def __init__(self,
                 COMPLEXES=None,
                 PDB_PATH=None,
                 WORKING_PATH=None,
                 OVERWRITE_PATH=False,
                 JOB_CODE=None,
                 DB_CONF_FILE=None,
                 PDB_TEMPLATE=None,
                 ATOM_PROP_FILE=DEFAULT_ATOM_PROP_FILE,
                 INTERACTION_CONF=DEFAULT_INTERACTION_CONF,
                 FIlTERED_INTERACTIONS=True,
                 PH=None,
                 FINGERPRINT_FUNC="pharm2d_fp",
                 SIMILARITY_FUNC="BulkTanimotoSimilarity",
                 RUN_FROM_STEP=None,
                 RUN_UNTIL_STEP=None):

        self.COMPLEXES = COMPLEXES
        self.PDB_PATH = PDB_PATH
        self.WORKING_PATH = WORKING_PATH
        self.OVERWRITE_PATH = OVERWRITE_PATH
        self.JOB_CODE = JOB_CODE
        self.DB_CONF_FILE = DB_CONF_FILE
        self.PDB_TEMPLATE = PDB_TEMPLATE
        self.ATOM_PROP_FILE = ATOM_PROP_FILE
        self.INTERACTION_CONF = INTERACTION_CONF
        self.FIlTERED_INTERACTIONS = FIlTERED_INTERACTIONS
        self.PH = PH
        self.ADD_HYDROG = True
        self.FINGERPRINT_FUNC = FINGERPRINT_FUNC
        self.SIMILARITY_FUNC = SIMILARITY_FUNC
        self.RUN_FROM_STEP = RUN_FROM_STEP
        self.RUN_UNTIL_STEP = RUN_UNTIL_STEP

    def run(self):

        try:

            if self.JOB_CODE and not self.DB_CONF_FILE:
                raise IllegalArgumentError("You informed a project id, but "
                                           "none DB configuration file was "
                                           "informed. A project id should be "
                                           "defined only when it is necessary "
                                           "to save results at a database.")

            elif not self.JOB_CODE and not self.WORKING_PATH:
                raise IllegalArgumentError("None variable used to determine "
                                           "a working path was informed.")

            ##########################################################
            step = "Prepare project directory"
            if self.WORKING_PATH:
                workingPath = self.WORKING_PATH
            else:
                workingPath = "%s/projects/%s" % (NAPOLI_PATH, self.JOB_CODE)

            workingPdbPath = "%s/pdbs" % workingPath
            create_directory(workingPath, self.OVERWRITE_PATH)
            create_directory(workingPdbPath)

            logger.info("Results and temporary files will be saved at '%s'."
                        % workingPath)

            # Set the LOG file at the working path
            logfile = "%s/napoli.log" % workingPath
            filehandler = logging.FileHandler(logfile, 'a')
            formatter = logging.Formatter('%(asctime)s - %(name)s - '
                                          '%(levelname)s - %(filename)s - '
                                          '%(funcName)s - %(lineno)d => '
                                          '%(message)s')
            filehandler.setFormatter(formatter)
            # Remove the existing file handlers
            for hdlr in logger.handlers[:]:
                if isinstance(hdlr, logger.FileHander):
                    logger.removeHandler(hdlr)
            logger.addHandler(filehandler)      # set the new handler
            # Set the log level to INFO, DEBUG as the default is ERROR
            logger.setLevel("INFO")

            logger.info("Logs will be saved in '%s'." % logfile)

            if self.JOB_CODE:
                logger.info("The project id '%s' was informed. "
                            "nAPOLI will try to save all results "
                            "in the DB using such id." % self.JOB_CODE)

            # Define which PDB path to use.
            if self.PDB_PATH is None:
                # If none PDB path was defined and none DB conf was defined
                # it is better to create a new directory to avoid
                # overwriting Public directories
                if self.DB_CONF_FILE is None:
                    pdbPath = workingPdbPath
                else:
                    pdbPath = DEFAULT_PDB_PATH
            else:
                pdbPath = self.PDB_PATH

            logger.info("PDB files will be obtained from and/or downloaded "
                        "to '%s'." % pdbPath)


            ##########################################################
            if self.DB_CONF_FILE:
                step = "Preparing database"
                logger.info("A DB configuration file was defined. "
                            "Starting a new DB connection...")

                config = Config(self.DB_CONF_FILE)
                dbConf = config.get_section_map("database")
                db = DBLoader(**dbConf)

                dbInter = InteractionManager(db)

                # TODO: create all mappers at once
                db.new_mapper(CompTypeCount, "comp_type_count")
                db.new_mapper(InterTypeCount, "inter_type_count")
                db.new_mapper(Project, "project")
                db.new_mapper(Complex, "complex")

                # TODO: Get the project ID in the DB.
                projectRow = db.session.query(Project).filter(Project.job_code == self.JOB_CODE).one()
                projectId = projectRow.id

                interTypeRows = db.session.query(InterType).all()
                interIdByType = {r.type: r.id for r in interTypeRows}
                interFilters = default_interaction_filters(interIdByType, self.INTERACTION_CONF)

                compTypeRows = db.session.query(CompoundType).all()
                compIdByType = {r.type: r.id for r in compTypeRows}

                if self.COMPLEXES is None:
                    logger.info("No complex was defined. "
                                "It will try to recover from the database.")

                    dbComplexes = ComplexManager(db).get_complexes(projectId)
                    self.COMPLEXES = prepare_complex_entries(dbComplexes)

            logger.info("Number of complexes to be processed: %d." % len(self.COMPLEXES))

            if self.ADD_HYDROG:
                if self.PH:
                    logger.info("Hydrogen atoms will be added to the "
                                "structures according to the pH %s." % self.PH)
                else:
                    logger.info("A pH value was not defined. "
                                "Hydrogen atoms will be added to the "
                                "structures without considering a pH "
                                "value.")
            else:
                logger.info("Option ADD_HYDROG disabled. nAPOLI will not add "
                            "hydrogen atoms to the structures.")

            ##########################################################
            # Set some variables
            validator = PliEntryValidator()
            pdbParser = PDBParser(PERMISSIVE=True,
                                  QUIET=True,
                                  FIX_ATOM_NAME_CONFLICT=True,
                                  FIX_OBABEL_FLAGS=True)

            # featFactory = ChemicalFeatures.BuildFeatureFactory(self.ATOM_PROP_FILE)
            # sigFactory = SigFactory(featFactory, minPointCount=2,
            #                         maxPointCount=3,
            #                         trianglePruneBins=False)
            # sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
            # sigFactory.Init()

            # featExtractor = FeatureExtractor(featFactory)

            # boundaryConf = InteractionConf({"boundary_cutoff": 7})

            fingerprints = []

            # The calculus of interactions depend on the pharmacophore
            # definition (ATOM_PROP_FILE) and the pH (PH).
            # Thus, if the user has kept such values unchanged and has
            # not uploaded any PDB file, nAPOLI can try to select the
            # pre-computed interactions from a database.
            getInterFromDB = (self.PDB_PATH is None and
                              self.ATOM_PROP_FILE == DEFAULT_ATOM_PROP_FILE and
                              self.PH is None and
                              self.DB_CONF_FILE is not None)
            print(getInterFromDB)
            exit()
            ##########################################################

            for pliComplex in self.COMPLEXES:

                try:
                    entry = pliComplex.to_string(ENTRIES_SEPARATOR)
                    myBioLigand = ("H_%s" % pliComplex.ligName,
                                   pliComplex.ligNumber, " ")

                    if self.DB_CONF_FILE:
                        step = "Check complex existance"
                        if isinstance(pliComplex, DBComplex) is False:
                            raise InvalidNapoliEntry("Entry '%s' is not a DBComplex "
                                                     "object. So, nAPOLI cannot "
                                                     "obtain an id information "
                                                     "for this complex. No database "
                                                     "update is possible." % entry)
                        elif pliComplex.id is None:
                            raise InvalidNapoliEntry("An invalid id for the entry '%s'"
                                                     "was defined." % entry)
                        else:
                            #TODO: Test if this ID exists in the DB.
                            dbComplex = (db.session
                                         .query(Complex)
                                         .filter(Complex.id == pliComplex.id)
                                         .first())

                            if dbComplex is None:
                                raise InvalidNapoliEntry("Entry '%s' with complex id equal to '%d' "
                                                         "does not exist in the database."
                                                         % (entry, pliComplex.id))

                    ##########################################################
                    step = "Validate nAPOLI entry"
                    if not validator.is_valid(entry):
                        raise InvalidNapoliEntry("Entry '%s' does not match "
                                                 "the nAPOLI entry format."
                                                 % entry)

                    ##########################################################
                    #TODO: option to change name of default PDB name.
                    #      BioPython cria arquivo com .env
                    step = "Check PDB file existance"
                    # User has defined a specific directory.
                    if self.PDB_PATH:
                        pdbFile = "%s/%s.pdb" % (pdbPath, pliComplex.pdb)
                        if is_file_valid(pdbFile) is False:
                            raise FileNotFoundError("The PDB file '%s' was "
                                                    "not found." % pdbFile)
                    # If none DB conf. was defined try to download the PDB.
                    elif self.DB_CONF_FILE is None:
                        pdbFile = "%s/pdb%s.ent" % (pdbPath,
                                                    pliComplex.pdb.lower())
                        download_pdb(pliComplex.pdb, pdbPath)
                    else:
                        # TODO: check if PDB is updated on the DB
                        pdbFile = "%s/pdb%s.ent" % (pdbPath,
                                                    pliComplex.pdb.lower())
                        pass

                    #####################################
                    step = "Prepare protein structure"
                    if self.ADD_HYDROG:
                        prepPdbFile = "%s/%s.H.pdb" % (workingPdbPath,
                                                       pliComplex.pdb)

                        if not exists(prepPdbFile):
                            if self.PH:
                                obOpt = {"p": 7, "error-level": 5}
                            else:
                                obOpt = {"h": None, "error-level": 5}

                            convert_molecule(pdbFile, prepPdbFile, opt=obOpt)

                        pdbFile = prepPdbFile

                    ##########################################################
                    step = "Parse PDB file"
                    structure = pdbParser.get_structure(pliComplex.pdb,
                                                        pdbFile)
                    ligand = structure[0][pliComplex.chain][myBioLigand]

                    ##########################################################
                    step = "Generate structural fingerprint"
                    ligSel = ResidueSelector({ligand})
                    ligBlock = pdb_object_2block(structure, ligSel,
                                                 write_conects=False)
                    rdLig = MolFromPDBBlock(ligBlock)
                    rdLig.SetProp("_Name", entry)
                    fpOpt = {"sigFactory": sigFactory}
                    fp = generate_fp_for_mols([rdLig], self.FINGERPRINT_FUNC,
                                              fpOpt=fpOpt, critical=True)[0]
                    fingerprints.append(fp)

                    ##########################################################
                    step = "Calculate interactions"
                    nearbyResidues = get_contacts_for_entity(structure[0],
                                                             ligand,
                                                             level='R')

                    compounds = set([x[1] for x in nearbyResidues])
                    groupsByCompounds = {}
                    for comp in compounds:
                        groups = find_compound_groups(comp, featExtractor)
                        groupsByCompounds[comp] = groups

                    targetGroups = [groupsByCompounds[x]
                                    for x in groupsByCompounds
                                    if x.get_id()[0] != " "]
                    nearbyGroups = [groupsByCompounds[x]
                                    for x in groupsByCompounds
                                    if x != ligand]

                    #TODO: Selecionar da base de dados se o usuario nao submeteu PDBs.

                    if self.FIlTERED_INTERACTIONS:
                        # First it calculates all interactions using a
                        # boundary limit. When using a database, it is
                        # useful to save all potential interactions.
                        # So, it is possible to filter interactions
                        # faster than recalculate them for each
                        # modification in the interaction criteria.
                        allInter = calc_all_interactions(targetGroups,
                                                         nearbyGroups,
                                                         conf=boundaryConf)

                        # Then it applies a filtering function.
                        filteredInter = filter_interactions(allInter,
                                                            conf=self.INTERACTION_CONF)
                        dbInteractions = allInter
                    else:
                        # It filters the interactions by using a cutoff at once.
                        filteredInter = calc_all_interactions(targetGroups,
                                                              nearbyGroups,
                                                              conf=self.INTERACTION_CONF)
                        dbInteractions = filteredInter

                    ##########################################################
                    if self.DB_CONF_FILE:
                        #TODO: Remover comentario
                        #TODO: Vincular interação a um PDB
                        # dbInter.insert_interactions(dbInteractions, dbComplex)
                        pass
                    else:
                        #TODO: Criar arquivo de saida formato .csv
                        pass

                    ##########################################################
                    step = "Calculate atom type statistics"
                    ligandGroups = groupsByCompounds[ligand]
                    groupTypes = count_group_types(ligandGroups, compIdByType)

                    if self.DB_CONF_FILE:
                        for typeId in groupTypes:
                            db.session.add(CompTypeCount(pliComplex.id, projectId,
                                                         typeId, groupTypes[typeId]))
                        #TODO: Remover comentario
                        # db.approve_session()

                    ##########################################################
                    step = "Calculate interaction type statistics"
                    interactionTypes = count_interaction_types(filteredInter, interIdByType)

                    if self.DB_CONF_FILE:
                        for typeId in interactionTypes:
                            db.session.add(InterTypeCount(pliComplex.id, projectId,
                                                          typeId, interactionTypes[typeId]))

                        #TODO: Remover comentario
                        # db.approve_session()

                    print("Last step: '%s'" % step)

                except Exception as e:
                    logger.exception(e)

                    # TODO: Inform that the Entry failed

            ##########################################################
            # TODO: after clusterize molecules to align the ligands in a same way.

            # step = "Generate 2D ligand structure"
            # # pdbExtractor = Extractor(structure[0][pliComplex.chain])
            # # ligandPdbFile = "%s/%s.pdb" % (workingPdbPath, entry)
            # # pdbExtractor.extract_residues([ligandBio], ligandPdbFile)
            # ligSel = ResidueSelector({ligand})
            # ligBlock = pdb_object_2block(structure, ligSel,
            #                              write_conects=False)
            # rdLig = MolFromPDBBlock(ligBlock)
            # Compute2DCoords(rdLig)

        except Exception as e:
            logger.exception(e)
