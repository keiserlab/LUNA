from util.entry import (DBLigandEntry, PLIEntryValidator)

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
                                           apply_interaction_criteria,
                                           InteractionConf)

from mol.groups import find_compound_groups
from mol.fingerprint import generate_fp_for_mols
from mol.chemical_feature import FeatureExtractor
from mol.obabel import convert_molecule

from file.validator import is_file_valid

from util.exceptions import (InvalidNapoliEntry, DuplicateEntry)
from util.config_parser import Config

from database.loader import *
from database.napoli_model import *
from database.helpers import *
from database.util import (get_ligand_tbl_join_filter,
                           default_interaction_filters,
                           format_db_ligand_entries,
                           format_db_interactions)

from data.summary import *

from os.path import exists

import logging
logger = logging.getLogger(__name__)


class InteractionsProject:

    def __init__(self,
                 entries,
                 pdb_path,
                 working_path):

        pass


    def perceive_chemical_groups():
        pass

    def perceive_interactions():
        pass


# POPULATE RCSB.
class RCSB_PLI_Population(InteractionsProject):
    def __init__(self,
                 pdb_path=None):
        pass


class DB_PLI_Project(InteractionsProject):

    def __init__(self,
                 job_code,
                 db_conf):
        pass

class Local_PLI_Project(InteractionsProject):

    def __init__(self,
                 working_path,
                 db_conf=None):
        pass

class NAPOLI_PLI:

    def __init__(self,
                entries=None,
                pdb_path=None,
                 working_path=None,
                 overwrite_path=False,
                 job_code=None,
                 db_conf_file=None,
                 populate_rcsb_tables=False,
                 pdb_template=None,
                 atom_prop_file=DEFAULT_ATOM_PROP_FILE,
                 interaction_conf=DEFAULT_INTERACTION_CONF,
                 force_calc_interactions=False,
                 save_all_interactions=True,
                 ph=None,
                 fingerprint_func="pharm2d_fp",
                 similarity_func="BulkTanimotoSimilarity",
                 run_from_step=None,
                 run_until_step=None):

        self.entries = entries

        self.pdb_path = pdb_path
        self.working_path = working_path
        self.overwrite_path = overwrite_path
        self.job_code = job_code
        self.db_conf_file = db_conf_file

        self.populate_rcsb_tables = populate_rcsb_tables

        self.pdb_template = pdb_template
        self.atom_prop_file = atom_prop_file
        self.interaction_conf = interaction_conf
        self.save_all_interactions = save_all_interactions
        self.ph = ph
        self.add_hydrog = True
        self.fingerprint_func = fingerprint_func
        self.similarity_func = similarity_func
        self.run_from_step = run_from_step
        self.run_until_step = run_until_step

    def run(self):
        # TODO: Create functions for reduce the amount of code in the run()
        try:
            if self.job_code and self.populate_rcsb_tables:
                logger.info("You informed a job code and set RCSB's population mode to true. "
                            "However, a job code has a higher prority over the latter. "
                            "So, the RCSB tables will not be populated.")
                self.populate_rcsb_tables = False

            if self.job_code and not self.db_conf_file:
                raise IllegalArgumentError("You informed a job code, but "
                                           "none database configuration file "
                                           "was informed. So, it would not "
                                           "be possible to store interactions "
                                           "and related information into the "
                                           "database.")
            elif self.populate_rcsb_tables and not self.db_conf_file:
                raise IllegalArgumentError("You activate the RCSB's "
                                           "population mode, but none "
                                           "database configuration file "
                                           "was informed. So, it would not "
                                           "be possible to store interactions "
                                           "and related information into the "
                                           "database.")
            elif not self.job_code and not self.working_path:
                raise IllegalArgumentError("Neither a job code or a working "
                                           "path was defined."
                                           "Without them, it is not possible "
                                           "to decide which output directory "
                                           "to use.")

            if self.populate_rcsb_tables:
                logger.info("RCSB's population mode active. "
                            "Interactions and related information will "
                            "be stored into the defined database without "
                            "relating them to a job code.")

            ##########################################################
            step = "Prepare project directory"
            if self.working_path:
                working_path = self.working_path
            else:
                working_path = "%s/projects/%s" % (DEFAULT_NAPOLI_PATH,
                                                   self.job_code)

            working_pdb_path = "%s/pdbs" % working_path
            create_directory(working_path, self.overwrite_path)
            create_directory(working_pdb_path)

            logger.info("Results and temporary files will be saved at '%s'."
                        % working_path)

            # Set the LOG file at the working path
            # logfile = "%s/napoli.log" % working_path
            # filehandler = logging.FileHandler(logfile, 'a')
            # formatter = logging.Formatter('%(asctime)s - %(name)s - '
            #                               '%(levelname)s - %(filename)s - '
            #                               '%(funcName)s - %(lineno)d => '
            #                               '%(message)s')
            # filehandler.setFormatter(formatter)
            # # Remove the existing file handlers
            # for hdlr in logger.handlers[:]:
            #     if isinstance(hdlr, logger.FileHander):
            #         logger.removeHandler(hdlr)
            # logger.addHandler(filehandler)      # set the new handler
            # # Set the log level to INFO, DEBUG as the default is ERROR
            # logger.setLevel("INFO")

            # logger.info("Logs will be saved in '%s'." % logfile)

            # TODO: avisar se eh csv ou salvar no BD.
            # OUTPUT MODE:
            if self.job_code:
                logger.info("A job code ('%s') was informed. "
                            "nAPOLI will try to save all results "
                            "in the database using such id." % self.job_code)

            # Define which PDB path to use.
            if self.pdb_path is None:
                # If none PDB path was defined and none DB conf was defined
                # it is better to create a new directory to avoid
                # overwriting Public directories
                if self.db_conf_file is None:
                    pdb_path = working_pdb_path
                else:
                    pdb_path = DEFAULT_PDB_PATH
            else:
                pdb_path = self.pdb_path

            logger.info("PDB files will be obtained from and/or downloaded "
                        "to '%s'." % pdb_path)

            ##########################################################
            if self.db_conf_file:
                step = "Preparing database"
                logger.info("A database configuration file was defined. "
                            "Starting a new database connection...")

                config = Config(self.db_conf_file)
                dbConf = config.get_section_map("database")
                db = DBLoader(**dbConf)

                rcsb_inter_manager = RCSBInteractionManager(db)
                proj_inter_manager = ProjectInteractionManager(db)

                # TODO: create all mappers at once
                db.new_mapper(CompTypeCount, "comp_type_count")
                db.new_mapper(InterTypeCount, "inter_type_count")
                db.new_mapper(Project, "project")
                db.new_mapper(LigandEntry, "ligand_entry")
                db.new_mapper(Status, "status")

                if self.job_code:
                    projectRow = (db.session
                                  .query(Project)
                                  .filter(Project.job_code ==
                                          self.job_code).one())
                    projectId = projectRow.id

                interTypeRows = db.session.query(InterType).all()
                interIdByType = {r.type: r.id for r in interTypeRows}
                interFilters = default_interaction_filters(interIdByType,
                                                           self.interaction_conf)

                compTypeRows = db.session.query(CompoundType).all()
                compIdByType = {r.type: r.id for r in compTypeRows}

                status_rows = db.session.query(Status).all()
                status_by_id = {r.id: r.title for r in status_rows}

                if self.entries is None:
                    if self.job_code is None:
                        raise IllegalArgumentError("No ligand entry list was "
                                                   "defined. "
                                                   "The program could try to "
                                                   "recover a list from the "
                                                   "database, but no job code "
                                                   "was defined either.")
                    else:
                        logger.info("No ligand entry was defined. "
                                    "It will try to recover a ligand entry "
                                    "list from the database.")
                        db_ligand_entries = (LigandEntryManager(db)
                                             .get_entries(projectId))
                        self.entries = format_db_ligand_entries(db_ligand_entries)

            logger.info("Number of ligand entries to be processed: %d." %
                        len(self.entries))

            if self.add_hydrog:
                if self.ph:
                    logger.info("Hydrogen atoms will be added to the "
                                "structures according to the pH %s." % self.ph)
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
            validator = PLIEntryValidator()
            pdb_parser = PDBParser(PERMISSIVE=True,
                                   QUIET=True,
                                   FIX_ATOM_NAME_CONFLICT=True,
                                   FIX_OBABEL_FLAGS=True)

            feat_factory = (ChemicalFeatures
                            .BuildFeatureFactory(self.atom_prop_file))

            sigFactory = SigFactory(feat_factory, minPointCount=2,
                                    maxPointCount=3,
                                    trianglePruneBins=False)
            sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
            sigFactory.Init()

            feat_extractor = FeatureExtractor(feat_factory)

            boundary_conf = InteractionConf({"boundary_cutoff": 7})

            fingerprints = []

            # The calculus of interactions depend on the pharmacophore
            # definition (ATOM_PROP_FILE) and the pH (PH).
            # Thus, if the user has kept such values unchanged and has
            # not uploaded any PDB file, nAPOLI can try to select the
            # pre-computed interactions from a database.
            is_inter_params_default = (self.pdb_path is None and
                                       self.atom_prop_file == DEFAULT_ATOM_PROP_FILE and
                                       self.ph is None and
                                       self.db_conf_file is not None)

            if is_inter_params_default:
                logger.info("Default configuration was kept unchanged. "
                            "nAPOLI will try to select pre-computed "
                            "interactions from the defined database")

            ##########################################################

            # Loop over each entry.
            for ligand_entry in self.entries:
                try:
                    entry_str = ligand_entry.to_string(ENTRIES_SEPARATOR)
                    myBioLigand = ("H_%s" % ligand_entry.lig_name,
                                   ligand_entry.lig_num,
                                   ligand_entry.lig_icode)

                    # TODO: Verificar se ligante existe no modelo (BioPDB)

                    db_ligand_entity = None
                    if self.db_conf_file and self.job_code:
                        step = "Check ligand entry existance"
                        if isinstance(ligand_entry, DBLigandEntry) is False:
                            message = ("Entry '%s' is not a DBLigandEntry "
                                       "object. So, nAPOLI cannot obtain an "
                                       "id information for this ligand entry "
                                       "and no database update will be "
                                       "possible." % entry_str)
                            raise InvalidNapoliEntry(message)
                        elif ligand_entry.id is None:
                            message = ("An invalid id for the entry '%s'"
                                       "was defined." % entry_str)
                            raise InvalidNapoliEntry(message)
                        else:
                            db_ligand_entity = (db.session
                                                .query(LigandEntry)
                                                .filter(LigandEntry.id == ligand_entry.id)
                                                .first())
                            if db_ligand_entity is None:
                                message = ("Entry '%s' with id equal "
                                           "to '%d' does not exist in the "
                                           "database."
                                           % (entry_str, ligand_entry.id))
                                raise InvalidNapoliEntry(message)
                    elif self.populate_rcsb_tables:
                        step = "Check ligand existance in RCSB tables"
                        db_ligand_entity = (db.session
                                            .query(Ligand)
                                            .filter(Ligand.pdb_id == ligand_entry.pdb_id and
                                                    Ligand.chain_id == ligand_entry.chain_id and
                                                    Ligand.lig_name == ligand_entry.lig_name and
                                                    Ligand.lig_num == ligand_entry.lig_num and
                                                    Ligand.lig_icode == ligand_entry.lig_icode)
                                            .first())

                        # If this entry does not exist in the database,
                        # raise an error.
                        if db_ligand_entity is None:
                            message = ("Entry '%s' does not exist in the "
                                       "table 'ligand'." % entry_str)
                            raise InvalidNapoliEntry(message)
                        # If there are already interactions to this entry,
                        # raise an error.
                        elif status_by_id[db_ligand_entity.status_id] == "AVAILABLE":
                            raise DuplicateEntry("Interactions to the entry "
                                                 "'%s' already exists in "
                                                 "the database."
                                                 % entry_str)

                    ##########################################################
                    step = "Validate nAPOLI entry"
                    if not validator.is_valid(entry_str):
                        raise InvalidNapoliEntry("Entry '%s' does not match "
                                                 "the nAPOLI entry format."
                                                 % entry_str)

                    ##########################################################
                    # TODO: option to change name of default PDB name.
                    #      BioPython cria arquivo com .env
                    step = "Check PDB file existance"
                    # User has defined a specific directory.
                    if self.pdb_path:
                        pdbFile = "%s/%s.pdb" % (pdb_path, ligand_entry.pdb_id)
                        # If the file does not exist or is invalid.
                        if is_file_valid(pdbFile) is False:
                            raise FileNotFoundError("The PDB file '%s' was "
                                                    "not found." % pdbFile)
                    # If none DB conf. was defined try to download the PDB.
                    # Here, pdb_path is equal to the working_pdb_path
                    elif self.db_conf_file is None:
                        pdbFile = "%s/pdb%s.ent" % (working_pdb_path,
                                                    ligand_entry.pdb_id.lower())
                        download_pdb(ligand_entry.pdb_id, working_pdb_path)
                    else:
                        pdbFile = "%s/pdb%s.ent" % (pdb_path,
                                                    ligand_entry.pdb_id.lower())

                        # If the file does not exist or is invalid, it is
                        # better to download a new PDB file at the working path.
                        if not is_file_valid(pdbFile):
                            pdbFile = "%s/pdb%s.ent" % (working_pdb_path,
                                                        ligand_entry.pdb_id.lower())
                            download_pdb(ligand_entry.pdb_id, working_pdb_path)
                        else:
                            # TODO: check if PDB is updated on the DB
                            pass

                    #####################################
                    step = "Prepare protein structure"
                    # TODO: Distinguish between files that need to remove hydrogens (NMR, for example).
                    # TODO: Add parameter to force remove hydrogens (user marks it)
                    if self.add_hydrog:
                        prepPdbFile = "%s/%s.H.pdb" % (working_pdb_path,
                                                       ligand_entry.pdb_id)

                        if not exists(prepPdbFile):
                            if self.ph:
                                obOpt = {"p": 7, "error-level": 5}
                            else:
                                obOpt = {"h": None, "error-level": 5}

                            convert_molecule(pdbFile, prepPdbFile, opt=obOpt)

                        pdbFile = prepPdbFile

                    ##########################################################
                    step = "Parse PDB file"
                    structure = pdb_parser.get_structure(ligand_entry.pdb_id,
                                                         pdbFile)
                    ligand = structure[0][ligand_entry.chain_id][myBioLigand]

                    ##########################################################
                    step = "Generate structural fingerprint"
                    lig_sel = ResidueSelector({ligand})
                    lig_block = pdb_object_2block(structure, lig_sel,
                                                  write_conects=False)
                    rdLig = MolFromPDBBlock(lig_block)
                    rdLig.SetProp("_Name", entry_str)
                    fpOpt = {"sigFactory": sigFactory}
                    fp = generate_fp_for_mols([rdLig], self.fingerprint_func,
                                              fpOpt=fpOpt, critical=True)[0]
                    fingerprints.append(fp)

                    ##########################################################

                    # TODO: Preparar estrutura somente se não for selecionar
                    # da base de dados

                    step = "Calculate interactions"

                    calculate_interactions = True
                    if (is_inter_params_default and
                            not self.populate_rcsb_tables):
                        logger.info("Trying to select pre-computed interactions for "
                                    "the entry '%s'." % entry_str)

                        db_ligand_entity = (db.session
                                            .query(Ligand)
                                            .filter(Ligand.pdb_id == ligand_entry.pdb_id and
                                                    Ligand.chain_id == ligand_entry.chain_id and
                                                    Ligand.lig_name == ligand_entry.lig_name and
                                                    Ligand.lig_num == ligand_entry.lig_num and
                                                    Ligand.lig_icode == ligand_entry.lig_icode)
                                            .first())

                        if db_ligand_entity:
                            if status_by_id[db_ligand_entity.status_id] == "AVAILABLE":
                                join_filter = get_ligand_tbl_join_filter(ligand_entry, Ligand)
                                db_interactions = (rcsb_inter_manager
                                                   .select_interactions(join_filter,
                                                                        interFilters))
                                filtered_inter = format_db_interactions(structure, db_interactions)

                                logger.info("%d pre-computed interaction(s) found in the "
                                               "database for the entry '%s'."
                                               % (len(filtered_inter), entry_str))

                                calculate_interactions = False
                            else:
                                logger.info("The entry '%s' exists in the database, but "
                                            "there is no pre-computed interaction available. "
                                            "So, nAPOLI will calculate the interactions to this "
                                            "ligand." % entry_str)
                        else:
                            logger.info("The entry '%s' does not exist in the "
                                        "database." % entry_str)

                    if not self.populate_rcsb_tables:
                        nb_compounds = get_contacts_for_entity(structure[0],
                                                               ligand,
                                                               level='R')

                        compounds = set([x[1] for x in nb_compounds])
                        grps_by_compounds = {}
                        for comp in compounds:
                            # TODO: verificar se estou usando o ICODE nos residuos
                            groups = find_compound_groups(comp, feat_extractor)
                            grps_by_compounds[comp] = groups

                        trgt_grps = [grps_by_compounds[x]
                                     for x in grps_by_compounds
                                     if x.get_id()[0] != " "]
                        nb_grps = [grps_by_compounds[x]
                                   for x in grps_by_compounds
                                   if x != ligand]

                    if calculate_interactions:
                        if self.save_all_interactions:
                            # First it calculates all interactions using a
                            # boundary limit. When using a database, it is
                            # useful to save all potential interactions.
                            # So, it is possible to filter interactions
                            # faster than recalculate them for each
                            # modification in the interaction criteria.
                            all_inter = calc_all_interactions(trgt_grps,
                                                              nb_grps,
                                                              conf=boundary_conf)

                            # Then it applies a filtering function.
                            filtered_inter = apply_interaction_criteria(all_inter,
                                                                        conf=self.interaction_conf)

                            db_interactions = all_inter
                        else:
                            # It filters the interactions by using a cutoff at once.
                            filtered_inter = calc_all_interactions(trgt_grps,
                                                                   nb_grps,
                                                                   conf=self.interaction_conf)
                            db_interactions = filtered_inter

                        ##########################################################
                        if self.db_conf_file and self.job_code:
                            #TODO: Remover comentario
                            #TODO: Adicionar ICODE nas chaves de entrada que o usuario submete
                            # proj_inter_manager.insert_interactions(db_interactions, db_ligand_entity)
                            pass
                        elif self.populate_rcsb_tables:
                            # rcsb_inter_manager.insert_interactions(db_interactions, db_ligand_entity)
                            pass
                        else:
                            #TODO: Criar arquivo de saida formato .tsv
                            pass

                    if not self.populate_rcsb_tables:
                        ##########################################################
                        step = "Calculate atom type statistics"
                        ligand_grps = grps_by_compounds[ligand]
                        grp_types = count_group_types(ligand_grps, compIdByType)
                        if self.db_conf_file:
                            for typeId in grp_types:
                                db.session.add(CompTypeCount(ligand_entry.id, projectId,
                                                             typeId, grp_types[typeId]))
                            #TODO: Remover comentario
                            db.approve_session()

                        ##########################################################
                        step = "Calculate interaction type statistics"
                        summary_trgts = {ligand}
                        summary_trgts = None
                        interaction_types = count_interaction_types(filtered_inter,
                                                                    summary_trgts,
                                                                    interIdByType)
                        if self.db_conf_file:
                            for typeId in interaction_types:
                                db.session.add(InterTypeCount(ligand_entry.id, projectId,
                                                              typeId, interaction_types[typeId]))

                            #TODO: Remover comentario
                            db.approve_session()

                        print("Last step: '%s'" % step)

                except Exception as e:
                    logger.exception(e)

                    # TODO: Inform that the Entry failed

            ##########################################################
            # TODO: after clusterize molecules to align the ligands in a same way.

            # step = "Generate 2D ligand structure"
            # # pdbExtractor = Extractor(structure[0][ligand_entry.chain])
            # # ligandPdbFile = "%s/%s.pdb" % (working_pdb_path, entry)
            # # pdbExtractor.extract_residues([ligandBio], ligandPdbFile)
            # lig_sel = ResidueSelector({ligand})
            # lig_block = pdb_object_2block(structure, lig_sel,
            #                              write_conects=False)
            # rdLig = MolFromPDBBlock(lig_block)
            # Compute2DCoords(rdLig)

        except Exception as e:
            logger.exception(e)
