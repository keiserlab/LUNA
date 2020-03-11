from math import ceil
from os.path import exists
from collections import defaultdict
from functools import wraps
import time
import multiprocessing as mp
import logging

# Open Babel
from openbabel.pybel import readfile
from openbabel.pybel import informats as OB_FORMATS
# RDKit
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem import MolFromPDBBlock, MolFromSmiles
# SQLAlchemy
from sqlalchemy.orm.exc import NoResultFound

# Local modules
from analysis.residues import InteractingResidues, get_interacting_residues
from analysis.summary import *
from database.loader import *
from database.luna_model import *
from database.helpers import *
from database.util import (get_ligand_tbl_join_filter, format_db_ligand_entries,
                           format_db_interactions, object_as_dict, get_default_mappers_list)
from mol.depiction import ligand_pharm_figure
from mol.clustering import cluster_fps_butina
from mol.features import FeatureExtractor
from mol.fingerprint import generate_fp_for_mols
from mol.entry import DBEntry, MolEntry
from mol.groups import AtomGroupPerceiver
from mol.interaction.contact import get_contacts_for_entity
from mol.interaction.calc import InteractionCalculator
from mol.interaction.conf import InteractionConf
from mol.interaction.fp.shell import ShellGenerator
from mol.wrappers.base import MolWrapper
from mol.wrappers.rdkit import RDKIT_FORMATS, read_multimol_file
from mol.amino_features import DEFAULT_AMINO_ATM_FEATURES
from util.default_values import *
from util.exceptions import *
from util.file import create_directory, get_file_format, get_unique_filename
from util.logging import new_logging_file, load_default_logging_conf
from util.config_parser import Config
from util import iter_to_chunks
from util.multiprocessing_logging import start_mp_handler

from MyBio.PDB.PDBParser import PDBParser
from MyBio.selector import ResidueSelector
from MyBio.util import download_pdb, entity_to_string, get_entity_from_entry


PDB_PARSER = PDBParser(PERMISSIVE=True, QUIET=True, FIX_ATOM_NAME_CONFLICT=True, FIX_OBABEL_FLAGS=False)


class StepControl:

    def __init__(self, step_id, num_subtasks, num_executed_subtasks, is_complete=False):

        self.step_id = step_id
        self.num_subtasks = num_subtasks
        self.num_executed_subtasks = num_executed_subtasks
        self.is_complete = is_complete

    @property
    def has_subtasks(self):
        return self.num_subtasks > 1

    @property
    def progress(self):
        if self.has_subtasks:
            return (self.num_executed_subtasks / self.num_subtasks)
        else:
            return 100 if self.is_complete else 0

    def update_progress(self, step=1):
        if self.has_subtasks:
            self.num_executed_subtasks += step
            if (self.num_subtasks ==
                    self.num_executed_subtasks):
                self.is_complete = True
        else:
            self.is_complete = True

    def __repr__(self):
        return '<StepControl: [Step id=%d, Progress=%.2f%%]>' % (self.step_id, self.progress)


class ExceptionWrapper(object):

    def __init__(self, step_id, has_subtasks, is_critical, desc=None):
        self.step_id = step_id
        self.has_subtasks = has_subtasks
        self.is_critical = is_critical
        self.desc = desc

    def __call__(self, func):
        @wraps(func)
        def callable(*args, **kwargs):
            # Project class
            proj_obj = args[0]

            # TODO: update if it is critical based on the database if it is available
            task = proj_obj.get_or_create_task(self.step_id, self.has_subtasks)
            error_message = None
            try:
                # Print a description of the function to be executed
                if self.desc:
                    desc = self.desc
                else:
                    get_default_message = True

                    # If the object is an instance of a project.
                    if hasattr(proj_obj, 'job_code'):
                        if self.step_id in proj_obj.step_details:
                            desc = proj_obj.step_details[self.step_id]["description"]
                            get_default_message = False

                    if get_default_message:
                        desc = ("Running function: %s." % func_call_2str(func, *args, **kwargs))
                proj_obj.logger.info(desc)

                # If the object is an instance of a project.
                if hasattr(proj_obj, 'job_code'):
                    proj_obj.update_step_details(task)

                result = func(*args, **kwargs)

                proj_obj.logger.warning("The function '%s' finished successfully." % func.__name__)

                return
            except Exception as e:
                proj_obj.logger.warning("The function '%s' failed." % func.__name__)
                proj_obj.logger.exception(e)

                # TODO: se chegar uma mensagem de erro não user friendly,
                # coloque uma mensagem generica
                error_message = e.args[0]
                if self.is_critical:
                    proj_obj.logger.warning("As the called function was critical, the program will be aborted.")
                    raise
            finally:
                task.update_progress()

                # If the object is an instance of a project,
                # update the step details in the project related tables.
                if hasattr(proj_obj, 'job_code') is False:
                    proj_obj.update_step_details(task, error_message)

                if task.is_complete:
                    proj_obj.logger.warning("The step %d has been completed." % task.step_id)

        return callable


class Project:

    def __init__(self,
                 entries,
                 working_path,
                 pdb_path=PDB_PATH,
                 has_local_files=False,
                 overwrite_path=False,
                 db_conf_file=None,
                 pdb_template=None,
                 atom_prop_file=ATOM_PROP_FILE,
                 try_h_addition=True,
                 ph=7.4,
                 amend_mol=True,
                 mol_obj_type='rdkit',
                 inter_conf=INTERACTION_CONF,
                 inter_calc=None,
                 calc_mfp=True,
                 mfp_opts=None,
                 calc_ifp=True,
                 ifp_num_levels=7,
                 ifp_radius_step=1,
                 ifp_length=IFP_LENGTH,
                 ifp_count=False,
                 similarity_func="BulkTanimotoSimilarity",
                 preload_mol_files=False,
                 default_properties=DEFAULT_AMINO_ATM_FEATURES,
                 butina_cutoff=0.2,
                 run_from_step=None,
                 run_until_step=None):

        if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
            raise IllegalArgumentError("Objects of type '%s' are not currently accepted. "
                                       "The available options are: %s." % (mol_obj_type, ", ".join(ACCEPTED_MOL_OBJ_TYPES)))

        if inter_conf is not None and isinstance(inter_conf, InteractionConf) is False:
            raise IllegalArgumentError("The informed interaction configuration must be an instance of '%s'." % InteractionConf)

        if inter_calc is not None and isinstance(inter_calc, InteractionCalculator) is False:
            raise IllegalArgumentError("The informed interaction configuration must be an instance of '%s'." % InteractionCalculator)

        self.entries = entries
        self.working_path = working_path
        self.pdb_path = pdb_path
        self.has_local_files = has_local_files
        self.overwrite_path = overwrite_path
        self.db_conf_file = db_conf_file
        self.pdb_template = pdb_template
        self.atom_prop_file = atom_prop_file
        self.ph = ph
        self.amend_mol = amend_mol
        self.mol_obj_type = mol_obj_type
        self.try_h_addition = try_h_addition

        # Interaction calculator parameters.
        self.inter_conf = inter_conf

        if inter_calc is None:
            inter_calc = InteractionCalculator(inter_conf=self.inter_conf)
        self.inter_calc = inter_calc

        # Fingerprint parameters.
        self.calc_mfp = calc_mfp
        self.mfp_opts = mfp_opts
        self.calc_ifp = calc_ifp
        self.ifp_num_levels = ifp_num_levels
        self.ifp_radius_step = ifp_radius_step
        self.ifp_length = ifp_length
        self.ifp_count = ifp_count

        self.similarity_func = similarity_func
        self.butina_cutoff = butina_cutoff
        self.default_properties = default_properties
        self.run_from_step = run_from_step
        self.run_until_step = run_until_step
        self.preload_mol_files = preload_mol_files
        self.step_controls = {}

        load_default_logging_conf()

        self.log_preferences()

    def __call__(self):
        raise NotImplementedError("This class is not callable. Use a class that implements this method.")

    def run(self):
        self()

    def log_preferences(self):
        logger.info("New project initialized...")
        params = ["\t\t-- %s = %s" % (key, str(self.__dict__[key])) for key in sorted(self.__dict__)]
        logger.info("Preferences:\n%s" % "\n".join(params))

    # @ExceptionWrapper(step_id=1, has_subtasks=False, is_critical=True)
    def prepare_project_path(self):
        logger.warning("Initializing path '%s'..." % self.working_path)

        create_directory(self.working_path, self.overwrite_path)
        create_directory("%s/pdbs" % self.working_path)
        create_directory("%s/figures" % self.working_path)
        create_directory("%s/results" % self.working_path)
        create_directory("%s/logs" % self.working_path)
        create_directory("%s/tmp" % self.working_path)

        logger.warning("Path '%s' created successfully!!!" % self.working_path)

    def init_logging_file(self, logging_filename=None, use_mp_handler=True):
        if not logging_filename:
            logging_filename = get_unique_filename(TMP_FILES)

        try:
            new_logging_file(logging_filename)

            if use_mp_handler:
                start_mp_handler()

            logger.warning("Logging file initialized successfully.")

            # Print preferences at the new logging file.
            self.log_preferences()
        except Exception as e:
            logger.exception(e)
            raise FileNotCreated("Logging file could not be created.")

    def init_db_connection(self):
        logger.info("A database configuration file was defined. Starting a new database connection...")

        config = Config(self.db_conf_file)
        dbConf = config.get_section_map("database")
        self.db = DBLoader(**dbConf)

    def init_common_tables(self):
        for conf in get_default_mappers_list(self.db):
            self.db.new_mapper(conf.table_class,
                               conf.table_name,
                               conf.properties)

    def get_or_create_task(self, step_id, has_subtasks):
        if step_id in self.step_controls:
            task = self.step_controls[step_id]
        else:
            num_subtasks = len(self.entries) if has_subtasks else 0
            num_executed_subtasks = 0
            task = StepControl(step_id, num_subtasks, num_executed_subtasks)
        return task

    def get_status_id(self, status_name):
        if status_name not in self.status_control:
            raise KeyError("The Status table does not have an %s entry for project control." % status_name)

        return self.status_control[status_name].id

    def update_step_details(self, task, error_message=None):
        status_id = self.get_status_id("RUNNING")
        warning = None
        if task.has_subtasks:
            if task.is_complete:
                status_id = self.get_status_id("DONE")

            if error_message:
                warning = "One or more entries failed."
                db_message = ProjectStepMessage(self.project_id, task.step_id, self.get_status_id("WARNING"), error_message)
                self.db.session.add(db_message)

                db_target_entry = (self.db.session
                                   .query(LigandEntryTable)
                                   .filter(LigandEntryTable.id == self.current_entry.id)
                                   .first())

                db_target_entry.step_messages.append(db_message)
        else:
            if task.is_complete:
                if error_message:
                    status_id = self.get_status_id("FAILED")
                else:
                    status_id = self.get_status_id("DONE")

        # TODO: verificar se temos um banco de dados para ser atualizado ou não
        db_step = ProjectStepDetail(self.project_id, task.step_id, status_id, warning, task.progress)

        self.db.session.merge(db_step)
        self.db.approve_session()

    def validate_entry_format(self, target_entry):
        if not target_entry.is_valid():
            raise InvalidEntry("Entry '%s' does not match a LUNA's entry format." % target_entry.to_string())

    def get_pdb_file(self, pdb_id):
        pdb_file = "%s/%s.pdb" % (self.pdb_path, pdb_id)

        if self.has_local_files:
            if not exists(pdb_file):
                raise FileNotFoundError("The PDB file '%s' was not found." % pdb_file)
        elif not exists(pdb_file):
            working_pdb_path = "%s/pdbs" % self.working_path
            pdb_file = "%s/%s.pdb" % (working_pdb_path, pdb_id)

            try:
                download_pdb(pdb_id=pdb_id, output_path=working_pdb_path)
            except Exception as e:
                logger.exception(e)
                raise FileNotCreated("PDB file '%s' was not created." % pdb_file) from e

        return pdb_file

    def decide_hydrogen_addition(self, pdb_header):
        if self.try_h_addition:
            if "structure_method" in pdb_header:
                method = pdb_header["structure_method"]
                # If the method is not a NMR type does not add hydrogen as it usually already has hydrogens.
                if method.upper() in NMR_METHODS:
                    logger.exception("The structure related to the entry '%s' was obtained by NMR, so it will "
                                     "not add hydrogens to it." % self.current_entry)
                    return False
            return True
        return False

    def perceive_chemical_groups(self, entity, ligand, add_h=False):
        perceiver = AtomGroupPerceiver(self.feature_extractor, add_h=add_h, ph=self.ph, amend_mol=self.amend_mol,
                                       mol_obj_type=self.mol_obj_type, default_properties=self.default_properties,
                                       tmp_path="%s/tmp" % self.working_path)

        radius = self.inter_conf.boundary_cutoff or BOUNDARY_CONF.boundary_cutoff
        nb_compounds = get_contacts_for_entity(entity, ligand, level='R', radius=radius)

        mol_objs_dict = {}
        if isinstance(self.current_entry, MolEntry):
            mol_objs_dict[self.current_entry.get_biopython_key()] = self.current_entry.mol_obj

        atm_grps_mngr = perceiver.perceive_atom_groups(set([x[1] for x in nb_compounds]), mol_objs_dict=mol_objs_dict)

        logger.info("Chemical group perception finished!!!")

        return atm_grps_mngr

    def set_pharm_objects(self):
        feature_factory = ChemicalFeatures.BuildFeatureFactory(self.atom_prop_file)

        # TODO: It should'be created unless pharm2d_fp is chosen.
        sig_factory = SigFactory(feature_factory, minPointCount=2, maxPointCount=3, trianglePruneBins=False)
        sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
        sig_factory.Init()
        self.sig_factory = sig_factory

        self.feature_extractor = FeatureExtractor(feature_factory)

    def get_fingerprint(self, rdmol):

        fp_opt = self.mfp_opts or {"fp_function": "pharm2d_fp", "sigFactory": self.sig_factory}
        fp_opt["critical"] = True

        return next(generate_fp_for_mols([rdmol], **fp_opt))

    def recover_rcsb_interactions(self, target_entry, manager):
        entry_str = target_entry.to_string()

        logger.info("Trying to select pre-computed interactions for the entry '%s'." % entry_str)

        db_ligand_entity = (self.db.session
                            .query(Ligand)
                            .filter(Ligand.pdb_id == target_entry.pdb_id and
                                    Ligand.chain_id == target_entry.chain_id and
                                    Ligand.lig_name == target_entry.lig_name and
                                    Ligand.lig_num == target_entry.lig_num and
                                    Ligand.lig_icode == target_entry.lig_icode)
                            .first())

        db_interactions = None
        if db_ligand_entity:
            # TODO: Mapear ligand_entity com status
            if db_ligand_entity.status.title == "AVAILABLE":
                join_filter = get_ligand_tbl_join_filter(target_entry, Ligand)
                db_interactions = (manager.select_interactions(join_filter,
                                                               interFilters))

                logger.info("%d pre-computed interaction(s) found in the "
                            "database for the entry '%s'."
                            % (len(filtered_inter), entry_str))
            else:
                logger.info("The entry '%s' exists in the database, but "
                            "there is no pre-computed interaction available. "
                            "So, LUNA will calculate the interactions to "
                            "this ligand." % entry_str)
        else:
            logger.info("The entry '%s' does not exist in the "
                        "database." % entry_str)

        return db_interactions

    def set_step_details(self):
        db_proj_type = (self.db.session.query(ProjectType).filter(ProjectType.type == self.project_type).first())

        if not db_proj_type:
            raise NoResultFound("No step details found for the project "
                                "type '%s'." % self.job_code)

        step_details = {}
        for db_step in db_proj_type.steps:
            step_details[db_step.id] = object_as_dict(db_step)

        self.step_details = step_details

    # @ExceptionWrapper(step_id=1, has_subtasks=False, is_critical=True)
    def set_project_id(self):
        db_proj = (self.db.session
                   .query(Project)
                   .filter(Project.job_code == self.job_code).one_or_none())

        if not db_proj:
            raise NoResultFound("Job code '%s' does not exist in the database."
                                % self.job_code)
        else:
            self.project_id = db_proj.id

    def recover_all_entries(self):
        db_ligand_entries = LigandEntryTableManager(self.db).get_entries(self.project_id)

        if not db_ligand_entries:
            raise NoResultFound("No ligand entries for the informed job code "
                                "'%s' was found." % self.job_code)

        return db_ligand_entries

    def set_status_control(self):
        db_status_entries = (self.db.session.query(Status).all())

        if not db_status_entries:
            raise NoResultFound("No status entries was found.")

        self.status_control = {r.name: r for r in db_status_entries}

    def get_rdkit_mol(self, entity, target, mol_name="Mol0"):
        target_sel = ResidueSelector({target})
        pdb_block = entity_to_string(entity, target_sel, write_conects=False)
        rdmol = MolFromPDBBlock(pdb_block)
        rdmol.SetProp("_Name", mol_name)

        return rdmol

    def add_mol_obj_to_entries(self):
        mol_files = defaultdict(dict)
        for entry in self.entries:
            if not entry.is_mol_obj_loaded():
                entry.mol_obj_type = self.mol_obj_type
                mol_files[(entry.mol_file, entry.mol_file_ext)][entry.mol_id] = entry
            else:
                logger.info("Molecular object in entry %s was manually informed and will not be reloaded." % entry)

        tool = "Open Babel" if self.mol_obj_type == "openbabel" else "RDKit"
        logger.info("It will try to preload the molecular objects using %s. Total of files to be read: %d."
                    % (tool, len(mol_files)))

        try:
            for mol_file, mol_file_ext in mol_files:
                key = (mol_file, mol_file_ext)
                ext = mol_file_ext or get_file_format(mol_file)

                available_formats = OB_FORMATS if self.mol_obj_type == "openbabel" else RDKIT_FORMATS
                if ext not in available_formats:
                    raise IllegalArgumentError("Extension '%s' informed or assumed from the filename is not a format "
                                               "recognized by %s." % (ext, tool))

                if not exists(mol_file):
                    raise FileNotFoundError("The file '%s' was not found." % mol_file)

                logger.info("Reading the file '%s'. The number of target entries located in this file is %d."
                            % (mol_file, len(mol_files[key])))

                try:
                    if self.mol_obj_type == "openbabel":
                        for ob_mol in readfile(ext, mol_file):
                            mol_id = ob_mol.OBMol.GetTitle().strip()
                            if mol_id in mol_files[key]:
                                entry = mol_files[key][mol_id]
                                entry.mol_obj = ob_mol
                                del(mol_files[key][mol_id])
                                logger.info("A structure to the entry '%s' was found in the file '%s' and loaded "
                                            "into the entry." % (entry, mol_file))

                                # If there is no other molecules to search, just stop the loop.
                                if len(mol_files[key]) == 0:
                                    break
                        else:
                            logger.info("All target ligands located in the file '%s' were successfully loaded." % mol_file)
                    else:
                        targets = list(mol_files[key].keys())
                        for (i, rdk_mol) in enumerate(read_multimol_file(mol_file, ext, targets=targets, removeHs=False)):
                            mol_id = targets[i]
                            entry = mol_files[key][mol_id]
                            # It returns None if the molecule parsing generated errors.
                            if rdk_mol:
                                entry.mol_obj = rdk_mol
                                del(mol_files[key][mol_id])
                                logger.info("A structure to the entry '%s' was found in the file '%s' and loaded "
                                            "into the entry." % (entry, mol_file))
                            else:
                                logger.warning("The structure related to the entry '%s' was found in the file '%s', but it could "
                                               "not be loaded as errors were found while parsing it." % (entry, mol_file))
                except Exception as e:
                    logger.exception(e)
                    raise MoleculeObjectError("An error occurred while parsing the molecular file '%s' with %s and the molecule "
                                              "objects could not be created. Check the logs for more information." %
                                              (mol_file, tool))
        except Exception as e:
            logger.exception(e)
            raise

        invalid_entries = [e for m in mol_files for e in mol_files[m].values()]
        if invalid_entries:
            entries = set(self.entries)
            for entry in invalid_entries:
                entries.remove(entry)
                logger.warning("Entry '%s' was not found or generated errors, so it will be removed from the entries list."
                               % entry)
            logger.warning("%d invalid entries removed." % len(invalid_entries))
            self.entries = entries

    def generate_ligand_figure(self, rdmol, group_types):
        atm_types = defaultdict(set)
        atm_map = {}
        for atm in rdmol.GetAtoms():
            atm_map[atm.GetPDBResidueInfo().GetSerialNumber()] = atm.GetIdx()

        for grp in group_types.atm_grps:
            for atm in grp.atoms:
                atm_id = atm_map[atm.serial_number]
                atm_types[atm_id].update(set(grp.chemicalFeatures))

        output = "%s/figures/%s.svg" % (self.working_path,
                                        rdmol.GetProp("_Name"))
        ligand_pharm_figure(rdmol, atm_types, output, ATOM_TYPES_COLOR)

    def clusterize_ligands(self, fingerprints):
        fps_only = [x["fp"] for x in fingerprints]

        try:
            clusters = cluster_fps_butina(fps_only, cutoff=self.butina_cutoff, similarity_func=self.similarity_func)
        except Exception:
            raise ProcessingFailed("Clustering step failed.")

        lig_clusters = {}
        for i, cluster in enumerate(clusters):
            for mol_id in cluster:
                lig_clusters[fingerprints[mol_id]["mol"]] = i

        return lig_clusters

    def store_compound_statistics(self, lig_entry_id, grp_types_count):
        for type_id in grp_types_count:
            self.db.session.merge(CompTypeCount(lig_entry_id, self.project_id, type_id, grp_types_count[type_id]))
        # TODO: Remover comentario
        # self.db.approve_session()

    def store_interaction_statistics(self, lig_entry_id, inter_type_count):
        for type_id in inter_type_count:
            self.db.session.merge(InterTypeCount(lig_entry_id, self.project_id, type_id, inter_type_count[type_id]))
        # TODO: Remover comentario
        # self.db.approve_session()

    def update_freq_by_cluster(self, count_dict, cluster_id, freq_by_cluster):
        for type_id in count_dict:
            key = (cluster_id, type_id, count_dict[type_id])
            freq_by_cluster[key] += 1

    def store_inter_type_freq(self, inter_type_freq_by_cluster):
        cumulative_num = defaultdict(int)
        keys = sorted(inter_type_freq_by_cluster)
        for key in keys:
            (cluster_id, inter_type_id, count) = key
            freq = inter_type_freq_by_cluster[key]
            cumulative_num[(cluster_id, inter_type_id)] += freq
            db_freq = InterTypeFreqByCluster(self.project_id, inter_type_id,
                                             cluster_id, count, freq,
                                             cumulative_num[(cluster_id, inter_type_id)])
            self.db.session.merge(db_freq)
        # TODO: Remover comentario
        # self.db.approve_session()

    def store_comp_type_freq(self, comp_type_freq_by_cluster):
        cumulative_num = defaultdict(int)
        keys = sorted(comp_type_freq_by_cluster)
        for key in keys:
            (cluster_id, comp_type_id, count) = key
            freq = comp_type_freq_by_cluster[key]
            cumulative_num[(cluster_id, comp_type_id)] += freq
            db_freq = CompTypeFreqByCluster(self.project_id, comp_type_id,
                                            cluster_id, count, freq,
                                            cumulative_num[(cluster_id, comp_type_id)])
            self.db.session.merge(db_freq)
        # TODO: Remover comentario
        # self.db.approve_session()


# POPULATE RCSB.
class RCSB_PLI_Population(Project):
    def __init__(self, entries, working_path, db_conf_file, **kwargs):

        super().__init__(entries=entries, working_path=working_path, db_conf_file=db_conf_file,
                         has_local_files=False, **kwargs)

    def __call__(self):
        self.init_logging_file()
        self.init_db_connection()

        self.prepare_project_path()
        self.init_common_tables()

        self.set_pharm_objects()

        rcsb_inter_manager = RCSBInteractionManager(self.db)

        # Loop over each entry.
        for target_entry in self.entries:
            mybio_ligand = ("H_%s" % target_entry.lig_name,
                            target_entry.lig_num,
                            target_entry.lig_icode)

            self.check_entry_existance(target_entry)
            self.validate_entry_format(target_entry)

            pdb_file = self.get_pdb_file(target_entry.pdb_id)
            pdb_file = self.add_hydrogen(pdb_file)

            structure = PDB_PARSER.get_structure(target_entry.pdb_id, pdb_file)
            ligand = structure[0][target_entry.chain_id][mybio_ligand]

            grps_by_compounds = self.perceive_chemical_groups(structure[0], ligand)
            src_grps = [grps_by_compounds[x] for x in grps_by_compounds
                        if x.get_id()[0] != " "]
            trgt_grps = [grps_by_compounds[x] for x in grps_by_compounds
                         if x != ligand]

            all_inter = calc_all_interactions(src_grps, trgt_grps, conf=BOUNDARY_CONF)

            rcsb_inter_manager.insert_interactions(all_inter, db_ligand_entity)

    def check_entry_existance(self, target_entry):
        db_ligand_entity = (self.db.session
                            .query(Ligand)
                            .filter(Ligand.pdb_id == target_entry.pdb_id and
                                    Ligand.chain_id == target_entry.chain_id and
                                    Ligand.lig_name == target_entry.lig_name and
                                    Ligand.lig_num == target_entry.lig_num and
                                    Ligand.lig_icode == target_entry.lig_icode)
                            .first())

        # If this entry does not exist in the database,
        # raise an error.
        if db_ligand_entity is None:
            message = ("Entry '%s' does not exist in the table 'ligand'." % target_entry.to_string())
            raise InvalidEntry(message)

        # If there are already interactions to this entry,
        # raise an error.
        elif status_by_id[db_ligand_entity.status_id] == "AVAILABLE":
            raise DuplicateEntry("Interactions to the entry '%s' already exists in the database."
                                 % target_entry.to_string())


class DB_PLI_Project(Project):

    def __init__(self, db_conf_file, job_code, working_path=None, entries=None, keep_all_potential_inter=True,
                 has_local_files=False, **kwargs):

        self.job_code = job_code

        if not working_path:
            working_path = "%s/projects/%s" % (LUNA_PATH, job_code)

        # If entries were not informed, get it from the database.
        self.get_entries_from_db = (entries is None)

        self.keep_all_potential_inter = keep_all_potential_inter

        self.has_local_files = has_local_files

        super().__init__(db_conf_file=db_conf_file, working_path=working_path, entries=entries, **kwargs)

    def __call__(self):

        # TODO: verificar o numero de entradas

        # TODO: Definir as mensagens de descrição para o usuário
        # TODO: colocar wrapers em cada função

        # TODO: add a new parameter: hydrogen addition mode to define how hydrogens will be added OR:
        #       add a new parameter: force hydrogen addition.

        self.init_logging_file()
        self.init_db_connection()
        self.init_common_tables()

        self.set_project_id()
        self.set_step_details()
        self.set_status_control()

        # Recover and set the entries related to this project
        if self.get_entries_from_db:
            db_lig_entries = self.recover_all_entries()
            db_lig_entries_by_id = {x.id: x for x in db_lig_entries}
            self.entries = format_db_ligand_entries(db_lig_entries)

        self.prepare_project_path()

        self.set_pharm_objects()

        rcsb_inter_manager = RCSBInteractionManager(self.db)
        proj_inter_manager = ProjectInteractionManager(self.db)

        db_comp_types = self.db.session.query(CompoundType).all()
        comp_type_id_map = {r.type: r.id for r in db_comp_types}
        db_inter_types = self.db.session.query(InterType).all()
        inter_type_id_map = {r.type: r.id for r in db_inter_types}

        # The calculus of interactions depend on the pharmacophore definition (ATOM_PROP_FILE) and the pH (PH).
        # Thus, if the user has kept such values unchanged and has not uploaded any PDB file, LUNA can try to
        # select the pre-computed interactions from a database.
        is_inter_params_default = (self.pdb_path is PDB_PATH and
                                   self.atom_prop_file == ATOM_PROP_FILE and
                                   self.ph == 7.4 and self.db_conf_file is not None)
        if is_inter_params_default:
            logger.info("Default configuration was kept unchanged. LUNA will try to select pre-computed "
                        "interactions from the defined database")

        fingerprints = []

        comp_type_count_by_entry = {}
        inter_type_count_by_entry = {}
        inter_res_by_entry = {}

        # Loop over each entry.
        for target_entry in self.entries:
            self.current_entry = target_entry

            # User has informed entries manually.
            # Check if the entries exist in the database.
            if self.get_entries_from_db is False:
                db_ligand_entity = self.recover_target_entry(target_entry)
            else:
                db_ligand_entity = db_lig_entries_by_id[target_entry.id]

            # Check if the entry is in the correct format.
            # It also accepts entries whose pdb_id is defined by
            # the filename.
            self.validate_entry_format(target_entry)

            pdb_file = self.get_pdb_file(target_entry.pdb_id)

            # TODO: resolver o problema dos hidrogenios.
            #       Vou adicionar hidrogenio a pH 7?
            #       Usuario vai poder definir pH?
            #       Obabel gera erro com posicoes alternativas

            structure = PDB_PARSER.get_structure(target_entry.pdb_id, pdb_file)
            ligand = get_entity_from_entry(structure, target_entry)

            # to_add_hydrogen = self.decide_hydrogen_addition(PDB_PARSER.get_header())
            # print(to_add_hydrogen)
            # exit()

            # exit()

            grps_by_compounds = self.perceive_chemical_groups(structure[0], ligand)
            src_grps = [grps_by_compounds[x] for x in grps_by_compounds if x.get_id()[0] != " "]
            trgt_grps = [grps_by_compounds[x] for x in grps_by_compounds if x != ligand]

            calc_interactions = True
            if is_inter_params_default:
                db_interactions = self.recover_rcsb_interactions(target_entry,
                                                                 rcsb_inter_manager)

                if db_interactions is not None:
                    filtered_inter = format_db_interactions(structure,
                                                            db_interactions)
                    calc_interactions = False

            if calc_interactions:
                if self.keep_all_potential_inter:
                    # TODO: calcular ponte de hidrogenio fraca
                    # TODO: detectar interações covalentes
                    # TODO: detectar complexos metalicos
                    # TODO: detectar CLASH

                    all_inter = calc_all_interactions(src_grps, trgt_grps, conf=BOUNDARY_CONF)

                    # Then it applies a filtering function.
                    filtered_inter = apply_interaction_criteria(all_inter, conf=self.inter_conf)
                    inter_to_be_stored = all_inter
                else:
                    # It filters the interactions by using a cutoff at once.
                    filtered_inter = calc_all_interactions(src_grps, trgt_grps, conf=self.inter_conf)
                    inter_to_be_stored = filtered_inter

                # TODO: Remover comentario
                # proj_inter_manager.insert_interactions(inter_to_be_stored, db_ligand_entity)

            # TODO: Ou uso STORE ou INSERT. Tenho que padronizar

            #
            # Count the number of compound types for each ligand
            #
            comp_type_count = count_group_types(grps_by_compounds[ligand], comp_type_id_map)
            # self.store_compound_statistics(target_entry.id,
            #                                comp_type_count)
            comp_type_count_by_entry[target_entry] = comp_type_count

            #
            # Count the number of interactions for each type
            #
            inter_type_count = count_interaction_types(filtered_inter, {ligand}, inter_type_id_map)
            # Store the count into the DB
            # self.store_interaction_statistics(target_entry.id, inter_type_count)
            inter_type_count_by_entry[target_entry] = inter_type_count

            # TODO: Atualizar outros scripts que verificam se um residuo é proteina ou nao.
            #       Agora temos uma função propria dentro de Residue.py
            interacting_residues = get_interacting_residues(filtered_inter, {ligand})
            inter_res_by_entry[target_entry] = interacting_residues

            # PAREI AQUI:
            # PROXIMO PASSO: contabilizar pra cada cluster o numero de ligantes que interagem com o residuo

            # TODO: Definir uma regĩao em volta do ligante pra ser o sitio ativo
            # TODO: Contar pra cada sitio o numero de cada tipo de grupo
            # TODO: Contar pra cada sitio o numero de interacoes

            # Select ligand and read it in a RDKit molecule object
            rdmol_lig = self.get_rdkit_mol(structure[0], ligand, target_entry.to_string())

            # Generate fingerprint for the ligand
            fp = self.get_fingerprint(rdmol_lig)
            fingerprints.append(fp)

            # self.generate_ligand_figure(rdmol_lig, grps_by_compounds[ligand])

        # TODO: Remover entradas repetidas
        # TODO: remover Entradas que derem erro

        # Clusterize ligands by using the Botina clustering method.
        clusters = self.clusterize_ligands(fingerprints)

        comp_type_freq_by_cluster = defaultdict(int)
        inter_type_freq_by_cluster = defaultdict(int)
        inter_res_freq_by_cluster = defaultdict(InteractingResidues)

        for target_entry in self.entries:
            cluster_id = clusters[target_entry.to_string()]

            # TODO: REMOVER
            cluster_id = 1

            # Get the ligand entry to update the cluster information
            if self.get_entries_from_db is False:
                db_ligand_entity = self.recover_target_entry(target_entry)
            else:
                db_ligand_entity = db_lig_entries_by_id[target_entry.id]

            db_ligand_entity.cluster = cluster_id

            comp_type_count = comp_type_count_by_entry[target_entry]
            # self.update_freq_by_cluster(comp_type_count, cluster_id, comp_type_freq_by_cluster)

            inter_type_count = inter_type_count_by_entry[target_entry]
            # self.update_freq_by_cluster(inter_type_count, cluster_id, inter_type_freq_by_cluster

            # TODO: Pra continuar fazendo isso vou precisar do alinhamento.
            # Vou fazer um alinhamento ficticio
            # self.update_res_freq_by_cluster(interacting_residues, inter_res_freq_by_cluster[cluster_id])

            # print("#######")
            # print(target_entry)
            # for k in interacting_residues.level1:
            #     print(k, k.__str__)
            # print()
            # for k in inter_res_freq_by_cluster[cluster_id].level1:
            #     print(k, inter_res_freq_by_cluster[cluster_id].level1[k])
            # print()
            # print()

        # Approve the cluster updates.
        self.db.approve_session()

        # Save the compound and interaction type frequency into the DB
        self.store_comp_type_freq(comp_type_freq_by_cluster)
        self.store_inter_type_freq(inter_type_freq_by_cluster)

        # TODO: alinhar proteína com o template
        # TODO: inserir alinhamento no BD

        # TODO: Calcular a frequencia dos residuos
        # TODO: Adicionar informações de frequencia no banco

    def update_res_freq_by_cluster(self, interacting_residues, inter_res_freq):
        inter_res_freq.update(interacting_residues)

    def recover_target_entry(self, target_entry):
        if isinstance(target_entry, DBEntry) is False:
            raise InvalidEntry("Entry '%s' is not a DBEntry object. So, LUNA cannot obtain "
                                     "an id information for this ligand entry and no database update will be "
                                     "possible." % entry_str)

        db_ligand_entity = (self.db.session
                            .query(LigandEntryTable)
                            .filter(LigandEntryTable.id == target_entry.id and
                                    LigandEntryTable.inter_proj_project_id == self.project_id and
                                    LigandEntryTable.pdb_id == target_entry.pdb_id and
                                    LigandEntryTable.chain_id == target_entry.chain_id and
                                    LigandEntryTable.lig_name == target_entry.lig_name and
                                    LigandEntryTable.lig_num == target_entry.lig_num and
                                    LigandEntryTable.lig_icode == target_entry.lig_icode)
                            .first())

        # If this entry does not exist in the database, raise an error.
        if db_ligand_entity is None:
            message = ("Entry '%s' does not exist in the database." % target_entry.to_string())
            raise InvalidNapoliEntry(message)

        return db_ligand_entity


class FingerprintProject(Project):

    def __init__(self, entries, working_path, mfp_output=None, ifp_output=None, **kwargs):
        self.mfp_output = mfp_output
        self.ifp_output = ifp_output

        super().__init__(entries=entries, working_path=working_path, **kwargs)

    def _process_entries(self, entries):

        for target_entry in entries:
            try:
                logger.info("Starting processing entry %s." % target_entry)
                self.current_entry = target_entry

                # # Check if the entry is in the correct format.
                # # It also accepts entries whose pdb_id is defined by the filename.
                if isinstance(target_entry, MolEntry) is False:
                    self.validate_entry_format(target_entry)

                if self.calc_ifp:
                    # # TODO: allow the person to pass a pdb_file into entries.
                    pdb_file = self.get_pdb_file(target_entry.pdb_id)

                    # # # TODO: resolver o problema dos hidrogenios.
                    # # #       Vou adicionar hidrogenio a pH 7?
                    # # #       Usuario vai poder definir pH?
                    structure = PDB_PARSER.get_structure(target_entry.pdb_id, pdb_file)
                    add_hydrogen = self.decide_hydrogen_addition(PDB_PARSER.get_header())

                    if isinstance(target_entry, MolEntry):
                        structure = target_entry.get_biopython_structure(structure, PDB_PARSER)

                    ligand = get_entity_from_entry(structure, target_entry)
                    ligand.set_as_target(is_target=True)

                    atm_grps_mngr = self.perceive_chemical_groups(structure[0], ligand, add_hydrogen)

                    #
                    # Calculate interactions
                    #
                    interactions_mngr = self.inter_calc.calc_interactions(atm_grps_mngr.atm_grps)

                    # Create hydrophobic islands.
                    atm_grps_mngr.merge_hydrophobic_atoms(interactions_mngr)

                result = {"id": (target_entry.to_string())}
                if self.calc_mfp:
                    if isinstance(target_entry, MolEntry):
                        rdmol_lig = MolFromSmiles(MolWrapper(target_entry.mol_obj).to_smiles())
                        rdmol_lig.SetProp("_Name", target_entry.mol_id)
                        result["mfp"] = generate_fp_for_mols([rdmol_lig], "morgan_fp")[0]["fp"]
                    else:
                        # Read from PDB.
                        pass

                if self.calc_ifp:
                    shells = ShellGenerator(self.ifp_num_levels, self.ifp_radius_step)
                    sm = shells.create_shells(atm_grps_mngr)

                    unique_shells = not self.ifp_count
                    result["ifp"] = sm.to_fingerprint(fold_to_size=self.ifp_length, unique_shells=unique_shells, count_fp=self.ifp_count)

                if isinstance(target_entry, MolEntry):
                    result["smiles"] = MolWrapper(target_entry.mol_obj).to_smiles()
                else:
                    # TODO: Get Smiles from PDB
                    result["smiles"] = ""
                    pass
                self.result.append(result)
                logger.warning("Processing of entry %s finished successfully." % target_entry)
            except Exception as e:
                logger.exception(e)
                logger.warning("Processing of entry %s failed. Check the logs for more information." % target_entry)

    def __call__(self):
        start = time.time()

        if not self.calc_mfp and not self.calc_ifp:
            logger.critical("Both molecular and interaction fingerprints were set off. So, there is nothing to be done...")
            return

        self.prepare_project_path()
        self.init_logging_file("%s/logs/project.log" % self.working_path)
        self.set_pharm_objects()

        fingerprints = []
        pli_fingerprints = []

        comp_type_count_by_entry = {}
        inter_type_count_by_entry = {}
        inter_res_by_entry = {}

        if self.preload_mol_files:
            self.add_mol_obj_to_entries()

        manager = mp.Manager()
        self.result = manager.list()

        start = time.time()
        chunk_size = ceil(len(self.entries) / (mp.cpu_count() - 1))
        chunks = iter_to_chunks(self.entries, chunk_size)

        processes = []
        for (i, l) in enumerate(chunks):
            p = mp.Process(name="Chunk %d" % i, target=self._process_entries, args=(l,))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()

        if self.calc_mfp:
            self.mfp_output = self.mfp_output or "%s/results/mfp.csv" % self.working_path
            with open(self.mfp_output, "w") as OUT:
                OUT.write("ligand_id,smiles,on_bits\n")
                for r in self.result:
                    if "mfp" in r:
                        fp_str = "\t".join([str(x) for x in r["mfp"].GetOnBits()])
                        OUT.write("%s,%s,%s\n" % (r["id"], r["smiles"], fp_str))

        if self.calc_ifp:
            self.ifp_output = self.ifp_output or "%s/results/ifp.csv" % self.working_path

            with open(self.ifp_output, "w") as OUT:

                if self.ifp_count:
                    OUT.write("ligand_id,smiles,on_bits,count\n")
                else:
                    OUT.write("ligand_id,smiles,on_bits\n")

                for r in self.result:
                    if "ifp" in r:

                        if self.ifp_count:
                            fp_bits_str = "\t".join([str(idx) for idx in r["ifp"].counts.keys()])
                            fp_count_str = "\t".join([str(count) for count in r["ifp"].counts.values()])
                            OUT.write("%s,%s,%s,%s\n" % (r["id"], r["smiles"], fp_bits_str, fp_count_str))
                        else:
                            fp_bits_str = "\t".join([str(x) for x in r["ifp"].get_on_bits()])
                            OUT.write("%s,%s,%s\n" % (r["id"], r["smiles"], fp_bits_str))

        end = time.time()
        logger.info("Processing finished!!!")
        logger.info("Check the results at '%s'." % self.working_path)
        logger.info("Processing time: %.2fs." % (end - start))


class LocalProject(Project):

    def __init__(self, entries, working_path, has_local_files=False, **kwargs):
        super().__init__(entries=entries, working_path=working_path, has_local_files=has_local_files, **kwargs)

    def __call__(self):
        start = time.time()

        # TODO: verificar o numero de entradas

        # TODO: Definir as mensagens de descrição para o usuário
        # TODO: colocar wrapers em cada função

        # TODO: Remove duplicated entries.

        # TODO: add a new parameter: hydrogen addition mode to define how hydrogens will be added OR:
        #       add a new parameter: force hydrogen addition.

        # TODO: logging is not working. The logs are not being directed to the logging.log file.

        self.prepare_project_path()
        self.init_logging_file("%s/logs/project.log" % self.working_path)
        self.set_pharm_objects()

        self.mfps = []
        self.ifps = []
        self.interactions = []
        self.neighborhoods = []

        comp_type_count_by_entry = {}
        inter_type_count_by_entry = {}
        inter_res_by_entry = {}

        if self.preload_mol_files:
            self.add_mol_obj_to_entries()

        # Loop over each entry.
        for target_entry in self.entries:
            logger.info("Processing entry: %s." % target_entry)

            self.current_entry = target_entry

            # Check if the entry is in the correct format.
            # It also accepts entries whose pdb_id is defined by the filename.
            if isinstance(target_entry, MolEntry) is False:
                self.validate_entry_format(target_entry)

            # TODO: allow the user to pass a pdb_file into entries.
            pdb_file = self.get_pdb_file(target_entry.pdb_id)
            target_entry.pdb_file = pdb_file

            # # TODO: resolver o problema dos hidrogenios.
            # #       Vou adicionar hidrogenio a pH 7?
            # #       Usuario vai poder definir pH?
            structure = PDB_PARSER.get_structure(target_entry.pdb_id, pdb_file)
            add_hydrogen = self.decide_hydrogen_addition(PDB_PARSER.get_header())

            if isinstance(target_entry, MolEntry):
                structure = target_entry.get_biopython_structure(structure, PDB_PARSER)

            ligand = get_entity_from_entry(structure, target_entry)
            ligand.set_as_target(is_target=True)

            atm_grps_mngr = self.perceive_chemical_groups(structure[0], ligand, add_hydrogen)

            self.neighborhoods.append((target_entry, atm_grps_mngr))

            #
            # Calculate interactions
            #
            interactions_mngr = self.inter_calc.calc_interactions(atm_grps_mngr.atm_grps)

            # Create hydrophobic islands.
            atm_grps_mngr.merge_hydrophobic_atoms(interactions_mngr)

            self.interactions.append((target_entry, interactions_mngr))
            continue

            inter_file = "%s/results/%s.tsv" % (self.working_path, target_entry.to_string())
            with open(inter_file, "w") as OUT:
                OUT.write("atom_group1\tfeatures1\tatom_group2\tfeatures2\tinteraction_type\n")
                for i in interactions:
                    src_grp = ", ".join([x.full_atom_name for x in i.src_grp.atoms])
                    trgt_grp = ", ".join([x.full_atom_name for x in i.trgt_grp.atoms])
                    feat1 = ", ".join([x.name for x in i.src_grp.features if x.name != "Atom"])
                    feat2 = ", ".join([x.name for x in i.trgt_grp.features if x.name != "Atom"])
                    OUT.write("%s\t%s\t%s\t%s\t%s\n" % (src_grp, feat1, trgt_grp, feat2, i.type))

            shells = ShellGenerator(self.ifp_num_levels, self.ifp_radius_step)
            sm = shells.create_shells(neighborhood)
            fp = sm.to_fingerprint(fold_to_size=self.ifp_length)

            self.ifps.append((target_entry, fp))

            #
            # Count the number of compound types for each ligand
            #
            # comp_type_count = count_group_types(grps_by_compounds[ligand])
            # comp_type_count_by_entry[target_entry] = comp_type_count

            #
            # Count the number of interactions for each type
            #
            # inter_type_count = count_interaction_types(filtered_inter,
            #                                            {ligand})
            # inter_type_count_by_entry[target_entry] = inter_type_count

            # TODO: Atualizar outros scripts que verificam se um residuo é proteina ou nao.
            #       Agora temos uma função propria dentro de Residue.py
            # interacting_residues = get_interacting_residues(filtered_inter, {ligand})
            # inter_res_by_entry[target_entry] = interacting_residues

            # PAREI AQUI:
            # PROXIMO PASSO: contabilizar pra cada cluster o numero de ligantes que interagem com o residuo

            # TODO: Definir uma regĩao em volta do ligante pra ser o sitio ativo
            # TODO: Contar pra cada sitio o numero de cada tipo de grupo
            # TODO: Contar pra cada sitio o numero de interacoes

            # Select ligand and read it in a RDKit molecule object.

            # TODO: I cannot use RDKit directly in the PDB file or it will crash the molecule.
            #       First I need to convert it to MOL.
            # rdmol_lig = self.get_rdkit_mol(structure[0], ligand, target_entry.to_string())

            # Generate fingerprint for the ligand
            # fp = self.get_fingerprint(rdmol_lig)
            # fingerprints.append(fp)

            # Generate fingerprint for the ligand.

            # if isinstance(target_entry, MolEntry):
            #     rdmol_lig = MolFromSmiles(str(target_entry.mol_obj).split("\t")[0])
            #     rdmol_lig.SetProp("_Name", target_entry. mol_id)
            #     fps = generate_fp_for_mols([rdmol_lig], "morgan_fp")
            #     fingerprints.append((target_entry, fps[0]["fp"]))

            # self.generate_ligand_figure(rdmol_lig, grps_by_compounds[ligand])

        # with open("ecfp4_fingerprints.csv", "w") as OUT:
        #     OUT.write("id,smarts,fp\n")
        #     for (entry, fp) in fingerprints:
        #         fp_str = "\t".join([str(x) for x in fp.GetOnBits()])
        #         OUT.write("%s,%s,%s\n" % (entry.mol_id, str(entry.mol_obj).split("\t")[0], fp_str))

        return

        with open("%s/results/ifp_fingerprints.csv" % self.working_path, "w") as OUT:
            OUT.write("id,fp\n")
            for (entry, fp) in pli_fingerprints:
                fp_str = "\t".join([str(x) for x in fp.get_on_bits()])
                OUT.write("%s,%s\n" % (entry.to_string(), fp_str))

        end = time.time()
        logger.info("Processing finished!!!")
        logger.info("Check the results at '%s'." % self.working_path)
        logger.info("Processing time: %.2fs." % (end - start))
        exit()

        # TODO: Remover entradas repetidas
        # TODO: remover Entradas que derem erro

        # Clusterize ligands by using the Botina clustering method.
        # clusters = self.clusterize_ligands(fingerprints)

        # comp_type_freq_by_cluster = defaultdict(int)
        # inter_type_freq_by_cluster = defaultdict(int)
        # inter_res_freq_by_cluster = defaultdict(InteractingResidues)

        # for target_entry in self.entries:
        #     cluster_id = clusters[target_entry.to_string()]

        #     # TODO: REMOVER
        #     cluster_id = 1

        #     # Get the ligand entry to update the cluster information
        #     if self.get_entries_from_db is False:
        #         db_ligand_entity = self.recover_target_entry(target_entry)
        #     else:
        #         db_ligand_entity = db_lig_entries_by_id[target_entry.id]

        #     db_ligand_entity.cluster = cluster_id

        #     comp_type_count = comp_type_count_by_entry[target_entry]
        #     # self.update_freq_by_cluster(comp_type_count, cluster_id,
        #                                 # comp_type_freq_by_cluster)

        #     inter_type_count = inter_type_count_by_entry[target_entry]
            # self.update_freq_by_cluster(inter_type_count, cluster_id,
                                        # inter_type_freq_by_cluster

            # TODO: Pra continuar fazendo isso vou precisar do alinhamento.
            # Vou fazer um alinhamento ficticio
            # self.update_res_freq_by_cluster(interacting_residues,
            #                                 inter_res_freq_by_cluster[cluster_id])

            # print("#######")
            # print(target_entry)
            # for k in interacting_residues.level1:
            #     print(k, k.__str__)
            # print()
            # for k in inter_res_freq_by_cluster[cluster_id].level1:
            #     print(k, inter_res_freq_by_cluster[cluster_id].level1[k])
            # print()
            # print()

        # TODO: alinhar proteína com o template
        # TODO: inserir alinhamento no BD

        # TODO: Calcular a frequencia dos residuos
        # TODO: Adicionar informações de frequencia no banco

    def update_res_freq_by_cluster(self, interacting_residues, inter_res_freq):
        inter_res_freq.update(interacting_residues)
