from os import remove
import math
import glob
import subprocess
import os.path
import warnings

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_protein

from luna.MyBio.util import parse_from_file, save_to_file
from luna.util.exceptions import InvalidSuperpositionFileError
from luna.util.file import is_directory_valid


import logging
logger = logging.getLogger()


TMALIGN = "/bin/tmalign"


class TMAlignment(MultipleSeqAlignment):

    """ Store TM-align results, including sequence alignment and TM-score.

    Parameters
    ----------
    score : float
        TM-score.
    records : iterable of :class:`Bio.SeqRecord.SeqRecord`
        Same length protein sequences. This may be an empty list.
    **kwargs : dict, optional
        Extra arguments to `TMAlignment`.
        Refer to :class:`Bio.Align.MultipleSeqAlignment` documentation for a
        list of all possible arguments.
    """

    def __init__(self, score, records, **kwargs):
        self.score = score
        super().__init__(records, **kwargs)


def align_structures(pdb_to_align, ref_pdb, output_path=None, tmalign=None):
    """Align two PDB files with TM-align.

    .. warning::

        TM-align performs pair-wise structural alignments.
        Therefore, if your PDB file contains multiple structures, it is
        recommended you extract the chains first and align them separately,
        otherwise the alignment may not produce the expected results.

        To extract chains you can use :class:`~luna.MyBio.extractor.Extractor`.

    Parameters
    ----------
    pdb_to_align : str
        The PDB structure that will be aligned to ``ref_pdb``.
    ref_pdb : str
        Reference PDB file, i.e., align ``pdb_to_align`` to ``ref_pdb``.
    output_path : str
        Where to save TM-align output structures.
    tmalign : str
        Pathname to TM-align binary. The default value is '/bin/tmalign'.

    Returns
    -------
     : `TMAlignment`
     Results of TM-align, including sequence alignment and TM-score.

    Examples
    --------

    >>> from luna.util.default_values import LUNA_PATH
    >>> from luna.align.tmalign import align_structures
    >>> pdb1 = "{LUNA_PATH}/tutorial/inputs/3QQF.pdb"
    >>> pdb2 = "{LUNA_PATH}/tutorial/inputs/3QQK.pdb"
    >>> alignment = align_structures(pdb1, pdb2, "./", \
tmalign="/media/data/Workspace/Softwares/TMalignc/TMalign")
    >>> print(alignment.score)
    0.98582
    """

    tmalign = tmalign or TMALIGN

    logger.info("TM-align will try to align the files %s and %s."
                % (pdb_to_align, ref_pdb))

    output = _run_tmalign(pdb_to_align, ref_pdb, output_path, tmalign)

    seq_pair, tm_score = get_seq_records(output, pdb_to_align, ref_pdb)

    alignment = TMAlignment(tm_score,
                            records=seq_pair,
                            alphabet=generic_protein)

    logger.info("TM-align finished successfully. %s."
                % alignment[0].description)

    return alignment


def _run_tmalign(file1, file2, output_path, tmalign):

    logger.debug("It will try to execute the command: '%s %s %s'."
                 % (tmalign, file1, file2))

    for fname in (file1, file2):
        if not os.path.isfile(fname):
            logging.error("Missing file: %s", fname)
            raise FileNotFoundError("Missing file: %s", fname)

    try:
        if output_path is not None and output_path.strip() != "":
            if is_directory_valid(output_path):
                logger.debug("The superposition files will be saved "
                             "at the directory '%s'" % output_path)

                filename = os.path.split(os.path.basename(file1))[1]
                output_file = "%s/%s.sup" % (output_path, filename)
                args = [tmalign, file1, file2, "-o", output_file]
        else:
            args = [tmalign, file1, file2]

        output = subprocess.check_output(args)
    except subprocess.CalledProcessError as e:
        logger.exception(e)

        raise RuntimeError("TMalign failed for PDB files: %s and %s"
                           % (file1, file2))

    return output.decode()


def get_seq_records(tm_output, aligned_pdb_id, ref_pdb_id):
    """Create a pair of :class:`Bio.SeqRecord.SeqRecord` from a
    TM-align output.

    Parameters
    ----------
    tm_output : str
        The output produced by TM-align.
    aligned_pdb_id : str
        An identifier for the aligned PDB file.
    ref_pdb_id : str
        An identifier for the reference PDB file.
    """
    logger.debug("Parsing the TMalign output.")

    lines = tm_output.splitlines()

    # Extract the TM-score (measure of structure similarity)
    # Take the mean of the (two) given TM-scores -- not sure which is reference
    tm_scores = []
    for line in lines:
        if line.startswith('TM-score'):
            # TMalign v. 2012/05/07 or earlier
            tm_scores.append(float(line.split(None, 2)[1]))
        elif 'TM-score=' in line:
            # TMalign v. 2013/05/11 or so
            tokens = line.split()
            for token in tokens:
                if token.startswith('TM-score='):
                    _key, _val = token.split('=')
                    tm_scores.append(float(_val.rstrip(',')))
                    break

    tm_score = math.fsum(tm_scores) / len(tm_scores)
    # Extract the sequence alignment
    last_lines = lines[-7:]

    assert last_lines[0].startswith('(":"')  # (":" denotes the residues pairs
    assert last_lines[-1].startswith('Total running time is')

    aligned_seq, ref_seq = last_lines[1].strip(), last_lines[3].strip()

    return (SeqRecord(Seq(aligned_seq),
                      id=aligned_pdb_id,
                      description="TM-score=%f" % tm_score),
            SeqRecord(Seq(ref_seq),
                      id=ref_pdb_id,
                      description="TM-score=%f" % tm_score)), tm_score


def extract_chain_from_sup(sup_file,
                           extract_chain,
                           output_file,
                           new_chain_id=None,
                           QUIET=True):
    """ Extract a chain from the superposition file generated by TM-align.

        .. warning::
            TM-align modifies the id of the original chains, so it is highly
            recommended to rename it to match the original chain id.
            To do so, use the parameter ``new_chain_id``.

        Parameters
        ----------
        sup_file: str
            A superposition file generated by TM-align.
        extract_chain: {'A', 'B'}
            Target chain id to be extracted. TM-align performs pair-wise
            structural alignment, so the output file will always contain
            two chains 'A' and 'B', where 'A' and 'B' are the aligned and
            reference structures, respectively.
        output_file: str
            Save the extracted chain to this file.
        new_chain_id: str, optional
            The new chain id of the extracted chain.
            If not provided, the chain ids will be maintained as it is in the
            TM-align output.
        QUIET: bool
            If True (the default), mute warning messages generated by
            Biopython.
    """
    if QUIET:
        try:
            warnings.filterwarnings("ignore")
            logger.debug("Quiet mode activated. From now on, no warning will "
                         "be printed.")
        except Exception:
            logger.warning("Quiet mode could not be activated.")

    if extract_chain not in ["A", "B"]:
        raise ValueError("Valid values for 'extract_chain' are 'A' and 'B'.")

    try:
        logger.debug("Trying to parse the file '%s'." % sup_file)

        structure = parse_from_file("SUP", sup_file)

        model = structure[0]
        if (len(model.child_list) != 2):
            error_msg = ("This structure has %d chains. The file generated "
                         "by TM-align must have 2 chains."
                         % len(model.child_list))
            raise InvalidSuperpositionFileError(error_msg)

        chain_to_remove = 'B' if extract_chain == "A" else "A"
        model.detach_child(chain_to_remove)

        if new_chain_id is not None and extract_chain != new_chain_id:
            model[extract_chain].id = new_chain_id
            logger.debug("Modifications completed.")

        save_to_file(structure, output_file)

        logger.debug("File '%s' created successfully." % output_file)
    except Exception as e:
        logger.exception(e)
        raise


def remove_sup_files(path):
    """ Remove all superposition files created by TM-align at a defined
    directory ``path``."""
    try:
        if (is_directory_valid(path)):
            targets = ['*.sup', '*.sup_atm', '*.sup_all',
                       '*.sup_all_atm', '*.sup_all_atm_lig']

            for target in targets:
                files = glob.glob('%s/%s' % (path, target))
                for file in files:
                    try:
                        remove(file)
                    except Exception as e:
                        logger.exception(e)
                        logger.warning("File %s not removed." % file)
    except Exception as e:
        logger.exception(e)
        raise
