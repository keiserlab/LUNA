from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_protein

import math
import os.path
import subprocess

import logging
logger = logging.getLogger(__name__)

TMALIGN = "/bin/tmalign"


def run_tmalign(file1, file2, tmalign):

    logger.info("Trying to execute the command: '%s %s %s'",
                tmalign, file1, file2)

    for fname in (file1, file2):
        if not os.path.isfile(fname):
            logging.error("Missing file: %s", fname)
            raise FileExistsError("Missing file: %s", fname)

    try:
        output = subprocess.check_output([tmalign, file1, file2])
    except subprocess.CalledProcessError as e:
        logger.exception("%s TMalign failed (returned %s):\n%s"
                         % (e.returncode, e.output))

        raise RuntimeError("TMalign failed for PDB files: %s %s"
                           % (file1, file2))

    return output.decode()


def get_seq_records(tmOutput, refId, eqvId):
    """Create a pair of SeqRecords from TMalign output."""

    logger.info("Parsing the TMalign output")

    lines = tmOutput.splitlines()

    # Extract the TM-score (measure of structure similarity)
    # Take the mean of the (two) given TM-scores -- not sure which is reference
    tmScores = []
    for line in lines:
        if line.startswith('TM-score'):
            # TMalign v. 2012/05/07 or earlier
            tmScores.append(float(line.split(None, 2)[1]))
        elif 'TM-score=' in line:
            # TMalign v. 2013/05/11 or so
            tokens = line.split()
            for token in tokens:
                if token.startswith('TM-score='):
                    _key, _val = token.split('=')
                    tmScores.append(float(_val.rstrip(',')))
                    break

    tmScore = math.fsum(tmScores) / len(tmScores)
    # Extract the sequence alignment
    lastLines = lines[-7:]

    assert lastLines[0].startswith('(":"')  # (":" denotes the residues pairs
    assert lastLines[-1].startswith('Total running time is')

    refSeq, eqvSeq = lastLines[1].strip(), lastLines[3].strip()

    return (SeqRecord(Seq(refSeq), id=refId,
                      description="TMalign TM-score=%f" % tmScore),
            SeqRecord(Seq(eqvSeq), id=eqvId,
                      description="TMalign TM-score=%f" % tmScore),
            )


def align_2struct(file1, file2, tmalign=None):

    if (tmalign is None):
        tmalign = TMALIGN

    tmOutput = run_tmalign(file1, file2, tmalign)

    tmSeqPair = get_seq_records(tmOutput, file1, file2)

    alignment = MultipleSeqAlignment(tmSeqPair, generic_protein)

    return alignment
