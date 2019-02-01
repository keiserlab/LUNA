from util.exceptions import (BitsValueError, InvalidFingerprintType, IllegalArgumentError, FingerprintCountsError)

from rdkit.DataStructs.cDataStructs import (ExplicitBitVect, SparseBitVect)
from scipy.sparse import (issparse, csr_matrix)

from collections import defaultdict

import numpy as np
import logging


logger = logging.getLogger(__name__)

DEFAULT_SHELL_NBITS = 2**32
DEFAULT_FP_LENGTH = 1024
DEFAULT_FP_DTYPE = np.int32


class Fingerprint:

    def __init__(self, indices, fp_length=DEFAULT_SHELL_NBITS, unfolded_fp=None, unfolding_map=None, props=None):

        indices = np.asarray(indices, dtype=np.long)

        if np.any(np.logical_or(indices < 0, indices >= fp_length)):
            logger.exception("Provided indices are in a different bit scale.")
            raise BitsValueError("Provided indices are in a different bit scale.")

        self._indices = np.unique(indices)
        self._fp_length = fp_length
        self._unfolded_fp = unfolded_fp
        self._unfolding_map = unfolding_map or {}
        self._props = props or {}

    @classmethod
    def from_indices(cls, indices, fp_length=DEFAULT_SHELL_NBITS, **kwargs):
        return cls(indices, fp_length, **kwargs)

    @classmethod
    def from_vector(cls, vector, fp_length=None, **kwargs):
        if fp_length is None:
            try:
                fp_length = vector.shape[1]
            except IndexError:
                fp_length = vector.shape[0]

        if issparse(vector):
            indices = vector.indices.astype(np.long)
        else:
            indices = np.asarray(np.where(vector), dtype=np.long).flatten()

        return cls.from_indices(indices, fp_length, **kwargs)

    @classmethod
    def from_bit_string(cls, bit_string, fp_length=None, **kwargs):
        indices = [i for i, char in enumerate(bit_string) if char != '0']
        if fp_length is None:
            fp_length = len(bit_string)

        return cls.from_indices(indices, fp_length, **kwargs)

    @classmethod
    def from_rdkit(cls, rdkit_fp, **kwargs):
        if not (isinstance(rdkit_fp, ExplicitBitVect) or isinstance(rdkit_fp, SparseBitVect)):
            logger.exception("Invalid fingerprint type. RDKit only accepts a SparseBitVect or ExplicitBitVect object.")
            raise TypeError("Invalid fingerprint type. RDKit only accepts a SparseBitVect or ExplicitBitVect object.")

        fp_length = rdkit_fp.GetNumBits()
        indices = np.asarray(rdkit_fp.GetOnBits(), dtype=np.long)

        return cls.from_indices(indices, fp_length, **kwargs)

    @classmethod
    def from_fingerprint(cls, fp, **kwargs):
        if not isinstance(fp, Fingerprint):
            logger.exception("Informed fingerprint is not an instance of %s." % (cls.__class__))
            raise InvalidFingerprintType("Informed fingerprint is not an instance of %s." % (cls.__class__))

        unfolded_fp = fp.__class__.from_fingerprint(fp.unfolded_fp) if fp.unfolded_fp is not None else None

        unfolding_map = dict(fp.unfolding_map)
        props = dict(fp.props)

        new_fp = cls.from_indices(fp.indices, fp.fp_length, unfolded_fp=unfolded_fp,
                                  unfolding_map=unfolding_map, props=props)

        return new_fp

    @property
    def indices(self):
        return self._indices

    @property
    def bit_count(self):
        return self.indices.shape[0]

    @property
    def density(self):
        return self.bit_count / self.fp_length

    @property
    def counts(self):
        return dict([(k, 1) for k in self.indices])

    @property
    def fp_length(self):
        return self._fp_length

    @property
    def unfolded_fp(self):
        if self._unfolded_fp is None:
            logger.warning("This fingerprint was not previously folded.")
            return None
        return self._unfolded_fp

    @property
    def unfolded_indices(self):
        if self._unfolding_map is None:
            logger.warning("This fingerprint was not previously folded.")
            return None
        return self.unfolded_fp.indices

    @property
    def unfolding_map(self):
        if self._unfolding_map is None:
            logger.warning("This fingerprint was not previously folded.")
            return None
        return self._unfolding_map

    @property
    def props(self):
        return self._props

    @property
    def name(self):
        return self.props.get("name", "")

    @name.setter
    def name(self, name):
        self.props["name"] = str(name)

    @property
    def num_levels(self):
        return self.props.get("num_levels", None)

    @num_levels.setter
    def num_levels(self, num_levels):
        self.props["num_levels"] = str(num_levels)

    @property
    def num_shells(self):
        return self.props.get("num_shells", None)

    @num_shells.setter
    def num_shells(self, num_shells):
        self.props["num_shells"] = str(num_shells)

    def get_prop(self, key):
        try:
            return self.props[key]
        except KeyError:
            logger.warning("Key '%s' does not exist." % key)
            return None

    def set_prop(self, key, value):
        self.props[key] = value

    def get_num_bits(self):
        return self.fp_length

    def get_num_on_bits(self):
        return self.bit_count

    def get_num_off_bits(self):
        return self.get_num_bits() - self.get_num_on_bits()

    def get_bit(self, index):
        if index in self.counts:
            return self.counts[index]
        elif index >= 0 and index < self.fp_length:
            return 0
        else:
            logger.exception("Provided indices are in a different bit scale.")
            raise BitsValueError("Provided index is in a different bit scale.")

    def get_on_bits(self):
        return np.array([k for (k, v) in self.counts.items() if v > 0])

    def to_vector(self, compressed=True, dtype=DEFAULT_FP_DTYPE):
        data = [self.counts[i] for i in self.indices]
        if compressed:
            try:
                row = np.zeros(self.bit_count)
                col = self.indices
                vector = csr_matrix((data, (row, col)), shape=(1, self.fp_length), dtype=dtype)
            except ValueError as e:
                logger.exception(e)
                raise BitsValueError("Sparse matrix construction failed. Invalid indices or input data.")
        else:
            try:
                # This function is causing a MemoryError exception when using a 2**32 vector.
                vector = np.zeros(self.fp_length, dtype=dtype)
            except MemoryError as e:
                logger.exception(e)
                raise MemoryError("Huge indices vector detected. An operation ran out of memory. "
                                  "Use a different data type or apply a folding operation.")

            try:
                vector[self.indices] = data
            except IndexError as e:
                logger.exception(e)
                raise BitsValueError("Some of the provided indices are greater than the fingerprint length.")

        return vector

    def to_bit_vector(self, compressed=True):
        return self.to_vector(compressed=compressed, dtype=np.bool_).astype(np.int8)

    def to_bit_string(self):
        try:
            # This function is causing a MemoryError exception when using a 2**32 vector.
            bit_vector = self.to_bit_vector(compressed=False).astype(np.int8)
            return "".join(map(str, bit_vector))
        except MemoryError as e:
            logger.exception(e)
            raise MemoryError("Huge indices vector detected. An operation ran out of memory. "
                              "Use a different data type or apply a folding operation.")

    def to_rdkit(self, rdkit_fp_cls=None):

        if rdkit_fp_cls is None:
            # Classes to store explicit bit vectors: ExplicitBitVect or SparseBitVect.
            # ExplicitBitVect is most useful for situations where the size of the vector is
            # relatively small (tens of thousands or smaller).
            # For larger vectors, use the _SparseBitVect_ class instead.
            if self.fp_length < 1e5:
                rdkit_fp_cls = ExplicitBitVect
            else:
                rdkit_fp_cls = SparseBitVect

        # RDKit data structure defines fingerprints as a std:set composed of ints (signed int).
        # Since we always have values higher than 0 and since the data structure contains only signed ints,
        # then the max length for a RDKit fingerprint is 2^31 - 1.
        # C signed int (32 bit) ranges: [-2^31, 2^31-1].
        max_rdkit_fp_length = 2**31 - 1
        fp_length = self.fp_length
        if max_rdkit_fp_length < fp_length:
            logger.warning("The current fingerprint will be folded as its size is higher than the maximum "
                           "size accepted by RDKit, which is 2**31 - 1.")
            fp_length = max_rdkit_fp_length
        indices = self.indices % max_rdkit_fp_length

        rdkit_fp = rdkit_fp_cls(fp_length)
        rdkit_fp.SetBitsFromList(indices.tolist())
        return rdkit_fp

    def fold(self, new_fp_length=DEFAULT_FP_LENGTH):

        if new_fp_length > self.fp_length:
            error_msg = ("Fold operation requires the new fingerprint length (%d) "
                         "to be greater than the current one (%d)." % (new_fp_length, self.fp_length))
            logger.exception(error_msg)
            raise BitsValueError(error_msg)

        if not np.log2(self.fp_length / new_fp_length).is_integer():
            error_msg = ("It is not possible to fold the current fingerprint into the informed new length."
                         "The current length divided by the new one is not a power of 2 number.")
            logger.exception(error_msg)
            raise BitsValueError()

        folded_indices = self.indices % new_fp_length

        unfolding_map = defaultdict(set)
        for k, v in sorted(zip(folded_indices, self.indices)):
            unfolding_map[k].add(v)

        props = dict(self.props)
        new_fp = self.__class__(indices=folded_indices, fp_length=new_fp_length,
                                unfolded_fp=self, unfolding_map=unfolding_map, props=props)

        return new_fp

    def unfold(self):
        return self.unfolded_fp

    def union(self, other):
        if not isinstance(other, Fingerprint):
            logger.exception("Informed fingerprint is not an instance of %s." % (other.__class__))
            raise InvalidFingerprintType("Informed fingerprint is not an instance of %s." % (other.__class__))

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.union1d(self.indices, other.indices)

    def intersection(self, other):
        if not isinstance(other, Fingerprint):
            logger.exception("Informed fingerprint is not an instance of %s." % (other.__class__))
            raise InvalidFingerprintType("Informed fingerprint is not an instance of %s." % (other.__class__))

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.intersect1d(self.indices, other.indices, assume_unique=True)

    def difference(self, other):
        if not isinstance(other, Fingerprint):
            logger.exception("Informed fingerprint is not an instance of %s." % (other.__class__))
            raise InvalidFingerprintType("Informed fingerprint is not an instance of %s." % (other.__class__))

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.setdiff1d(self.indices, other.indices, assume_unique=True)

    def symmetric_difference(self, other):
        if not isinstance(other, Fingerprint):
            logger.exception("Informed fingerprint is not an instance of %s." % (other.__class__))
            raise InvalidFingerprintType("Informed fingerprint is not an instance of %s." % (other.__class__))

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.setxor1d(self.indices, other.indices, assume_unique=True)

    def __repr__(self):
        return ("<%s: indices=%s length=%d>" %
                (self.__class__, repr(self.indices).replace('\n', '').replace(' ', ''), self.fp_length))

    def __eq__(self, other):
        if isinstance(other, Fingerprint):
            return (self.__class__ == other.__class__ and
                    self.fp_length == other.fp_length and
                    np.all(np.in1d(self.indices, other.indices, assume_unique=True)))
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __or__(self, other):
        return self.union(other)

    def __and__(self, other):
        return self.intersection(other)

    def __sub__(self, other):
        return self.difference(other)

    def __xor__(self, other):
        return self.symmetric_difference(other)


class CountFingerprint(Fingerprint):

    def __init__(self, indices=None, counts=None, fp_length=DEFAULT_SHELL_NBITS,
                 unfolded_fp=None, unfolding_map=None, props=None):

        if indices is None and counts is None:
            logger.exception("Indices or counts must be provided.")
            raise IllegalArgumentError("Indices or counts must be provided.")

        if indices is not None:
            indices = np.asarray(indices, dtype=np.long)

            if np.any(np.logical_or(indices < 0, indices >= fp_length)):
                logger.exception("Provided indices are in a different bit scale.")
                raise BitsValueError("Provided indices are in a different bit scale.")

            if counts is None:
                indices, counts = np.unique(indices, return_counts=True)
                counts = dict(zip(indices, counts))
            else:
                indices = np.unique(indices)
                if not np.all([x in indices for x in counts]):
                    logger.exception("At least one index from 'counts' is not in 'indices'.")
                    raise FingerprintCountsError("At least one index from 'counts' is not in 'indices'.")
                if len(set(indices).symmetric_difference(counts)) > 0:
                    logger.exception("At least one index in 'indices' is not in 'counts'.")
                    raise FingerprintCountsError("At least one index in 'indices' is not in 'counts'.")
        else:
            indices = np.asarray(sorted(counts.keys()), dtype=np.long)

            if np.any(np.logical_or(indices < 0, indices >= fp_length)):
                logger.exception("Provided indices are in a different bit scale.")
                raise BitsValueError("Provided indices are in a different bit scale.")

        self._counts = counts
        super().__init__(indices, fp_length, unfolded_fp, unfolding_map, props)

    @classmethod
    def from_indices(cls, indices, counts=None, fp_length=DEFAULT_SHELL_NBITS, **kwargs):
        return cls(indices=indices, counts=counts, fp_length=fp_length, **kwargs)

    @classmethod
    def from_counts(cls, counts, fp_length=DEFAULT_SHELL_NBITS, **kwargs):
        return cls(counts=counts, fp_length=fp_length, **kwargs)

    @classmethod
    def from_bit_string(cls, bit_string, counts=None, fp_length=None, **kwargs):
        indices = [i for i, char in enumerate(bit_string) if char != '0']
        if fp_length is None:
            fp_length = len(bit_string)

        return cls.from_indices(indices, counts, fp_length, **kwargs)

    @classmethod
    def from_vector(cls, vector, fp_length=None, **kwargs):
        if fp_length is None:
            try:
                fp_length = vector.shape[1]
            except IndexError:
                fp_length = vector.shape[0]

        if issparse(vector):
            indices = vector.indices.astype(np.long)
            counts = vector.data
        else:
            indices = np.asarray(np.where(vector), dtype=np.long).flatten()
            counts = vector[indices]
        counts = dict(zip(indices, counts))

        return cls.from_indices(indices, counts, fp_length, **kwargs)

    @classmethod
    def from_fingerprint(cls, fp, **kwargs):
        if not isinstance(fp, Fingerprint):
            logger.exception("Informed fingerprint is not an instance of %s." % (cls.__class__))
            raise InvalidFingerprintType("Informed fingerprint is not an instance of %s." % (cls.__class__))

        counts = dict([(i, c) for i, c in fp.counts.items() if c > 0])
        unfolded_fp = fp.__class__.from_fingerprint(fp.unfolded_fp) if fp.unfolded_fp is not None else None
        unfolding_map = dict(fp.unfolding_map)
        props = dict(fp.props)

        new_fp = cls.from_counts(counts, fp.fp_length, unfolded_fp=unfolded_fp,
                                 unfolding_map=unfolding_map, props=props)

        return new_fp

    @property
    def counts(self):
        return self._counts

    def get_count(self, index):
        return self.counts.get(index, 0)

    def fold(self, new_fp_length=DEFAULT_FP_LENGTH):

        new_fp = super().fold(new_fp_length)

        new_fp._counts = dict([(folded_idx, sum([self.get_count(x) for x in unfolded_set]))
                              for folded_idx, unfolded_set in new_fp.unfolding_map.items()])

        return new_fp

    def __repr__(self):
        return ("<%s: counts={%s} length=%d>" %
                (self.__class__, tuple([(k, v) for k, v in self.counts.items()]), self.fp_length))

    def __eq__(self, other):
        if isinstance(other, Fingerprint):
            return (self.__class__ == other.__class__ and
                    self.counts == other.counts and
                    self.fp_length == other.fp_length and
                    np.all(np.in1d(self.indices, other.indices, assume_unique=True)))
        return False
