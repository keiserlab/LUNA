import numpy as np
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, SparseBitVect
from scipy.sparse import issparse, csr_matrix
from collections import defaultdict
from rdkit import DataStructs

from luna.util.exceptions import (BitsValueError, InvalidFingerprintType,
                                  IllegalArgumentError, FingerprintCountsError)
from luna.version import __version__


import logging

logger = logging.getLogger()

DEFAULT_FP_LENGTH = 2**32
DEFAULT_FOLDED_FP_LENGTH = 4096
DEFAULT_FP_DTYPE = np.int32


class Fingerprint:
    """A fingerprint that stores indices of "on" bits.

    Parameters
    ----------
    indices : array_like of int
        Indices of "on" bits.
    fp_length : int
        The fingerprint length (total number of bits).
        The default value is :math:`2^{32}`.
    unfolded_fp : `Fingerprint` or None
        The unfolded version of this fingerprint.
        If None, this fingerprint may have not been folded yet.
    unfolding_map : dict, optional
        A mapping between current indices and indices from the unfolded version
        of this fingerprint what makes it possible to trace folded bits back to
        the original shells (features).
    props: dict, optional
        Custom properties of the fingerprint, consisting of a string keyword
        and some value. It can be used, for instance, to save the ligand name
        and parameters used to generate shells (IFP features).
    """

    def __init__(self,
                 indices,
                 fp_length=DEFAULT_FP_LENGTH,
                 unfolded_fp=None,
                 unfolding_map=None,
                 props=None):

        indices = np.asarray(indices, dtype=np.long)

        if np.any(np.logical_or(indices < 0, indices >= fp_length)):
            error_msg = "Provided indices are in a different bit scale."
            logger.exception(error_msg)
            raise BitsValueError(error_msg)

        self._indices = np.unique(indices)
        self._fp_length = fp_length
        self._unfolded_fp = unfolded_fp
        self._unfolding_map = unfolding_map or {}
        self._props = props or {}

        self.version = __version__

    @classmethod
    def from_indices(cls, indices, fp_length=DEFAULT_FP_LENGTH, **kwargs):
        """Initialize from an array of indices.

        Parameters
        ----------
        indices : array_like of int
            Indices of "on" bits.
        fp_length : int
            The fingerprint length (total number of bits).
            The default value is :math:`2^{32}`.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `Fingerprint`

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import Fingerprint
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> on_bits = 8
        >>> fp_length = 32
        >>> indices = np.random.randint(0, fp_length, on_bits)
        >>> print(indices)
        [12 15 21  0  3 27  3  7]
        >>> fp = Fingerprint.from_indices(indices, fp_length=fp_length)
        >>> print(fp.indices)
        [ 0  3  7 12 15 21 27]
        >>> print(fp.to_vector(compressed=False))
        [1 0 0 1 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0]
        """
        return cls(indices, fp_length, **kwargs)

    @classmethod
    def from_vector(cls, vector, fp_length=None, **kwargs):
        """Initialize from a vector.

        Parameters
        ----------
        vector : :class:`numpy.ndarray` or `scipy.sparse.csr_matrix`
            Array of bits.
        fp_length : int, optional
            The fingerprint length (total number of bits).
            If not provided, the fingerprint length will be defined based on
            the ``vector`` shape.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `Fingerprint`

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import Fingerprint
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> fp_length = 32
        >>> vector = np.random.choice([0, 1], size=(fp_length,), p=[0.8, 0.2])
        >>> print(vector)
        [0 0 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0 1 0 1 1 0 0 0 0 0 0 1 0 0 0 0]
        >>> fp = Fingerprint.from_vector(vector)
        >>> print(fp.indices)
        [ 7  8 13 17 19 20 27]
        >>> print(fp.fp_length)
        32
        """
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
        """Initialize from a bit string (e.g. '0010100110').

        Parameters
        ----------
        bit_string : str
            String of 0s and 1s.
        fp_length : int, optional
            The fingerprint length (total number of bits).
            If not provided, the fingerprint length will be defined based on
            the string length.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `Fingerprint`

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import Fingerprint
        >>> fp = Fingerprint.from_bit_string("0010100110000010")
        >>> print(fp.indices)
        [ 2  4  7  8 14]
        >>> print(fp.fp_length)
        16
        """
        indices = [i for i, char in enumerate(bit_string) if char != '0']
        if fp_length is None:
            fp_length = len(bit_string)

        return cls.from_indices(indices, fp_length, **kwargs)

    @classmethod
    def from_rdkit(cls, rdkit_fp, **kwargs):
        """Initialize from an RDKit fingerprint.

        Parameters
        ----------
        rdkit_fp : :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` or\
                        :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`
            An existing RDKit fingerprint.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `Fingerprint`
        """
        if (not (isinstance(rdkit_fp, ExplicitBitVect)
                 or isinstance(rdkit_fp, SparseBitVect))):
            error_msg = ("Invalid fingerprint type. RDKit only accepts a "
                         "SparseBitVect or ExplicitBitVect object.")
            logger.exception(error_msg)
            raise TypeError(error_msg)

        fp_length = rdkit_fp.GetNumBits()
        indices = np.asarray(rdkit_fp.GetOnBits(), dtype=np.long)

        return cls.from_indices(indices, fp_length, **kwargs)

    @classmethod
    def from_fingerprint(cls, fp, **kwargs):
        """Initialize from an existing fingerprint.

        Parameters
        ----------
        fp : `Fingerprint`
            An existing fingerprint.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `Fingerprint`
        """
        if not isinstance(fp, Fingerprint):
            error_msg = ("Informed fingerprint is not an instance of %s."
                         % (cls.__class__))
            logger.exception(error_msg)
            raise InvalidFingerprintType(error_msg)

        unfolded_fp = (fp.__class__.from_fingerprint(fp.unfolded_fp)
                       if fp.unfolded_fp is not None else None)
        unfolding_map = dict(fp.unfolding_map)
        props = dict(fp.props)

        return cls.from_indices(fp.indices,
                                fp.fp_length,
                                unfolded_fp=unfolded_fp,
                                unfolding_map=unfolding_map,
                                props=props)

    @property
    def indices(self):
        """array_like of int, read-only: Indices of "on" bits."""
        return self._indices

    @property
    def bit_count(self):
        """int, read-only: Number of "on" bits."""
        return self.indices.shape[0]

    @property
    def density(self):
        """float, read-only: Proportion of "on" bits in fingerprint."""
        return self.bit_count / self.fp_length

    @property
    def counts(self):
        """dict, read-only: Mapping between each index in ``indices`` to the \
        number of counts, which is always 1 for bit fingerprints."""
        return dict([(k, 1) for k in self.indices])

    @property
    def fp_length(self):
        """int, read-only: The fingerprint length (total number of bits)."""
        return self._fp_length

    @property
    def unfolded_fp(self):
        """`Fingerprint` or None, read-only: The unfolded version of this \
        fingerprint. If None, this fingerprint may have not been folded yet."""
        if self._unfolded_fp is None:
            logger.warning("This fingerprint was not previously folded.")
            return None
        return self._unfolded_fp

    @property
    def unfolded_indices(self):
        """array_like of int, read-only: Indices of "on" bits in the unfolded \
        fingerprint."""
        if self._unfolding_map is None:
            logger.warning("This fingerprint was not previously folded.")
            return None
        return self.unfolded_fp.indices

    @property
    def unfolding_map(self):
        """dict, read-only: The mapping between current indices and indices \
        from the unfolded version of this fingerprint what makes it possible \
        to trace folded bits back to the original shells (features)."""
        if self._unfolding_map is None:
            logger.warning("This fingerprint was not previously folded.")
            return None
        return self._unfolding_map

    @property
    def props(self):
        """dict, read-only: The custom properties of the fingerprint."""
        return self._props

    @property
    def name(self):
        """str: The property 'name'. If it was not provided, then return an \
        empty string."""
        return self.props.get("name", "")

    @name.setter
    def name(self, name):
        self.props["name"] = str(name)

    @property
    def num_levels(self):
        """int: The property `num_levels` used to generate this fingerprint
        (see :class:`~luna.interaction.fp.shell.ShellGenerator`).
        If it was not provided, then return None."""
        return self.props.get("num_levels", None)

    @num_levels.setter
    def num_levels(self, num_levels):
        self.props["num_levels"] = str(num_levels)

    @property
    def radius_step(self):
        """float: The property 'radius_step' used to generate this
        fingerprint (see \
        :class:`~luna.interaction.fp.shell.ShellGenerator`). \
        If it was not provided, then return None."""
        return self.props.get("radius_step", None)

    @radius_step.setter
    def radius_step(self, radius_step):
        self.props["radius_step"] = str(radius_step)

    @property
    def num_shells(self):
        """int: The property `num_shells`
        (see :class:`~luna.interaction.fp.shell.ShellGenerator`).
        If it was not provided, then return None."""
        return self.props.get("num_shells", None)

    @num_shells.setter
    def num_shells(self, num_shells):
        self.props["num_shells"] = str(num_shells)

    def get_prop(self, key):
        """Get value of the property ``key``. If not set, raise KeyError."""
        try:
            return self.props[key]
        except KeyError:
            logger.warning("Key '%s' does not exist." % key)
            return None

    def set_prop(self, key, value):
        """Set value to the property ``key``."""
        self.props[key] = value

    def get_num_bits(self):
        """Get the fingerprint length (total number of bits)."""
        return self.fp_length

    def get_num_on_bits(self):
        """Get the number of "on" bits."""
        return self.bit_count

    def get_num_off_bits(self):
        """Get the number of "off" bits."""
        return self.get_num_bits() - self.get_num_on_bits()

    def get_bit(self, index):
        """Get the bit/count value at index ``index``.

        Raises
        ------
        BitsValueError
            If the provided index is in a different bit scale.
        """
        if index in self.counts:
            return self.counts[index]
        elif index >= 0 and index < self.fp_length:
            return 0
        else:
            error_msg = "The provided index is in a different bit scale."
            logger.exception(error_msg)
            raise BitsValueError(error_msg)

    def get_on_bits(self):
        """Get "on" bits.

        Returns
        -------
         : :class:`numpy.ndarray`
        """
        return np.array([k for (k, v) in self.counts.items() if v > 0])

    def to_vector(self, compressed=True, dtype=DEFAULT_FP_DTYPE):
        """Convert this fingerprint to a vector of bits/counts.

        .. warning::
            This function may raise a `MemoryError` exception when using
            huge indices vectors. If you found this issue, you may want
            to try a different data type or apply a folding operation
            before calling `to_vector`.

        Parameters
        -------
        compressed : bool
            If True, build a compressed sparse matrix
            (scipy.sparse.csr_matrix).
        dtype : data-type
            The default value is np.int32.

        Returns
        -------
         : :class:`numpy.ndarray` or `scipy.sparse.csr_matrix`
            Vector of bits/counts.
            Return a compressed sparse matrix (`scipy.sparse.csr_matrix`)
            if ``compressed`` is True. Otherwise, return a Numpy array
            (:class:`numpy.ndarray`)

        Raises
        ------
        BitsValueError
            If some of the fingerprint indices are greater than the
            fingerprint length.
        MemoryError
            If the operation ran out of memory.
        """
        data = [self.counts[i] for i in self.indices]
        if compressed:
            try:
                row = np.zeros(self.bit_count)
                col = self.indices
                vector = csr_matrix((data, (row, col)),
                                    shape=(1, self.fp_length),
                                    dtype=dtype)
            except ValueError as e:
                logger.exception(e)
                raise BitsValueError("Sparse matrix construction failed. "
                                     "Invalid indices or input data.")
        else:
            try:
                # This function is causing a MemoryError exception
                # when using a 2**32 vector.
                vector = np.zeros(self.fp_length, dtype=dtype)
            except MemoryError as e:
                logger.exception(e)
                raise MemoryError("Huge indices vector detected. "
                                  "An operation ran out of memory. "
                                  "Use a different data type or apply a "
                                  "folding operation.")

            try:
                vector[self.indices] = data
            except IndexError as e:
                logger.exception(e)
                raise BitsValueError("Some of the provided indices are "
                                     "greater than the fingerprint length.")

        return vector

    def to_bit_vector(self, compressed=True):
        """Convert this fingerprint to a vector of bits.

        .. warning::
            This function may raise a `MemoryError` exception when using huge
            indices vectors. If you found this issue, you may want to try a
            different data type or apply a folding operation before calling
            `to_bit_vector`.

        Parameters
        -------
        compressed : bool
            If True, build a compressed sparse matrix
            (scipy.sparse.csr_matrix).

        Returns
        -------
         : :class:`numpy.ndarray` or `scipy.sparse.csr_matrix`
            Vector of bits/counts.
            Return a compressed sparse matrix (`scipy.sparse.csr_matrix`)
            if ``compressed`` is True. Otherwise, return a Numpy array
            (:class:`numpy.ndarray`)

        Raises
        ------
        BitsValueError
            If some of the fingerprint indices are greater than the
            fingerprint length.
        MemoryError
            If the operation ran out of memory.
        """
        return self.to_vector(compressed=compressed,
                              dtype=np.bool_).astype(np.int8)

    def to_bit_string(self):
        """Convert this fingerprint to a string of bits.

        .. warning::
            This function may raise a `MemoryError` exception when using huge
            indices vectors. If you found this issue, you may want to try a
            different data type or apply a folding operation before calling
            `to_bit_string`.

        Returns
        -------
         : str

        Raises
        ------
        MemoryError
            If the operation ran out of memory.
        """
        try:
            # This function is causing a MemoryError exception
            # when using a 2**32 vector.
            bit_vector = self.to_bit_vector(compressed=False).astype(np.int8)
            return "".join(map(str, bit_vector))
        except MemoryError as e:
            logger.exception(e)
            raise MemoryError("Huge indices vector detected. An operation ran "
                              "out of memory. Use a different data type or "
                              "apply a folding operation.")

    def to_rdkit(self, rdkit_fp_cls=None):
        """Convert this fingerprint to an RDKit fingerprint.

        .. note::
            If the fingerprint length exceeds the maximum RDKit fingerprint
            length (:math:`2^{31} - 1`), this fingerprint will be folded to
            length :math:`2^{31} - 1` before conversion.

        Returns
        -------
         : :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` or \
                :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`
            If ``fp_length`` is less than :math:`1e5`,
            :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` is used.
            Otherwise, :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`
            is used.
        """
        if rdkit_fp_cls is None:
            # Classes to store explicit bit vectors: ExplicitBitVect or
            #    SparseBitVect.
            # ExplicitBitVect is most useful for situations where the size of
            # the vector is relatively small (tens of thousands or smaller).
            # For larger vectors, use the _SparseBitVect_ class instead.
            if self.fp_length < 1e5:
                rdkit_fp_cls = ExplicitBitVect
            else:
                rdkit_fp_cls = SparseBitVect

        # RDKit data structure defines fingerprints as a std:set composed of
        # ints (signed int). Since we always have values higher than 0 and
        # since the data structure contains only signed ints, then the max
        # length for a RDKit fingerprint is 2^31 - 1.
        # C signed int (32 bit) ranges: [-2^31, 2^31-1].
        max_rdkit_fp_length = 2**31 - 1
        fp_length = self.fp_length
        if max_rdkit_fp_length < fp_length:
            logger.warning("The current fingerprint will be folded as its "
                           "size is higher than the maximum size accepted "
                           "by RDKit, which is 2**31 - 1.")
            fp_length = max_rdkit_fp_length
        indices = self.indices % max_rdkit_fp_length

        rdkit_fp = rdkit_fp_cls(fp_length)
        rdkit_fp.SetBitsFromList(indices.tolist())
        return rdkit_fp

    def fold(self, new_length=DEFAULT_FOLDED_FP_LENGTH):
        """Fold this fingerprint to size ``new_length``.

        Parameters
        ----------
        new_length : int
            Length of the new fingerprint, ideally multiple of 2.
            The default value is 4096.

        Returns
        -------
         : `Fingerprint`
            Folded `Fingerprint`.

        Raises
        ------
        BitsValueError
            If the new fingerprint length is not a multiple of 2 or is greater
            than the existing fingerprint length.

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import Fingerprint
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> on_bits = 8
        >>> fp_length = 32
        >>> indices = np.random.randint(0, fp_length, on_bits)
        >>> print(indices)
        [12 15 21  0  3 27  3  7]
        >>> fp = Fingerprint.from_indices(indices, fp_length=fp_length)
        >>> print(fp.indices)
        [ 0  3  7 12 15 21 27]
        >>> print(fp.to_vector(compressed=False))
        [1 0 0 1 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0]
        >>> folded_fp = fp.fold(8)
        >>> print(folded_fp.indices)
        [0 3 4 5 7]
        >>> print(folded_fp.to_vector(compressed=False))
        [1 0 0 1 1 1 0 1]
        """
        if new_length > self.fp_length:
            error_msg = ("The new fingerprint length must be smaller than the "
                         "existing fingerprint length.")
            logger.exception(error_msg)
            raise BitsValueError(error_msg)

        if not np.log2(self.fp_length / new_length).is_integer():
            error_msg = ("It is not possible to fold the current fingerprint "
                         "into the informed new length. The current length "
                         "divided by the new one is not a power of 2 number.")
            logger.exception(error_msg)
            raise BitsValueError(error_msg)

        folded_indices = self.indices % new_length

        unfolding_map = defaultdict(set)
        for k, v in sorted(zip(folded_indices, self.indices)):
            unfolding_map[k].add(v)

        props = dict(self.props)
        if "fp_length" in props:
            props["fp_length"] = new_length
        new_fp = self.__class__(indices=folded_indices,
                                fp_length=new_length,
                                unfolded_fp=self,
                                unfolding_map=unfolding_map,
                                props=props)

        return new_fp

    def unfold(self):
        """Unfold this fingerprint and return its parent fingerprint.

        Returns
        -------
         : `Fingerprint`
        """
        return self.unfolded_fp

    def union(self, other):
        """Return the union of indices of two fingerprints.

        Returns
        -------
         : :class:`numpy.ndarray`

        Raises
        ------
        InvalidFingerprintType
            If the informed fingerprint is not an instance of `Fingerprint`.
        BitsValueError
            If the fingerprints have different lengths.
        """
        if not isinstance(other, Fingerprint):
            error_msg = ("The informed fingerprint is not an instance of %s."
                         % (other.__class__))
            logger.exception(error_msg)
            raise InvalidFingerprintType(error_msg)

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.union1d(self.indices, other.indices)

    def intersection(self, other):
        """Return the intersection between indices of two fingerprints.

        Returns
        -------
         : :class:`numpy.ndarray`

        Raises
        ------
        InvalidFingerprintType
            If the informed fingerprint is not an instance of `Fingerprint`.
        BitsValueError
            If the fingerprints have different lengths.
        """
        if not isinstance(other, Fingerprint):
            error_msg = ("Informed fingerprint is not an instance of %s."
                         % (other.__class__))
            logger.exception(error_msg)
            raise InvalidFingerprintType(error_msg)

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.intersect1d(self.indices, other.indices, assume_unique=True)

    def difference(self, other):
        """Return indices in this fingerprint but not in ``other``.

        Returns
        -------
         : :class:`numpy.ndarray`

        Raises
        ------
        InvalidFingerprintType
            If the informed fingerprint is not an instance of `Fingerprint`.
        BitsValueError
            If the fingerprints have different lengths.
        """
        if not isinstance(other, Fingerprint):
            error_msg = ("Informed fingerprint is not an instance of %s."
                         % (other.__class__))
            logger.exception(error_msg)
            raise InvalidFingerprintType(error_msg)

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.setdiff1d(self.indices, other.indices, assume_unique=True)

    def symmetric_difference(self, other):
        """Return indices in either this fingerprint or ``other`` but not both.

        Returns
        -------
         : :class:`numpy.ndarray`

        Raises
        ------
        InvalidFingerprintType
            If the informed fingerprint is not an instance of `Fingerprint`.
        BitsValueError
            If the fingerprints have different lengths.
        """
        if not isinstance(other, Fingerprint):
            error_msg = ("Informed fingerprint is not an instance of %s."
                         % (other.__class__))
            logger.exception(error_msg)
            raise InvalidFingerprintType(error_msg)

        if self.fp_length != other.fp_length:
            raise BitsValueError("Fingerprints are in a different bit scale")

        return np.setxor1d(self.indices, other.indices, assume_unique=True)

    def calc_similarity(self, other):
        """Calculates the Tanimoto similarity between this fingeprint
        and ``other``.

        Returns
        -------
         : float

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import Fingerprint
        >>> fp1 = Fingerprint.from_bit_string("0010101110000010")
        >>> fp2 = Fingerprint.from_bit_string("1010100110010010")
        >>> print(fp1.calc_similarity(fp2))
        0.625
        """
        return DataStructs.FingerprintSimilarity(self.to_rdkit(),
                                                 other.to_rdkit())

    def __repr__(self):
        return ("<%s: indices=%s length=%d>" %
                (self.__class__,
                 repr(self.indices).replace('\n', '').replace(' ', ''),
                 self.fp_length))

    def __eq__(self, other):
        if isinstance(other, Fingerprint):
            return (self.__class__ == other.__class__
                    and self.fp_length == other.fp_length
                    and np.all(np.in1d(self.indices,
                                       other.indices,
                                       assume_unique=True)))
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

    """A fingerprint that stores the number of occurrences of each index.

    Parameters
    ----------
    indices : array_like of int, optional
        Indices of "on" bits. It is optional if ``counts`` is provided.
    counts : dict, optional
        Mapping between each index in ``indices`` to the number of counts.
        If not provided, the default count value of 1 will be used instead.
    fp_length : int
        The fingerprint length (total number of bits).
        The default value is :math:`2^{32}`.
    unfolded_fp : `Fingerprint` or None
        The unfolded version of this fingerprint.
        If None, this fingerprint may have not been folded yet.
    unfolding_map : dict, optional
        A mapping between current indices and indices from the unfolded version
        of this fingerprint what makes it possible to trace folded bits back to
        the original shells (features).
    props: dict, optional
        Custom properties of the fingerprint, consisting of a string keyword
        and some value. It can be used, for instance, to save the ligand name
        and parameters used to generate shells (IFP features).
    """

    def __init__(self, indices=None, counts=None, fp_length=DEFAULT_FP_LENGTH,
                 unfolded_fp=None, unfolding_map=None, props=None):

        if indices is None and counts is None:
            logger.exception("Indices or counts must be provided.")
            raise IllegalArgumentError("Indices or counts must be provided.")

        if indices is not None:
            indices = np.asarray(indices, dtype=np.long)

            if np.any(np.logical_or(indices < 0, indices >= fp_length)):
                error_msg = "Provided indices are in a different bit scale."
                logger.exception(error_msg)
                raise BitsValueError(error_msg)

            if counts is None:
                indices, counts = np.unique(indices, return_counts=True)
                counts = dict(zip(indices, counts))
            else:
                indices = np.unique(indices)
                if not np.all([x in indices for x in counts]):
                    error_msg = ("At least one index from 'counts' is not in "
                                 "'indices'.")
                    logger.exception(error_msg)
                    raise FingerprintCountsError(error_msg)

                if len(set(indices).symmetric_difference(counts)) > 0:
                    error_msg = ("At least one index in 'indices' is not "
                                 "in 'counts'.")
                    logger.exception(error_msg)
                    raise FingerprintCountsError(error_msg)
        else:
            indices = np.asarray(sorted(counts.keys()), dtype=np.long)

            if np.any(np.logical_or(indices < 0, indices >= fp_length)):
                error_msg = ("Provided indices are in a different bit scale.")
                logger.exception(error_msg)
                raise BitsValueError(error_msg)

        self._counts = counts
        super().__init__(indices, fp_length, unfolded_fp, unfolding_map, props)

    @classmethod
    def from_indices(cls,
                     indices=None,
                     counts=None,
                     fp_length=DEFAULT_FP_LENGTH,
                     **kwargs):
        """Initialize from an array of indices.

        Parameters
        ----------
        indices : array_like of int, optional
            Indices of "on" bits. It is optional if ``counts`` is provided.
        counts : dict, optional
            Mapping between each index in ``indices`` to the number of counts.
            If not provided, the default count value of 1 will be used instead.
        fp_length : int
            The fingerprint length (total number of bits).
            The default value is :math:`2^{32}`.
        **kwargs : dict, optional
            Extra arguments to `CountFingerprint`. Refer to the documentation
            for a list of all possible arguments.

        Returns
        -------
         : `CountFingerprint`

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import CountFingerprint
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> on_bits = 8
        >>> fp_length = 32
        >>> indices, counts = np.unique(np.random.randint(0, \
fp_length, on_bits), return_counts=True)
        >>> counts = dict(zip(indices, counts))
        >>> print(counts)
        {0: 1, 3: 2, 7: 1, 12: 1, 15: 1, 21: 1, 27: 1}
        >>> fp = CountFingerprint.from_indices(indices, counts=counts, \
fp_length=fp_length)
        >>> print(fp.indices)
        [ 0  3  7 12 15 21 27]
        >>> print(fp.to_vector(compressed=False))
        [1 0 0 2 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0]
        """
        return cls(indices=indices, 
                   counts=counts,
                   fp_length=fp_length,
                   **kwargs)

    @classmethod
    def from_counts(cls, counts, fp_length=DEFAULT_FP_LENGTH, **kwargs):
        """Initialize from a counting map.

        Parameters
        ----------
        counts : dict
            Mapping between each index in ``indices`` to the number of counts.
        fp_length : int
            The fingerprint length (total number of bits).
            The default value is :math:`2^{32}`.
        **kwargs : dict, optional
            Extra arguments to `CountFingerprint`. Refer to the documentation
            for a list of all possible arguments.

        Returns
        -------
         : `CountFingerprint`

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import CountFingerprint
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> on_bits = 8
        >>> fp_length = 32
        >>> counts = dict(zip(*np.unique(np.random.randint(0, fp_length, \
on_bits),
        ...                              return_counts=True)))
        >>> print(counts)
        {0: 1, 3: 2, 7: 1, 12: 1, 15: 1, 21: 1, 27: 1}
        >>> fp = CountFingerprint.from_counts(counts=counts,
        ...                                   fp_length=fp_length)
        >>> print(fp.indices)
        [ 0  3  7 12 15 21 27]
        >>> print(fp.to_vector(compressed=False))
        1 0 0 2 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0]
        """
        return cls(counts=counts, fp_length=fp_length, **kwargs)

    @classmethod
    def from_bit_string(cls,
                        bit_string,
                        counts=None,
                        fp_length=None,
                        **kwargs):
        """Initialize from a bit string (e.g. '0010100110').

        Parameters
        ----------
        bit_string : str
            String of 0s and 1s.
        counts : dict, optional
            Mapping between each index in ``indices`` to the number of counts.
            If not provided, the default count value of 1 will be used instead.
        fp_length : int, optional
            The fingerprint length (total number of bits).
            If not provided, the fingerprint length will be defined based on
            the string length.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `CountFingerprint`

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import CountFingerprint
        >>> fp = CountFingerprint.from_bit_string("0010100110000010",
        ...                                       counts={2: 5, 4: 1, 7: 3, \
8: 1, 14: 2})
        >>> print(fp.indices)
        [ 2  4  7  8 14]
        >>> print(fp.counts)
        {2: 5, 4: 1, 7: 3, 8: 1, 14: 2}
        """
        indices = [i for i, char in enumerate(bit_string) if char != '0']
        if fp_length is None:
            fp_length = len(bit_string)

        return cls.from_indices(indices, counts, fp_length, **kwargs)

    @classmethod
    def from_vector(cls, vector, fp_length=None, **kwargs):
        """Initialize from a vector.

        Parameters
        ----------
        vector : :class:`numpy.ndarray` or `scipy.sparse.csr_matrix`
            Array of counts.
        fp_length : int, optional
            The fingerprint length (total number of bits).
            If not provided, the fingerprint length will be defined based on
            the ``vector`` shape.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `CountFingerprint`

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import CountFingerprint
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> fp_length = 32
        >>> vector = np.random.choice(5, size=(fp_length,),
        ...                           p=[0.76, 0.1, 0.1, 0.02, 0.02])
        >>> print(vector)
        [0 0 0 0 2 3 0 1 0 0 2 0 0 0 1 1 2 3 1 0 1 0 0 0 2 0 0 0 1 0 0 0]
        >>> fp = CountFingerprint.from_vector(vector)
        >>> print(fp.indices)
        [ 4  5  7 10 14 15 16 17 18 20 24 28]
        >>> print(fp.counts)
        {4: 2, 5: 3, 7: 1, 10: 2, 14: 1, 15: 1, 16: 2, 17: 3, 18: 1, 20: 1, \
24: 2, 28: 1}
        """
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
        """Initialize from an existing fingerprint.

        Parameters
        ----------
        fp : `Fingerprint`
            An existing fingerprint.
        **kwargs : dict, optional
            Extra arguments to `Fingerprint`. Refer to the documentation for a
            list of all possible arguments.

        Returns
        -------
         : `CountFingerprint`
        """
        if not isinstance(fp, Fingerprint):
            error_msg = ("Informed fingerprint is not an instance of %s."
                         % (cls.__class__))
            logger.exception(error_msg)
            raise InvalidFingerprintType(error_msg)

        counts = dict([(i, c) for i, c in fp.counts.items() if c > 0])
        unfolded_fp = (fp.__class__.from_fingerprint(fp.unfolded_fp)
                       if fp.unfolded_fp is not None else None)
        unfolding_map = dict(fp.unfolding_map)
        props = dict(fp.props)

        new_fp = cls.from_counts(counts, fp.fp_length, unfolded_fp=unfolded_fp,
                                 unfolding_map=unfolding_map, props=props)

        return new_fp

    @property
    def counts(self):
        """dict, read-only: Mapping between each index in ``indices`` \
        to the number of counts."""
        return self._counts

    def get_count(self, index):
        """Get the count value at index ``index``. Return 0 if index \
        is not in ``counts``."""
        return self.counts.get(index, 0)

    def fold(self, new_length=DEFAULT_FOLDED_FP_LENGTH):
        """Fold this fingerprint to size ``new_length``.

        Parameters
        ----------
        new_length : int
            Length of the new fingerprint, ideally multiple of 2.
            The default value is 4096.

        Returns
        -------
         : `Fingerprint`
            Folded `Fingerprint`.

        Raises
        ------
        BitsValueError
            If the new fingerprint length is not a multiple of 2 or is greater
            than the existing fingerprint length.

        Examples
        --------
        >>> from luna.interaction.fp.fingerprint import CountFingerprint
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> on_bits = 8
        >>> fp_length = 32
        >>> indices, counts = np.unique(np.random.randint(0, fp_length, \
on_bits),
        ...                             return_counts=True)
        >>> counts = dict(zip(indices, counts))
        >>> print(counts)
        {0: 1, 3: 2, 7: 1, 12: 1, 15: 1, 21: 1, 27: 1}
        >>> fp = CountFingerprint.from_indices(indices, counts=counts,
        ...                                    fp_length=fp_length)
        >>> print(fp.indices)
        [ 0  3  7 12 15 21 27]
        >>> print(fp.to_vector(compressed=False))
        [1 0 0 2 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0]
        >>> folded_fp = fp.fold(8)
        >>> print(folded_fp.indices)
        [0 3 4 5 7]
        >>> print(folded_fp.to_vector(compressed=False))
        [1 0 0 3 1 1 0 2]
        """
        new_fp = super().fold(new_length)

        new_fp._counts = dict([(folded_idx, sum([self.get_count(x)
                                                 for x in unfolded_set]))
                               for folded_idx, unfolded_set
                               in new_fp.unfolding_map.items()])

        return new_fp

    def __repr__(self):
        return ("<%s: counts={%s} length=%d>" %
                (self.__class__,
                 tuple([(k, v)
                        for k, v in self.counts.items()]),
                 self.fp_length))

    def __eq__(self, other):
        if isinstance(other, Fingerprint):
            return (self.__class__ == other.__class__
                    and self.counts == other.counts
                    and self.fp_length == other.fp_length
                    and np.all(np.in1d(self.indices,
                                       other.indices,
                                       assume_unique=True)))
        return False
