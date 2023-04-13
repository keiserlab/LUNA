from collections import defaultdict
import logging


logger = logging.getLogger()


def count_interaction_types(interactions,
                            must_have_target=True,
                            compounds=None,
                            key_map={}):
    """Count the number of each type of interaction in ``interactions``.

    Parameters
    ----------
    interactions : iterable
            An iterable object containing a sequence of interactions
            (``InteractionType``).
    must_have_target : bool
            If True, count only interactions involving the target ligand.
            The default value is True.
    compounds: iterable, optional
            Only count interactions involving the compounds in ``compounds``.
            The default value is None, which implies that all compounds will
            be considered.
    key_map: dictionary, optional
            A dictionary to control which interactions to count and how to
            aggregate them. The keys are the interaction types to be
            considered and the values are the final interaction type, which
            can be used to aggregate interactions. If a value is None, the
            interaction will be ignored.

            For example, to aggregate all covalent interactions, ``key_map``
            could be defined as follows:

               .. code-block:: python

                    key_map = {"Single bond": "Covalent bond",
                               "Double bond": "Covalent bond",
                               "Triple bond": "Covalent bond",
                               "Aromatic bond": "Covalent bond"}

            Now, if Ionic interactions should be ignored, ``key_map`` could
            be defined as follows:

                .. code-block:: python

                    key_map = {"Ionic": None}

    Returns
    -------
    interaction_types_count : dict
        The count of each interaction type.
    """
    if compounds:
        compounds = set(compounds)

    interaction_types_count = defaultdict(int)
    seen_pairs = set()
    for i in interactions:
        is_valid = False
        if not must_have_target and compounds is None:
            is_valid = True
        else:
            if must_have_target:
                if i.src_grp.has_target() or i.trgt_grp.has_target():
                    is_valid = True

            if compounds:
                if (len(i.src_grp.compounds & compounds) > 0
                        or len(i.trgt_grp.compounds & compounds)):
                    is_valid = True

        if is_valid:
            pair_key1 = (i.type, i.src_grp, i.trgt_grp)
            pair_key2 = (i.type, i.trgt_grp, i.src_grp)

            if pair_key1 in seen_pairs or pair_key2 in seen_pairs:
                continue

            seen_pairs.add(pair_key1)
            key = i.type
            if key_map:
                if i.type in key_map:
                    if key_map[i.type] is None:
                        continue
                    key = key_map[i.type]
                else:
                    key = i.type

            interaction_types_count[key] += 1
    return interaction_types_count
