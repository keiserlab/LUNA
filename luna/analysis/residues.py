from collections import defaultdict
import pandas as pd
from luna.mol.entry import MolFileEntry


def generate_residue_matrix(interactions_mngrs, by_interaction=True):
    """Generate a matrix to count interactions per residue.

    Parameters
    ----------
    interactions_mngrs : iterable of :class:`~luna.interaction.calc.InteractionsManager`
        A sequence of :class:`~luna.interaction.calc.InteractionsManager` objects
        from where interactions will be recovered.
    by_interaction : bool
        If True (the default), count the number of each interaction type per residue.
        Otherwise, count the overall number of interactions per residue.

    Returns
    -------
     : :class:`pandas.DataFrame`
    """

    data_by_entry = defaultdict(lambda: defaultdict(int))
    residues = set()

    for inter_mngr in interactions_mngrs:
        entry = inter_mngr.entry

        for inter in inter_mngr:
            # Continue if no target is in the interaction.
            if not inter.src_grp.has_target() and not inter.trgt_grp.has_target():
                continue
            # Ignore interactions involving the same compounds.
            if inter.src_grp.compounds == inter.trgt_grp.compounds:
                continue

            if inter.src_grp.has_hetatm():
                comp1 = sorted(inter.src_grp.compounds)
                comp2 = sorted(inter.trgt_grp.compounds)
            elif inter.trgt_grp.has_hetatm():
                comp1 = sorted(inter.trgt_grp.compounds)
                comp2 = sorted(inter.src_grp.compounds)
            else:
                comp1 = sorted(inter.src_grp.compounds)
                comp2 = sorted(inter.trgt_grp.compounds)
                comp1, comp2 = sorted([comp1, comp2])

            comp1 = ";".join(["%s/%s/%d%s" % (r.parent.id, r.resname, r.id[1], r.id[2].strip()) for r in comp1])
            comp2 = ";".join(["%s/%s/%d%s" % (r.parent.id, r.resname, r.id[1], r.id[2].strip()) for r in comp2])

            entry_id = entry.mol_id if isinstance(entry, MolFileEntry) else entry.to_string()
            if by_interaction:
                key = (entry_id, inter.type)
            else:
                key = entry_id

            data_by_entry[key][comp2] += 1

            residues.add(comp2)

    heatmap_data = defaultdict(list)

    if by_interaction:
        entries = set([k[0] for k in data_by_entry.keys()])
        interactions = set([k[1] for k in data_by_entry.keys()])

        for e in entries:
            for i in interactions:
                for res in residues:
                    heatmap_data["entry"].append(e)
                    heatmap_data["interaction"].append(i)
                    heatmap_data["residues"].append(res)
                    heatmap_data["frequency"].append(data_by_entry[(e, i)][res])

    else:
        for key in data_by_entry:
            for res in data_by_entry[key]:
                heatmap_data["entry"].append(key)

                heatmap_data["residues"].append(res)
                heatmap_data["frequency"].append(data_by_entry[key][res])

    df = pd.DataFrame.from_dict(heatmap_data)

    if by_interaction:
        return pd.pivot_table(df, index=['entry', 'interaction'], columns='residues', values='frequency', fill_value=0)
    else:
        return pd.pivot_table(df, index='entry', columns='residues', values='frequency', fill_value=0)
