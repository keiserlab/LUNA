import math
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from luna.mol.entry import MolFileEntry


def generate_residue_matrix(interactions_mngrs, by_interaction=True):
    """Generate a matrix to count interactions per residue.

    Parameters
    ----------
    interactions_mngrs : iterable of \
            :class:`~luna.interaction.calc.InteractionsManager`
        A sequence of :class:`~luna.interaction.calc.InteractionsManager`
        objects from where interactions will be recovered.
    by_interaction : bool
        If True (the default), count the number of each interaction type
        per residue. Otherwise, count the overall number of interactions
        per residue.

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
            if (not inter.src_grp.has_target()
                    and not inter.trgt_grp.has_target()):
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

            comp1 = ";".join(["%s/%s/%d%s" % (r.parent.id,
                                              r.resname,
                                              r.id[1],
                                              r.id[2].strip()) for r in comp1])
            comp2 = ";".join(["%s/%s/%d%s" % (r.parent.id,
                                              r.resname,
                                              r.id[1],
                                              r.id[2].strip()) for r in comp2])

            entry_id = (entry.mol_id if isinstance(entry, MolFileEntry)
                        else entry.to_string())

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
        return pd.pivot_table(df,
                              index=['entry', 'interaction'],
                              columns='residues',
                              values='frequency',
                              fill_value=0)
    else:
        return pd.pivot_table(df,
                              index='entry',
                              columns='residues',
                              values='frequency',
                              fill_value=0)


def heatmap(data_df,
            figsize=None,
            cmap="Blues",
            heatmap_kw=None,
            gridspec_kw=None):
    """ Plot a residue matrix as a color-encoded matrix.

    Parameters
    ----------
    data_df : :class:`pandas.DataFrame`
        A residue matrix produced
        with :func:`~luna.analysis.residues.generate_residue_matrix`.
    figsize : tuple, optional
        Size (width, height) of a figure in inches.
    cmap : str, iterable of str
        The mapping from data values to color space.
        The default value is 'Blues'.
    heatmap_kw : dict, optional
        Keyword arguments for :func:`seaborn.heatmap`.
    gridspec_kw : dict, optional
        Keyword arguments for :class:`matplotlib.gridspec.GridSpec`.
        Used only if the residue matrix (``data_df``) contains interactions.

    Returns
    -------
     : :class:`matplotlib.axes.Axes` or :class:`numpy.ndarray` \
            of :class:`matplotlib.axes.Axes`

    """
    data_df = data_df.reset_index()

    heatmap_kw = heatmap_kw or {}
    gridspec_kw = gridspec_kw or {}

    interactions = None
    if "interaction" in data_df.columns:
        interactions = sorted(data_df["interaction"].unique())
        max_value = data_df[data_df.columns[2:]].max().max()
    else:
        max_value = data_df[data_df.columns[1:]].max().max()

    if not interactions:
        data_df.set_index('entry', inplace=True)

        fig = plt.figure(figsize=figsize)
        ax = sns.heatmap(data_df, cmap=cmap,
                         vmax=max_value, vmin=0, **heatmap_kw)
        ax.set_xlabel("")
        ax.set_ylabel("")
        return ax
    else:
        ncols = 3
        if "ncols" in gridspec_kw:
            ncols = gridspec_kw["ncols"]
            del gridspec_kw["ncols"]
        nrows = math.ceil(len(interactions) / ncols)

        fig, axs = plt.subplots(nrows, ncols, 
                                figsize=figsize, gridspec_kw=gridspec_kw)

        row, col = 0, 0
        for i, interaction in enumerate(interactions):
            df = data_df[data_df["interaction"] == interaction].copy()
            df.drop(columns="interaction", inplace=True)
            df.set_index('entry', inplace=True)

            g = sns.heatmap(df, cmap=cmap, vmax=max_value, vmin=0,
                            ax=axs[row][col], **heatmap_kw)

            g.set_title(interaction)
            g.set_xlabel("")
            g.set_ylabel("")

            col += 1
            if col == ncols:
                row += 1
                col = 0

        if len(interactions) < nrows * ncols:
            diff = (nrows * ncols) - len(interactions)
            for i in range(1, diff + 1):
                axs[-1][-1 * i].axis('off')

        return axs
