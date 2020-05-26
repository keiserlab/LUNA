import math

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from collections import defaultdict


def res_freq_heatmap(num_entries, freqs, output_file, keep_only_target=True, ignore_intra_inter=True, freq_cutoff=0.1,
                     num_cols=3, fig_size=(60, 65), format="png"):

    data = defaultdict(list)

    for key, rf_obj in freqs.items():

        freq = rf_obj.freq / num_entries

        if freq >= freq_cutoff:
            # We can validate using only one of the interactions as the residues/ligands have the same names.
            inter = rf_obj.interactions[0]

            # Apply filters.
            if keep_only_target:
                # Continue if no target is in the interaction.
                if not inter.src_grp.has_target() and not inter.trgt_grp.has_target():
                    continue
            if ignore_intra_inter:
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

            data["src_grps"].append(comp1)
            data["trgt_grps"].append(comp2)
            data["interaction"].append(inter.type)
            data["frequency"].append(freq * 100)

    df = pd.DataFrame.from_dict(data)
    df = df.sort_values(['src_grps', 'trgt_grps', "interaction", 'frequency'])

    interactions = df["interaction"].unique()
    interactions.sort()
    num_rows = math.ceil(len(interactions) / num_cols)

    fig = plt.figure(figsize=fig_size)
    gs = gridspec.GridSpec(nrows=num_rows, ncols=num_cols, figure=fig,
                           width_ratios=[1] * num_cols, height_ratios=[3] * num_rows,
                           wspace=0.4, hspace=0.3)

    col = 0
    line = 0
    for i, inter in enumerate(interactions):
        df_filtered = df.query("interaction == '%s'" % inter)
        df_filtered = df_filtered.pivot("trgt_grps", "src_grps", "frequency")

        ax = fig.add_subplot(gs[line, col])

        ax.set_title(inter, fontsize=50, y=1.1)

        sns.heatmap(df_filtered, cmap="Blues", linewidths=2, vmax=100, vmin=0, ax=ax)

        ax.invert_yaxis()
        ax.set_ylabel("")
        ax.set_xlabel("")

        ax.set_xticklabels([])

        ax.tick_params(axis="y", labelsize=40, pad=50, rotation=0)
        ax.tick_params(labelsize=40, pad=70)

        col += 1
        if col == num_cols:
            line += 1
            col = 0

    fig.savefig(output_file, format=format, bbox_inches="tight")
