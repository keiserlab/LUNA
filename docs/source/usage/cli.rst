Command Line Interface
======================

Command line interfaces (CLI) are provided for the most common task:
calculate interactions for a series of protein-ligand complexes.

In the below examples, we assume the LUNA repository is located at $LUNA_REPO and that
the current working directory is $LUNA_REPO/example.


run.py
------

This CLI provides options to calculate protein-ligand interactions using LUNA default methods,
save Pymol sessions to analyze interactions, generate interaction fingerprints, and calculate the
similarity between ligand binding modes using these generated fingerprints.

To see all available options, run:

.. code:: bash

    $ python $LUNA_REPO/luna/run.py --help

    usage: run.py [-h] -p PDB_FILE -e ENTRIES_FILE -l LIGAND_FILE -w WORKING_PATH
                  [--out_ifp] [-L IFP_NUM_LEVELS] [-R IFP_RADIUS_STEP]
                  [-S IFP_LENGTH] [-T {EIFP,HIFP,FIFP}] [-B] [-O IFP_OUTPUT]
                  [--sim_matrix_output SIM_MATRIX_OUTPUT]
                  [--filter_binding_modes BINDING_MODES_FILE] [--out_pse]
                  [--pse_path PSE_PATH] [--overwrite] [--nproc NPROC]

    optional arguments:
      -h, --help            show this help message and exit
      -p PDB_FILE, --prot PDB_FILE
                            the protein PDB file
      -e ENTRIES_FILE, --entries ENTRIES_FILE
                            an input file containing a list of ligand ids to process
      -l LIGAND_FILE, --lig LIGAND_FILE
                            a molecular file containing 1 or more ligands
      -w WORKING_PATH       the path where the project and its results will be saved
      --out_ifp             defines whether it should generate LUNA interaction fingerprints
      -L IFP_NUM_LEVELS     the number of level defines the number of iterations to construct
                            the fingerprint. Default: 2
      -R IFP_RADIUS_STEP    the radius growth rate defines the multiplier to increase the sphere
                            size at each level. Default: 5.73171
      -S IFP_LENGTH         the fingerprint length. Default: 4096
      -T {EIFP,HIFP,FIFP}   the fingerprint type. Default: EIFP
      -B                    defines whether it should use bit fingerprints. The default value
                            is False, which implies that count fingerprints are used instead
      -O IFP_OUTPUT         the fingerprint output file.
                            Default: <WORKING_PATH>/results/fingerprints/ifp.csv
      --sim_matrix_output SIM_MATRIX_OUTPUT
                            the path where the similarity matrix will be saved. If not provided,
                            it won't be generated
      --filter_binding_modes BINDING_MODES_FILE
                            the path of a file containing binding modes to filter
      --out_pse             defines whether it should export interactions to Pymol
      --pse_path PSE_PATH   the path where Pymol sessions (PSE files) will be saved.
                            Default: <WORKING_PATH>/results/pse/
      --overwrite           defines whether it should overwrite an existing project
      -f FORK_PROJECT, --fork_project FORK_PROJECT
                            If provided, copy an existing project to <WORKING_PATH>
      --nproc NPROC         the number of processors to use



Calculate protein-ligand interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following example, we will calculate protein-ligand interactions for docked Dopamine D4 complexes.
The input files are located at $LUNA_REPO/examples/inputs.

.. code:: bash

  $ python $LUNA_REPO/luna/run.py -p inputs/protein.pdb -l inputs/ligands.mol2 -e inputs/entries.txt -w dopamine_results


Generate interaction fingerprints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To generate interaction fingerprints, you need to activate it with the flag ``--out_ifp`` and
modify fingerprint parameters (options ``-L``, ``-R``, ``-S``, ``-T``, ``-O``) as necessary.


.. code:: bash

    $ python $LUNA_REPO/luna/run.py -p inputs/protein.pdb -l inputs/ligands.mol2 -e inputs/entries.txt -w dopamine_results
                                    --out_ifp -L 2 -R 5.73 -S 4096 -T EIFP -O dopamine_results/results/new_fp.csv

.. note::

  Note that if you have an existing project, you can provide its working path (``-w``) and then LUNA will automatically load the
  entire project results, which allows you to generate different fingerprints without reprocessing everything from scratch.


Generate similarity matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~

To compute the Tanimoto similarity between interaction fingerprints (IFPs) and create a similarity matrix, you need to use the flag
``--sim_matrix_output`` in conjunction with ``--out_ifp``.

.. code:: bash

  $ python $LUNA_REPO/luna/run.py -p inputs/protein.pdb -l inputs/ligands.mol2 -e inputs/entries.txt -w dopamine_results
                                  --out_ifp --sim_matrix_output sim_matrix.csv


Visualize interactions on Pymol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To depict interactions as Pymol sessions, you need to activate the flag ``--out_pse``.

.. code:: bash

    $ python $LUNA_REPO/luna/run.py -p inputs/protein.pdb -l inputs/ligands.mol2
                                    -e inputs/entries.txt -w dopamine_results --out_pse


Filter interactions by binding mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To filter interactions based on binding modes, you can use the option ``--filter_binding_modes``.
This option expects a configuration file that defines how interactions should be filtered.
See an example from the configuration file ``$LUNA_REPO/example/inputs/binding_modes.cfg``::

  ; To configurate an interaction type, create a new line and define the interaction: [New interaction].
  ; Then you can define whether or not all interactions must be accepted by setting 'accept_only' to True or False.

  ; If you want to specify binding modes, use the variable 'accept_only', which expects a list of strings \
  in the format: <CHAIN ID>/<COMPOUND NAME>/<COMPOUND NUMBER>/<ATOM>
  ; Wildcards are accepted for the expected fields.
  ; For example, "*/HIS/*/*" represents all histidines' atoms from all chains.
  ;               "A/CBL/*/*" represents all ligands named CBL from chain A.
  ;               "B/HIS/*/N*" represents all histidines' nitrogens from chain B.

  [Hydrogen bond]
  accept_only=["A/LYS/245/*", "*/HIS/*/*"]

  [Hydrophobic]
  accept_all=True

  [Cation-pi]
  accept_only=["*"]
  accept_all=False

  [Weak hydrogen bond]
  accept_all=False
  accept_only=["*/THR/434/O*"]

  [Face-to-edge pi-stacking]
  accept_all=False

  [Aromatic stacking]
  accept_all=True

  [*]
  accept_all=False

.. code:: bash

    $ python $LUNA_REPO/luna/run.py -p inputs/protein.pdb -l inputs/ligands.mol2
                                    -e inputs/entries.txt -w dopamine_results --out_pse
                                    --filter_binding_modes inputs/binding_modes.cfg

.. warning::

  After executing the command above, existing results will be overwritten. If you want to keep the original
  results, you should fork the target project. To do so, use the option ``-f`` or ``--filter_binding_modes``.
  Thus, filterings will only have an effect on the forked project.

  See an example in the next section.


Fork an existing project
~~~~~~~~~~~~~~~~~~~~~~~~

This option allows you to fork an existing project to apply filterings without modifying the original project.
To do so, you should use the option ``-f`` or ``--filter_binding_modes``.

.. code:: bash

    $ python $LUNA_REPO/luna/run.py -p inputs/protein.pdb -l inputs/ligands.mol2 -e inputs/entries.txt
                                    -w filtered_dopamine_results -f dopamine_results
                                    --filter_binding_modes inputs/binding_modes.cfg --out_pse

.. note::

  If you fork a project without the filtering option, it will only create a copy of the original project.
  Later, if you decide to filter interactions, you don’t need to use the fork option again.
  Just use the filtering option directly on the forked directory and it will be overwritten.