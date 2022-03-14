Overview of LUNA
================

Introduction
------------

LUNA is an object-oriented Python 3 toolkit for drug design that makes it easy to analyze very large data sets of 3D molecular structures and complexes, and that allows identifying, filtering, and visualizing atomic interactions. LUNA implements several features geared towards the analysis of molecular complexes, such as: a) accepting any molecular complex type, including protein-ligand and protein-protein; b) accepting multiple file formats, including PDB, MOL, and MOL2; c) providing pre- and post-processing functions to control how interactions are calculated or selected based on geometric constraints; and d) providing several functions to summarize, characterize, and visualize molecular interactions in Pymol.

LUNA also implements three hashed interaction fingerprints (IFP): Extended Interaction FingerPrint (EIFP), Functional Interaction FingerPrint (FIFP), and Hybrid Interaction FingerPrint (HIFP) -- inspired by ECFP, FCFP [1]_, and E3FP [2]_. While EIFP encodes explicit atomic substructures, FIFP only encodes more ''coarse-grained" pharmacophoric properties. HIFP adopts a ''hybrid" approach, encoding pharmacophoric properties for atom groups and precise environment information for atoms. All these IFPs are RDKit_-compatible and their features represent interactions and contacts between protein residues, ligand atoms, and water molecules, by detecting their presence or absence (bit FPs), or their frequency (count FPs). Besides, these IFPs encode molecular interactions at different levels of detail and are fully interpretable, providing several functionalities to trace individual bits back to their original atomic substructures in the context of the binding site.

LUNA is developed by the `Keiser Lab`_ at UCSF_ and maintained primarily by Alexandre Fassio.

For a thorough description of LUNA and the three IFPs, please consult the original paper [3]_ and
`paper repository`_ or :ref:`Usage and Examples`.

.. Documentation is hosted by ReadTheDocs_.


Contributing
------------

Development occurs on GitHub_.
Contributions, feature requests, and bug reports are greatly appreciated.
Please consult the `issue tracker`_.


License
-------
LUNA is released under the |license|.

Briefly, this means LUNA can be used in any manner without modification, with proper attribution. However, if the source code is modified for an
application, this modified source must also be released under |license| so that the community may benefit.


Citing LUNA
-----------

To cite LUNA, please reference the original paper [3]_.

.. rubric:: References

.. [1] |rogers2010|
.. [2] |axen2017|
.. [3] |afassio2022|

.. include:: substitutions.rst

.. _GitHub: https://github.com/keiserlab/luna
