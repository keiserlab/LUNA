Developer notes
===============

.. _Ligands in PDB files:
Ligands in PDB files
--------------------

LUNA works with both `Open Babel`_ and RDKit_. You can switch between both libraries when setting a new project. So, the choice between `Open Babel`_ and RDKit_ depends on the type of molecular file you are working with. If you have ligands in MOL or MOL2 files, for instance, you can use RDKit_ for your convenience as you will be able to apply RDKit_ functions directly on the molecular object stored at :doc:`LUNA entries <../api/luna.mol.entry>`. However, if you frequently work with **ligands** in **PDB files**, you may be especially interested in reading this note.

Although each library has its own difficulties at parsing PDB files, we identified that `Open Babel`_ works better at parsing ligands from PDB files, especially what concerns ligands containing aromatic rings. That happens because PDB files don't contain information about atomic charges, valences, and bond types. Consequently, it may occur that some ligands are incorrectly perceived by `Open Babel`_ and RDKit_.

For this reason, when working with ligands in PDB files, no matter you choose `Open Babel`_ or RDKit_, we use `Open Babel`_ under the hood for parsing, amending, and converting ligands so that both final `Open Babel`_ and RDKit_ objects reflect at most the original and expected molecular structure.

However, in the past, I identified different errors in different `Open Babel`_ versions, which can be seen in this `issue thread <https://github.com/openbabel/openbabel/issues/1925>`_. During that time, I identified that `Open Babel`_ 2.3.2 works better for parsing PDB files. For the most common errors identified in this version, I implemented :doc:`amending solutions <../api/luna.mol.validator>` like charge correction based on the :doc:`OpenEye charge model <../api/luna.mol.charge_model>`, or valence correction. Since this solution was taken on basis of `Open Babel`_ 2.3.2, some errors from other versions may still persist.

Having that in mind, I prefer to use two `Open Babel`_ versions in my projects: **2.3.2** for parsing/amending/converting molecules and **3.1.1** as the standard Python library to wrap molecules as Python objects. However, by default, LUNA is configurated to use only the version 3.1.1 of `Open Babel`_. So, if you work with PDB files and are also interested to use this strategy, you can change the OPENBABEL default executable at **luna/util/default_values.py**. There, modify the OPENBABEL variable to the executable of your preference. For example, in my local environment, you will find OPENBABEL set to "/usr/bin/obabel" that is the version 2.3.3.

.. include:: ../substitutions.rst
