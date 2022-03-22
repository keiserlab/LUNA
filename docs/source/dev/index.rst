Developer notes
===============

.. _Ligands in PDB files:
Ligands in PDB files
--------------------

LUNA works with both `Open Babel`_ and RDKit_. You can switch between both libraries when setting a new project. So, the choice between `Open Babel`_ and RDKit_ depends on the type of molecular file you are working with. If you have ligands in MOL or MOL2 files, for instance, you can use RDKit_ for your convenience as you will be able to apply RDKit_ functions directly on the molecular object stored at :doc:`LUNA entries <../api/luna.mol.entry>`. However, if you frequently work with **ligands** in **PDB files**, you may be especially interested in reading this note.

Although each library has its own difficulties at parsing PDB files, we identified that `Open Babel`_ works better at parsing ligands from PDB files, especially what concerns ligands containing aromatic rings. That happens because PDB files don't contain information about atomic charges, valences, and bond types. Consequently, it may occur that some ligands are incorrectly perceived by `Open Babel`_ and RDKit_.

For this reason, when working with ligands in PDB files, no matter you choose `Open Babel`_ or RDKit_, we use `Open Babel`_ under the hood for parsing, amending, and converting ligands so that both final `Open Babel`_ and RDKit_ objects reflect at most the original and expected molecular structure.

However, in the past, I identified different errors in different `Open Babel`_ versions, which can be seen in this `issue thread <https://github.com/openbabel/openbabel/issues/1925>`_. During that time, I identified that `Open Babel`_ 2.3.2 works better for parsing PDB files. For the most common errors identified in this version, I implemented :doc:`amending solutions <../api/luna.mol.validator>` like charge correction based on the :doc:`OpenEye charge model <../api/luna.mol.charge_model>`, or valence correction. Since this solution was taken on basis of `Open Babel`_ 2.3.2, some errors from other versions may still persist.

Having that in mind, I usually prefer to use two `Open Babel`_ versions in LUNA projects: **2.3.2** for parsing/amending/converting molecules and **3.1.1** as the standard Python library to wrap molecules as Python objects. However, by default, LUNA is configurated to use only the version 3.1.1 of `Open Babel`_. So, if you work with PDB files and are also interested to use this strategy, you can change the OPENBABEL default executable at **luna/util/default_values.py**. There, modify the variable ``OPENBABEL``  to the executable of your preference. For example, in my local environment, you will find ``OPENBABEL`` set to "/usr/bin/obabel" that is the version 2.3.3.


Authoring code
--------------

We welcome contributions to LUNA!!! These notes are designed to help developers contribute code.

Code Formatting
~~~~~~~~~~~~~~~

LUNA's code should be *readable*. To ensure this, we rigorously follow the PEP8_ style conventions
and PEP257_ docstring conventions, which maximize readability of the code and ease of future development.
You may check your code for conformation to these conventions with the pycodestyle_ and pydocstyle_ utilities,
respectively. Where the code is necessarily complicated, inline comments should reorient the reader.

Errors
~~~~~~

LUNA-specific errors should inherit `luna.util.exceptions.LUNAError` base class.
However, several built-in exceptions are also included in `luna.util.exceptions`.


Contributing Code
~~~~~~~~~~~~~~~~~

Before contributing code to LUNA, it is advisable for major modifications to
submit an issue to the `issue tracker`_ to enable other developers to contribute
to the design of the code and to reduce the amount of work necessary to conform
the code to LUNA's standards. After writing the code, create a `pull request`_.
This is best even if you have push access to the LUNA repo, as it enables the test
suite to be run on the new code prior to merging it with the remaining code base.


Versioning
~~~~~~~~~~

LUNA versioning system is inspired by the `Semantic Versioning Specification (SemVer)`_, in which
the version number takes on the following pattern MAJOR.MINOR.PATCH.

Given a version number, increment:

   * **MAJOR** version when you make big scientific updates, which will make results produced by LUNA backward-incompatible;
   * **MINOR** version when you make backward-incompatible updates that do not involve big scientific updates \
     (e.g., modify class and function parameters);
   * **PATCH** version when you make backward-compatible bug fixes or add new functionalities.

Writing Tests
~~~~~~~~~~~~~

The standard in LUNA is to commit a test for new functionality simultaneously
with the new functionality or within the same pull request. While this slows
development, it prevents building a large backlog of untested methods and
classes.

These should ideally be unit tests, though for some complicated
functionalities, integration tests are also necessary.
For these complicated functions, specific units may still be tested using
:py:mod:`unittest.mock`.


Continuous Integration
~~~~~~~~~~~~~~~~~~~~~~

LUNA uses `Travis CI`_ for continuous integration. This ensures that each commit
and pull request passes all tests on a variety of a systems and for all
supported versions of Python. Additionally, Travis CI updates code coverage on
Coveralls_ and tests all usage examples in the documentation using `doctest`.




Documentation
-------------

In general, it is best to document the rationale and basic usage of a module,
class, or method in its docstring instead of in a separate documentation file.
See, for example, the docstring for `~luna.interaction.calc.InteractionCalculator`.
We use a variety of tools to ensure that our documentation is always
up-to-date. The official documentation is hosted on ReadtheDocs_ and is
automatically generated when new code is committed to the repository.

Documenting Code
~~~~~~~~~~~~~~~~

LUNA uses NumPy's `docstring conventions`_ for all docstrings. These are
parsed by Sphinx_ using Napoleon_.

The purpose of a docstring is to explain the purpose of a class/method, any
relevant implementation details, its parameters, its attributes, its outputs,
and its usage. The goal is clarity. For self-evident methods with descriptive
variables, a simple one-line summary is all that is needed. For complicated use
cases, often involving other methods/classes, it is better to document the
usage elsewhere in the documentation.



.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _PEP257: https://www.python.org/dev/peps/pep-0257/
.. _pycodestyle: http://pycodestyle.pycqa.org/en/latest/
.. _pydocstyle: http://pydocstyle.pycqa.org/en/latest/
.. _docstring conventions: https://numpydoc.readthedocs.io/en/latest/format.html
.. _Napoleon: http://www.sphinx-doc.org/en/stable/ext/napoleon.html
.. _Sphinx: http://www.sphinx-doc.org/en/stable/index.html
.. _pull request: https://help.github.com/articles/creating-a-pull-request/
.. _Semantic Versioning Specification (SemVer): https://semver.org/spec/v2.0.0.html
.. _Travis CI: https://travis-ci.org/keiserlab/luna
.. _Coveralls: https://coveralls.io/github/keiserlab/luna


.. include:: ../substitutions.rst

