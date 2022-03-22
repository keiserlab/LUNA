.. image:: https://github.com/keiserlab/LUNA/blob/master/docs/source/_static/luna.svg


|Docs Status| |Build Status| |Coveralls Status| |PyPi Version| |Conda Version| |License|

LUNA [1]_ is an object-oriented Python 3 toolkit for drug design that makes it easy to analyze very large data sets of 3D molecular structures and complexes, and that allows identifying, filtering, and visualizing atomic interactions.

LUNA also implements three hashed interaction fingerprints (IFP): Extended Interaction FingerPrint (EIFP), Functional Interaction FingerPrint (FIFP), and Hybrid Interaction FingerPrint (HIFP) -- inspired by ECFP [2]_, FCFP [2]_, and E3FP [3]_. These IFPs encode molecular interactions at different levels of detail, provide several functionalities to trace individual bits back to their original atomic substructures in the context of the binding site, and are RDKit_-compatible.

Documentation is hosted by ReadTheDocs_, and development occurs on GitHub_.


Installation and Usage
----------------------

For installation and usage instructions, see the `documentation <http://luna-toolkit.readthedocs.io>`_.

See the LUNA `paper repository`_ for an application of LUNA and all code used for the LUNA paper [1]_.


License
-------

LUNA is available under the |license|.



References
----------

.. [1] |afassio2022|
.. [2] |rogers2010|
.. [3] |axen2017|

.. substitutions

.. |license| replace:: `MIT License`_
.. _MIT License: https://github.com/keiserlab/LUNA/blob/master/LICENSE


.. _RDKit: http://www.rdkit.org
.. _GitHub: https://github.com/keiserlab/LUNA
.. _paper repository: https://github.com/keiserlab/luna-paper
.. _ReadTheDocs: http://luna-toolkit.readthedocs.io
.. |afassio2022| replace:: TO BE DEFINED (2022).
.. |axen2017_doi| image:: https://img.shields.io/badge/doi-10.1021/acs.jmedchem.7b00696-blue.svg
    :target: http://dx.doi.org/10.1021/acs.jmedchem.7b00696
    :alt: Access the paper
.. |axen2017| replace:: Axen SD, Huang XP, Caceres EL, Gendelev L, Roth BL, Keiser MJ. A Simple Representation Of Three-Dimensional Molecular Structure. *J. Med. Chem.* **60** (17): 7393–7409 (2017). |axen2017_doi| |bioRxiv| |F1000 recommended|
.. |rogers2010_doi| image:: https://img.shields.io/badge/doi-10.1021/ci100050t-blue.svg
    :target: http://dx.doi.org/10.1021/ci100050t
    :alt: Access the paper
.. |rogers2010| replace:: Rogers D & Hahn M. Extended-connectivity fingerprints. *J. Chem. Inf. Model.* **50**: 742-54 (2010). |rogers2010_doi|
.. |Build Status| image:: https://travis-ci.org/keiserlab/luna.svg?branch=master
   :target: https://travis-ci.org/keiserlab/luna
   :alt: Build Status
.. |Docs Status| image:: http://readthedocs.org/projects/luna/badge/?version=latest
   :target: http://luna-toolkit.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. |Coveralls Status| image:: https://coveralls.io/repos/github/keiserlab/luna/badge.svg?branch=master
   :target: https://coveralls.io/github/keiserlab/luna?branch=master
   :alt: Code Coverage
.. |PyPi Version| image:: https://img.shields.io/pypi/v/luna.svg
   :target: https://pypi.python.org/pypi/luna
   :alt: Package on PyPi
.. |Conda Version| image:: https://img.shields.io/conda/v/keiserlab/luna.svg
   :target: https://anaconda.org/keiserlab/luna
   :alt: Package on Anaconda
.. |License| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://github.com/keiserlab/LUNA/blob/master/LICENSE
.. |F1000 recommended| image:: http://cdn.f1000.com.s3.amazonaws.com/images/badges/badgef1000.gif
   :target: http://f1000.com/prime/727824514?bd=1
   :alt: Access the recommendation on F1000Prime
   :width: 120px
   :scale: 75 %
.. |bioRxiv| image:: https://img.shields.io/badge/bioRxiv-136705-blue.svg
    :target: https://doi.org/10.1101/136705
    :alt: Access the preprint on bioRxiv
