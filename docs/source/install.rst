============
Installation
============

.. contents::
   :local:


Dependencies
------------

LUNA is compatible with Python 3.x. It additionally has the following dependencies:

- Biopython_ == 1.72
- colorlog_
- Matplotlib_
- mmh3_ >= 2.5.1
- NetworkX_
- NumPy_
- `Open Babel`_
- pandas_
- Pymol_
- RDKit_
- SciPy_
- Seaborn_
- xopen_


Installation
------------

The following installation approaches are listed in order of recommendation.
These approaches require a prior installation of the dependencies listed above.

Option 1: Pip
~~~~~~~~~~~~~

To install with pip, run:

.. code:: bash

    pip install -U luna

We recommend using a virtual environment for this.

Option 2: Build from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download LUNA repository to your machine.

   -  Clone it to your machine with

      .. code:: bash

          $ git clone https://github.com/keiserlab/LUNA.git

   -  OR download an archive by navigating to the repository_ and clicking
      "Download ZIP". Extract the archive.

2. Install with

   .. code:: bash

       $ cd LUNA
       $ python setup.py build
       $ python setup.py install


We recommend using a virtual environment for this.

.. include:: substitutions.rst
