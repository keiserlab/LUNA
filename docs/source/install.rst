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

Option 1: Conda + Pip

1. Download LUNA's `environment.yml <https://github.com/keiserlab/LUNA/blob/master/luna-env.yml>`_ file.

2. Create the Conda environment using the downloaded file:

    .. code:: bash

        $ conda env create -f <LUNA-ENV-FILE>

3. After creating the Conda environment, activate it:

    .. code:: bash

        $ conda activate luna-env

4. Finally, install LUNA from Pip:

    .. code:: bash

        $ pip install luna


Option 2: Pip
~~~~~~~~~~~~~

To install with pip, run:

.. code:: bash

    pip install -U luna

We recommend using a virtual environment for this.

.. note::

    This approach requires a prior installation of the dependencies listed above.

Option 3: Build from source
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

.. note::

    This approach requires a prior installation of the dependencies listed above.


.. include:: substitutions.rst

