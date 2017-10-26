************
Installation
************

Requirements
============

The code for the X-ray polarimetry beamline depends on the following packages:

- `numpy <http://www.numpy.org/>`_
- `astropy`_
- `scipy <https://www.scipy.org/scipylib/index.html>`_

They are best installed with a package manager such as conda. See the `astropy installation instructions for a detailed discussion <https://astropy.readthedocs.io/en/stable/install.html>`_.

  
Install the python code
=======================

setup.py
--------

To download the latest version of the X-ray polarimetry code for the first time:

.. code-block:: bash

   $ git clone 
   $ cd ccdcamera-process

On the computers in the lab, we already have a directory that hold the git repository.
Simply go into that directory and do:

.. code-block:: bash

   $ git pull

Now you install, run tests or build the documentation:

.. code-block:: bash

   $ python setup.py install
   $ python setup.py test
   $ python setup.py build_docs
