Getting started
##################

.. highlight:: bash

Installing ``pyhma``
=====================

``pyhma`` can be directlly installed from `Python package index <https://pypi.python.org/pypi/pyhma>`_ using ``pip`` command (Python 3.x):

.. code-block:: bash

   pip install pyhma



Example (``pyhma`` package with VASP)
=====================================

* Using Python 3.x interpreter:

.. code-block:: python

   >>> import pyhma
   >>> r = pyhma.ReadVASP('OUTCAR1','OUTCAR2','OUTCAR3')
   >>> r.read()
   >>> import pyhma
   >>> sim = pyhma.Simulation()
   >>> sim.run() 
   >>> data = sim.get_statistics() 
   >>> print(' Anharmonic energy [eV/atom]:')
   >>> print(' conv.:', data['uahc'])
   >>> print(' hma  :', data['uahm'])
   >>> print(' Anharmonic pressure [GPa]:')
   >>> print(' conv.:', data['pahc'], '[GPa]')
   >>> print(' hma  :', data['pahm'], '[GPa]')

* Using pyhma script:

In the folder containing OUTCAR file/files, create a new file named pyhma.in which specifies the simulation temperature (K), excess harmonic pressure from lattice pressure (GPa), number of equilibrium steps, and blocksize respectively as shown below.

.. code-block:: bash

   TEMP 250
   PQH 5
   NEQ 1000
   BLOCK 1000
 
To calculate the HMA anharmonic energy and pressure:


.. code-block:: bash

   pyhma --readvasp OUTCAR1 OUTCAR2 OUTCAR3
   pyhma --compute 
