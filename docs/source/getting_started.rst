Getting Started
##################

.. highlight:: bash

Installing ``pyhma``
=====================

``pyhma`` can be directlly installed from `Python package index <https://pypi.python.org/pypi/pyhma>`_ using ``pip`` command:

.. code-block:: bash

   pip install pyhma



Example
========

.. code-block:: python

   >>> import pyhma
   >>> r = pyhma.ReadVASP('OUTCAR1','OUTCAR2','OUTCAR3')
   >>> r.read()
   >>> import pyhma
   >>> sim = pyhma.Simulation()
   >>> sim.run()  # computes conventional and hma estimates for each configuration
   >>> data = sim.get_statistics() # computes average, stochastic uncertainty (1:math:`\sigma`), and correlation
   >>> print(' Anharmonic energy [eV/atom]:')
   >>> print(' conv.:', data['uahc'])
   >>> print(' hma  :', data['uahm'])
   >>> print(' Anharmonic pressure [GPa]:')
   >>> print(' conv.:', data['pahc'], '[GPa]')
   >>> print(' hma  :', data['pahm'], '[GPa]')

