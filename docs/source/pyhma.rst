pyhma
######

``pyhma`` is an implementation of the `mapped-averaging <https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00018>`_ method for crystalline systems, where the harmonic behavior is used to define the mapping "velocity"; hence the name harmonically-mapped averaging (hma). The ensemble averages are obtained in ``pyhma`` via post-processing of outputs from molecular simulation codes (currently, VASP) in two steps:

**1. Reading raw data**
   Raw data (lattice vectors, configurations, forces, energies, and pressures) are extracted from molecular simulation code output and written to the following .dat files:

   * **lattice.dat [Å]:** contains lattice vectors and coordinates of the initial (minimized) system.
   * **posfor.dat [Å and eV/Å]:** contains the atomic positions and forces at each simulation step.
   * **u.dat  [eV]:** contains the total potential energy of the supercell at each step.
   * **p.dat [GPa]:** contains the virial pressure (:math:`-dU/dV`, with uniform scaling of coordinates) at each step.

For the case of VASP, this step is done using the :class:`pyhma.ReadVASP` class to read the OUTCAR files. Below is an example of reading two OUTCAR files of FCC Al system of 32 atoms, under high pressure (100 GPa) and temperature (1000 K). This example can be found in ``pyhma/examples/Al/`` directory.

.. code-block:: python

  >>> import pyhma
  >>> r = pyhma.ReadVASP('OUTCAR1','OUTCAR2')
  >>> r.read()

Or, it can be invoked using the ``pyhma`` script:

.. code-block:: bash

   pyhma --readvasp OUTCAR1 OUTCAR2


The output from this step should look like this: 

.. code-block:: python

 32 atoms (total)
 Lattice vectors (Angs)
 ----------------------
 6.861855122000  0.000000000000  0.000000000000
 0.000000000000  6.861855122000  0.000000000000
 0.000000000000  0.000000000000  6.861855122000

 Atom    initial xyz coordinates (Angs)                xyz forces (eV/Ang)
 --------------------------------------------------------------------------------
  1      0.000000   0.000000   0.000000        -0.0000000  -0.0000040   0.0000000
  2      3.430930   0.000000   0.000000         0.0000000  -0.0000040   0.0000000
  3      0.000000   3.430930   0.000000        -0.0000000   0.0000040   0.0000000
  4      3.430930   3.430930   0.000000         0.0000000   0.0000040   0.0000000
  5      0.000000   0.000000   3.430930        -0.0000000  -0.0000040  -0.0000000
  6      3.430930   0.000000   3.430930         0.0000000  -0.0000040  -0.0000000
  7      0.000000   3.430930   3.430930        -0.0000000   0.0000040  -0.0000000
  8      3.430930   3.430930   3.430930         0.0000000   0.0000040  -0.0000000
  9      1.715460   1.715460   0.000000        -0.0000070  -0.0000090   0.0000000
  10      5.146390   1.715460   0.000000         0.0000070  -0.0000090   0.0000000
  11      1.715460   5.146390   0.000000        -0.0000070   0.0000090   0.0000000
  12      5.146390   5.146390   0.000000         0.0000070   0.0000090   0.0000000
  13      1.715460   1.715460   3.430930        -0.0000070  -0.0000090  -0.0000000
  14      5.146390   1.715460   3.430930         0.0000070  -0.0000090  -0.0000000
  15      1.715460   5.146390   3.430930        -0.0000070   0.0000090  -0.0000000
  16      5.146390   5.146390   3.430930         0.0000070   0.0000090  -0.0000000
  17      0.000000   1.715460   1.715460        -0.0000000  -0.0000090   0.0000060
  18      3.430930   1.715460   1.715460         0.0000000  -0.0000090   0.0000060
  19      0.000000   5.146390   1.715460        -0.0000000   0.0000090   0.0000060
  20      3.430930   5.146390   1.715460         0.0000000   0.0000090   0.0000060
  21      0.000000   1.715460   5.146390         0.0000000  -0.0000090  -0.0000060
  22      3.430930   1.715460   5.146390         0.0000000  -0.0000090  -0.0000060
  23      0.000000   5.146390   5.146390        -0.0000000   0.0000090  -0.0000060
  24      3.430930   5.146390   5.146390         0.0000000   0.0000090  -0.0000060
  25      1.715460   0.000000   1.715460        -0.0000070  -0.0000040   0.0000060
  26      5.146390   0.000000   1.715460         0.0000070  -0.0000040   0.0000060
  27      1.715460   3.430930   1.715460        -0.0000070   0.0000040   0.0000060
  28      5.146390   3.430930   1.715460         0.0000070   0.0000040   0.0000060
  29      1.715460   0.000000   5.146390        -0.0000070  -0.0000040  -0.0000060
  30      5.146390   0.000000   5.146390         0.0000070  -0.0000040  -0.0000060
  31      1.715460   3.430930   5.146390        -0.0000070   0.0000040  -0.0000060
  32      5.146390   3.430930   5.146390         0.0000070   0.0000040  -0.0000060


 Lattice energy (eV/atom):  -2.93845111125
 Lattice pressure   (GPa):  99.464317



.. warning::

   The configuration of the first molecular simulation step must be the one with minimized energy (i.e., zero force on each).



**2. Computing averages**

In this step, the .dat files generated in the first step are used to compute (conventional and hma) anharmonic ensemble averages (currently, energy and pressure) using the :class:`pyhma.Simulation` class. Since the .dat files are universal, this step is independent on the molecular simulation code. This step uses a user-generated ``pyhma.in`` input file, which contains the following information (again, using the above Al system):::

   # pyhma.in input file for step #2
   TEMP    1000 # temperature (K)
   PQH     4.0  # quasiharmonic pressure (GPa/K)
   NEQ     100  # equilibration steps
   BLOCK   100  # blocksize (in steps)

Here is an example of how to use the :class:`pyhma.Simulation` class within the Python interpreter: 

.. code-block:: python

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

Or, from the command line using the pyhma script:

.. code-block:: bash

   pyhma --compute

The output from this step should be 

.. code-block:: python

   Run ....
   Done. Found 900 steps after 100 steps of equilibaration
   
   Temperature           (K): 1000.0
   Volume       (Ang^3/atom): 10.096588460306489
   Lattice energy  (eV/atom): -0.0918265972265625
   Harmonic energy (eV/atom): 0.1252205858238222
   Lattice pressure    (GPa): 99.464317
   Harmonic pressure   (GPa): 4.0 

   Block averaging statistics
   ==========================
   Using  9  blocks of size  100  MD steps 

   Anharmonic energy [eV/atom]:
    conv.: {'avg': -2.8459263453288584, 'err': 0.001774507316355033, 'corr': -0.07466066266723995}
    hma  : {'avg': -2.8477249829084768, 'err': 0.0002183231516704866, 'corr': -0.31459267645202793} 

   Anharmonic pressure [GPa]:
    conv.: {'avg': 1.1533150204254468, 'err': 0.04680482720426752, 'corr': -0.046065939075653525}
    hma  : {'avg': 1.1155015968330577, 'err': 0.01369902534482252, 'corr': 0.0741058952868103}





.. autoclass:: pyhma.ReadVASP
   :members:

.. autoclass:: pyhma.Simulation
   :members:


