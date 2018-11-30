pyhma package
#############

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
  >>> r = pyhma.ReadVASP('OUTCAR1','OUTCAR2','OUTCAR3')
  >>> r.read()

Or, it can be invoked using the ``pyhma`` script:

.. code-block:: bash

   pyhma --readvasp OUTCAR1 OUTCAR2 OUTCAR3


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
  1      0.000000   0.000000   0.000000         0.0000780  -0.0001250   0.0000010 
  2      3.430930   0.000000   0.000000        -0.0000790  -0.0001220   0.0000040 
  3      0.000000   3.430930   0.000000         0.0000760   0.0001240   0.0000030 
  4      3.430930   3.430930   0.000000        -0.0000790   0.0001250   0.0000020 
  5      0.000000   0.000000   3.430930         0.0000800  -0.0001250  -0.0000040 
  6      3.430930   0.000000   3.430930        -0.0000770  -0.0001260  -0.0000020 
  7      0.000000   3.430930   3.430930         0.0000770   0.0001230  -0.0000040 
  8      3.430930   3.430930   3.430930        -0.0000790   0.0001260   0.0000010 
  9      1.715460   1.715460   0.000000        -0.0001040  -0.0001170   0.0000020 
  10      5.146390   1.715460   0.000000         0.0001060  -0.0001140   0.0000050 
  11      1.715460   5.146390   0.000000        -0.0001060   0.0001140   0.0000020 
  12      5.146390   5.146390   0.000000         0.0001080   0.0001160   0.0000010 
  13      1.715460   1.715460   3.430930        -0.0001060  -0.0001170   0.0000010 
  14      5.146390   1.715460   3.430930         0.0001060  -0.0001150  -0.0000070 
  15      1.715460   5.146390   3.430930        -0.0001040   0.0001170  -0.0000050 
  16      5.146390   5.146390   3.430930         0.0001050   0.0001160   0.0000010 
  17      0.000000   1.715460   1.715460         0.0000750  -0.0001150   0.0000070 
  18      3.430930   1.715460   1.715460        -0.0000780  -0.0001170   0.0000070 
  19      0.000000   5.146390   1.715460         0.0000790   0.0001170   0.0000070 
  20      3.430930   5.146390   1.715460        -0.0000810   0.0001190   0.0000070 
  21      0.000000   1.715460   5.146390         0.0000830  -0.0001160  -0.0000080 
  22      3.430930   1.715460   5.146390        -0.0000800  -0.0001150  -0.0000040 
  23      0.000000   5.146390   5.146390         0.0000770   0.0001160  -0.0000060 
  24      3.430930   5.146390   5.146390        -0.0000750   0.0001150  -0.0000090 
  25      1.715460   0.000000   1.715460        -0.0001060  -0.0001280   0.0000080 
  26      5.146390   0.000000   1.715460         0.0001040  -0.0001230   0.0000060 
  27      1.715460   3.430930   1.715460        -0.0001060   0.0001250   0.0000060 
  28      5.146390   3.430930   1.715460         0.0001080   0.0001220   0.0000090 
  29      1.715460   0.000000   5.146390        -0.0001060  -0.0001200  -0.0000060 
  30      5.146390   0.000000   5.146390         0.0001050  -0.0001290  -0.0000080 
  31      1.715460   3.430930   5.146390        -0.0001040   0.0001230  -0.0000070 
  32      5.146390   3.430930   5.146390         0.0001030   0.0001250  -0.0000080 


 Lattice energy (eV/atom):  -2.9800631290625
 Lattice pressure   (GPa):  98.638139



.. warning::

   The configuration of the first molecular simulation step must be the one with minimized energy (i.e., zero force on each).



**2. Computing averages**

In this step, the .dat files generated in the first step are used to compute (conventional and hma) anharmonic ensemble averages (currently, energy and pressure) using the :class:`pyhma.Simulation` class. Since the .dat files are universal, this step is independent on the molecular simulation code. This step uses a user-generated ``pyhma.in`` input file, which contains the following information (again, using the above Al system):::

   # pyhma.in input file for step #2 (Simulation)
   TEMP    1000 # temperature (K)
   PQH     5.5  # quasiharmonic pressure (GPa)
   NEQ     500 # equilibration steps
   BLOCK   30  # blocksize (in steps)

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
   Done. Found 1500 steps after 500 steps of equilibaration

   Temperature           (K): 1000.0
   Volume       (Ang^3/atom): 10.096588460306489
   Lattice energy  (eV/atom): -2.9800631290625
   Harmonic energy (eV/atom): 0.1252205858238222
   Lattice pressure    (GPa): 98.638139
   Harmonic pressure   (GPa): 5.5

   Block averaging statistics
   ==========================
   Using  50  blocks of size  30  MD steps

   Anharmonic energy [eV/atom]:
    conv.: {'avg': -0.0014000991084052998, 'err': 0.0019232457429067629, 'corr': 0.5626751626440605}
    hma  : {'avg': -0.0014448343942825206, 'err': 0.00018553677013852636, 'corr': 0.3045234376482494}

   Anharmonic pressure [GPa]:
    conv.: {'avg': -0.31463552438915376, 'err': 0.06244769350788043, 'corr': 0.6339961621258104}
    hma  : {'avg': -0.31611188863687606, 'err': 0.013142430793349463, 'corr': -0.1737875001083584}



.. autoclass:: pyhma.ReadVASP
   :members:

.. autoclass:: pyhma.Simulation
   :members:


