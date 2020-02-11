.. _pyhma_package:

################
pyHMA
################

``pyHMA`` is a VASP post-processor (written in Python 3) for precise measurment of crystalline *anharmonic* properties using `Harmonically Mapped Averaging (HMA) <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.043303>`_ methd.
It is based on post-processing ``vasprun.xml`` output file(s) obtained from NVT Born-Oppenheimer *ab initio* molecular dynamics (AIMD) simulation. 
See :numref:`Table %s <pyhma_eqs>` as an example for HMA expressions for anharmonic energy and pressure, along with direct/conventional (Conv) counterpart.


``pyHMA`` is free software: you can modify and/or redistribute it under the terms of the Mozilla Public License.

Please cite this paper when using pyHMA package in your research:

  Sabry G. Moustafa, Apoorva Purohit, Andrew J. Schultz, and David A. Kofke, pyHMA: A VASP Post-processor for Precise Measurement of CrystallineAnharmonic Properties using Harmonically Mapped Averaging, Comput. Phys. Commun., 2020.


.. note::
   * The term *anharmonicity* is commonly used in literature to qualitatively descripe a system with no equilibrium configuration at :math:`0` K (i.e., imaginary frequencies); in other words, it referes to a "non-harmonic" potential energy surface. 
   * Here, however, we define *anharmonic contribution* of some property :math:`X` as the residual resulted caused by not using the harmonic approximation; i.e., :math:`X_{\rm ah} \equiv X - (X_{\rm lat} + X_{\rm qh})`. Therefore, this specific definition is meaningless if the system does not have equilibrium lattice configuration at :math:`T=0` K. For this reason, ``pyHMA`` checks forces on the first configuration to make sure the system has an equilibrium configuration (i.e., zero forces).



.. toctree::
   :maxdepth: 1
   
   pyhma_installation
   pyhma_usage
   pyhma_modules
