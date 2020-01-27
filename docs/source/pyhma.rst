.. _pyhma_package:

################
pyHMA package
################


``pyHMA`` is a Python implementation of the `Mapped-Averaging <https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00018>`_ method for *crystalline systems*, where the harmonic character is used to define the mapping "velocity"; hence the name harmonically-mapped averaging (HMA). The ensemble averages are obtained in ``pyHMA`` via postprocessing of ``vasprun.xml`` output file(s) of VASP AIMD simulation. See :ref:`application_crystal` for details on HMA formulas.

``pyHMA`` is free software: you can modify and/or redistribute it under the terms of the Mozilla Public License.


.. toctree::
   :maxdepth: 1
   
   pyhma_installation
   pyhma_structure
   pyhma_usage
