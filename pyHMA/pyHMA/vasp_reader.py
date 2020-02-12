########################################################################
# pyHMA: A Python Library for HMA 
# 
# Copyright (c) 2020 University at Buffalo
# 
# Authors: Sabry Moustafa, Andrew Schultz, and David Kofke 
# 
# pyHMA is free software: you can modify and/or redistribute it under 
# the terms of the Mozilla Public License.
#
# pyHMA is distributed in the hope that it will be useful, but without 
# any warranty. See the Mozilla Public License for more details.
#
########################################################################

"""
**Overview**

A module for extracting data from ``vasprun.xml`` output file(s) of VASP AIMD simulation so that it can be used by the :py:mod:`pyhma.processor` module to compute anharmonic properties, using Conv and HMA methods.


"""


import numpy as np
import lxml.etree 

def read(vasprun_files, force_tol=0.001, raw_files=False, verbose=False):
  """
  A function that uses LXML parser to extract raw data from ``vasprun.xml`` file(s).

  Parameters
  -----------
  vasprun_files : list
    List of vasprun.xml files of the same AIMD simulation.
  force_tol : float
    Force tolerance (in eV/Å) on initial configuration. *Default: 0.001*.
  raw_files : bool
    If True, the following raw data files will be generated: energy.dat, poscar_eq.dat, posfor.dat, and pressure.dat. *Default: False*.
  verbose : bool
    If True, pyHMA will print simulation details while reading data. *Default: False*


  Returns
  -------
  box_row_vecs : list
    Box edge (row) vectors in Å.
  num_atoms : int
    Total number of atoms.
  volume_atom : float
    Volume per atom in Å^3.
  basis : list
    List of atomic fractional positions of first configuration.
  position : list
    Instantaneous atomic fractional positions.
  force : list
    Instantaneous atomic forces in eV/Å.
  energy : list
    Instantaneous potential energy (E0) in eV/atom. 
  pressure : list
    Instantaneous pressure in GPa. 
  pressure_ig : float
    Ideal gas pressure in GPa.
  timestep : float
    MD timestep in fs.
  temperature : float
    NVT set temperature in K.  


  Notes
  -----
  * A list of ``vasprun.xml`` files of continued simulation can be passed; e.g., read(['vasprun1.xml', 'vasprun2.xml']).
  * The read() function handles incomplete XML file(s) generated from interrupted AIMD runs (by the user, or due to some time constraint). This was possible by using the recover capability of the LXML parser.


  Example
  --------

  .. code-block:: python

      >>> import pyhma
      >>> data = pyhma.read(['vasprun-1.xml', 'vasprun-2.xml'], raw_files=True, force_tol=0.002, verbose=False)

  .. warning::

    The initial configuration of the first ``vasprun.xml`` file **must be** the lattice one (i.e., forces < force_tol).

  """

  box_row_vecs = [] # box edge (row) vectors  
  basis        = [] # List of atomic fractional positions of first configuration 
  position     = [] # Instantaneous atomic fractional positions 
  force        = [] # Instantaneous atomic forces in eV/Å
  energy       = [] # Instantaneous potential energy (E0) in eV/atom
  pressure     = [] # Instantaneous pressure in GPa
  n_files      = len(vasprun_files) # number of vasprun.xml files
  parser       = lxml.etree.XMLParser(recover=True) # LXML parser with the capability to handle broken (incomplete) XML files
  list_len     = 0 # total number of complete scf steps found in vasprun.xml files

  for i, vasprun_file_i in enumerate(vasprun_files): 
    tree = lxml.etree.parse(vasprun_file_i, parser) # parsing the whole vasprun.xml file 
    if verbose: 
      if i == 0: 
        print('\nReading' , *vasprun_files) 
        print('==============================================') 
        print(' first configuration data from', vasprun_files[0]) 
        print(' ----------------------------------------------') 

    # extract first step information
    if i == 0: 
      num_atoms   = int(tree.find("./atominfo/atoms").text) # total number of atoms
      volume_atom = float(tree.find("./structure/crystal/i[@name='volume']").text)/num_atoms # average volume per atom
      # box edge (row) vectors
      for v in tree.find("./structure/crystal/varray[@name='basis']"): 
        box_row_vecs.append([float(x) for x in v.text.split()]) 
      # initial positions (must be equilibrium configuration that minimizes energy)
      for v in tree.find("./structure[@name='initialpos']/varray[@name='positions']"): 
        basis.append([float(x) for x in v.text.split()]) 

      # extract initial forces (must be smaller than the user-defined force tolerance, force_tol)
      force_0 = [] # forces of atoms in the first configuration 
      for v in tree.find("./calculation/varray[@name='forces']"): 
        force_0.append([float(x) for x in v.text.split()])
      # print initial configuration structure and forces 
      if verbose:
        print('',num_atoms, 'atoms (total)')
        print(' Box edge (row) vectors')
        for l in range(len(box_row_vecs)):     
          print('%12.8f %12.8f %12.8f ' % (*box_row_vecs[l],))
        print('\n atom       xyz (direct) coordinates (A)                   xyz forces (eV/A)')
        for l in range(len(force_0)):
          print(' %3d  %12.8f %12.8f %12.8f    %12.8f %12.8f %12.8f' % (l+1, *basis[l], *force_0[l]))
        print('')

      # check for equilibrium 
      if _is_large_force(force_0, force_tol) == True:
        print(' MUST START FROM MINIMIZED CONFIGURATION (ZERO FORCES).')
        print(' EXITING pyHMA!\n')
        return

    # extract timestep (fs) and temperature (K), and compute ideal-gas pressure (GPa).
    if i == n_files-1: 
      timestep = float(tree.find("./incar/i[@name='POTIM']").text) # timestep (fs)
      temperature  = float(tree.find("./incar/i[@name='TEBEG']").text) # temperature (K)
      kB = 0.0000861733063733830                     # Boltzmann's constant (eV/K)
      eV2J = 1.602176634e-19                         # eV to Joules conversion factor
      kBT_eV = kB*temperature                        # kT (eV)
      kBT_J  = kBT_eV*eV2J                           # kT (J)
      pressure_ig = kBT_J/(volume_atom*1e-30)*1e-9   # ideal gas pressure (GPa)

    # print files being read
    if verbose:
      print('  Reading' , vasprun_file_i ,' (', i+1,'out of' , n_files, ')')

    # extract positions
    for pos_elem in tree.findall("./calculation/structure/varray[@name='positions']"): 
      r = [] 
      for v in pos_elem: 
        r.append([float(x) for x in v.text.split()]) 
      position.append(r) 
   
    # extract forces
    for for_elem in tree.findall("./calculation/varray[@name='forces']"): 
      f = [] 
      for v in for_elem: 
        f.append([float(x) for x in v.text.split()]) 
      force.append(f) 
    
    # extract energies (E0)
    for e in tree.findall("./calculation"):
      ee = e.findall("./scstep/energy/i[@name='e_0_energy']")
      if len(ee) != 0:
        energy.append(float(ee[-1].text)/num_atoms) 
     
    # extract virial pressures (i.e., total - ideal gas)
    for pvir_elem in tree.findall("./calculation/varray[@name='stress']"): 
      pvir = 0 
      for n, v in enumerate(pvir_elem): 
        pvir += float(v.text.split()[n]) 
      pvir /= 3.0
      pvir /= 10.0 # convert kbar to GPa 
      pressure.append(pvir) 
    
    # determine the number of complete scf steps in each vasprun.xml file
    list_len +=  len(tree.findall("./calculation/energy/i[@name='total']"))
    position     = position[0:list_len]
    force        = force[0:list_len]
    energy      = energy[0:list_len]
    pressure  = pressure[0:list_len]

  # compute the total pressure (i.e., virial + ideal gas)
  for j in range(len(pressure)):
    pressure[j] += pressure_ig

  # if raw_files=True, generate the following raw data files: poscar_eq.dat, posfor.dat, energy.dat, and pressure.dat.
  if raw_files: 
    _make_raw_files(num_atoms, box_row_vecs, basis, position, force, energy, pressure)

  # return a dict of the extracted data
  return {'box_row_vecs': box_row_vecs, 'num_atoms': num_atoms, 'volume_atom': volume_atom, 'basis': basis, 'position': position,\
          'force': force, 'energy': energy, 'pressure': pressure, 'pressure_ig': pressure_ig, 'timestep': timestep, 'temperature': temperature, }


def _is_large_force(force, force_tol):

  """
  Check if the magnitude of force (force) on any atom is larger than some tolerance (force_tol).

  """

  is_large_force = False
  for i,f in enumerate(force):
    force_mag = np.linalg.norm(force)
    if force_mag > force_tol:
      print(' WARNING! Magnitude of force on atom ',(i+1),' is ', force_mag, ' (> force_tol = ', force_tol,').')
      is_large_force = True
  return is_large_force


def  _make_raw_files(num_atoms, box_row_vecs, basis, position, force, energy, pressure):
  """
  Generate the following raw data files: poscar_eq.dat, posfor.dat, energy.dat, and pressure.dat. 

  """
  
  file_poscar_eq = open('poscar_eq.dat', 'w') # initial POSCAR file (in fractional coordinates)
  file_posfor    = open('posfor.dat', 'w') # instantaneous 
  file_energy    = open('energy.dat', 'w')
  file_pressure  = open('pressure.dat', 'w')
  # generate poscar_eq file 
  print('Lattice vectors\n1.0 scaling factor', file = file_poscar_eq)
  for i in range(len(box_row_vecs)):
    print('%12.8f %12.8f %12.8f' % (*box_row_vecs[i],) , file=file_poscar_eq)
  print(num_atoms,' atoms (total)\nDirect', file=file_poscar_eq)
  for i in range(len(basis)):
    print('%12.8f %12.8f %12.8f' % (*basis[i],) , file=file_poscar_eq)
 

  # generate posfor file
  for i in range(len(position)):
    print(i, file=file_posfor)
    for j in range(len(position[i])):  
      print('%12.8f %12.8f %12.8f    %12.8f %12.8f %12.8f' % (*position[i][j], *force[i][j]), file=file_posfor)

  # generate energy file
  for i in range(len(energy)):
    print('%12.8f' % (energy[i]) , file=file_energy)
 
  # generate virial pressure file
  for i in range(len(pressure)):
    print('%12.8f' % (pressure[i]) , file=file_pressure)
 
  # close all files 
  file_poscar_eq.close()
  file_posfor.close()
  file_energy.close()
  file_pressure.close()
 
