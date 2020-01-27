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
A module for extracting data from vasprun.xml AIMD output file(s) so that it can 
be analyzed with the processor.py module.

Please cite this paper of you use pyHMA package in your research:

Sabry G. Moustafa, Apoorva Purohit, Andrew J. Schultz, and David A. Kofke,
"pyHMA: A Python Package for Precise Measurement of Crystalline Anharmonic 
Properties using Harmonically Mapped Averaging", Journal of XXXXXXX
"""


import numpy as np
import lxml.etree 

def read(vasprun_files, force_tol=0.001, raw_files=False, verbose=False):
  """
  Read AIMD data from vasprun.xml output file(s).

  Parameters
  -----------
   a : ww
   b : www

  Returns
  --------
  box_row_vecs : 2D list, float, shape=(3,3)
    hey man! I am here!!
  num_atoms : int
    hey!
  volume_atom : float
    hey!!
  basis : 2D list, float, shape=(num_atoms, 3)
    hey!!!
  """
 

  box_row_vecs = [] # box edge (row) vectors  
  basis        = [] 
  position     = [] 
  force        = [] 
  energy       = [] 
  pressure_vir = [] 
  n_files      = len(vasprun_files)
  parser       = lxml.etree.XMLParser(recover=True) 
  list_len     = 0  

  for i, vasprun_file_i in enumerate(vasprun_files): 
    tree = lxml.etree.parse(vasprun_file_i, parser) 
    if verbose: 
      if i == 0: 
        print('\nReading' , *vasprun_files) 
        print('==============================================') 
        print(' first configuration data from', vasprun_files[0]) 
        print(' ----------------------------------------------') 

    # Extract simulation parameters
    if i == n_files-1: 
      timestep = float(tree.find("./incar/i[@name='POTIM']").text) 
      temperature  = float(tree.find("./incar/i[@name='TEBEG']").text) 

    # Extract first step information
    if i == 0: 
      num_atoms   = int(tree.find("./atominfo/atoms").text) # number of atoms
      volume_atom = float(tree.find("./structure/crystal/i[@name='volume']").text)/num_atoms # specific volume
      # box edge (row) vectors
      for v in tree.find("./structure/crystal/varray[@name='basis']"): 
        box_row_vecs.append([float(x) for x in v.text.split()]) 
      # initial positions (must be equilibrium/lattice positions)
      for v in tree.find("./structure[@name='initialpos']/varray[@name='positions']"): 
        basis.append([float(x) for x in v.text.split()]) 

      # extract initial forces (must be smaller than the user-defined force tolerance, force_tol)
      force_0 = [] 
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

      # check for equilubrium
      if _is_large_force(force_0, force_tol) == True:
        print(' MUST START FROM MINIMIZED CONFIGURATION (ZERO FORCES).')
        print(' EXITING pyHMA!\n')
        return

    if verbose:
      print('  Reading' , vasprun_file_i ,' (', i+1,'out of' , n_files, ')')

    # Extract positions
    for pos_elem in tree.findall("./calculation/structure/varray[@name='positions']"): 
      r = [] 
      for v in pos_elem: 
        r.append([float(x) for x in v.text.split()]) 
      position.append(r) 
   
    # Extract forces
    for for_elem in tree.findall("./calculation/varray[@name='forces']"): 
      f = [] 
      for v in for_elem: 
        f.append([float(x) for x in v.text.split()]) 
      force.append(f) 
    
    # Extract energies (E0)
#    for e in tree.findall("./calculation/energy/i[@name='e_wo_entrp']"): correct for the wrong VASP!
#    for e in tree.findall("./calculation/energy/i[@name='e_0_energy']"): wrong for the correct VASP!!
#      energy.append(float(e.text)/num_atoms) 
    for e in tree.findall("./calculation"):
      ee = e.findall("./scstep/energy/i[@name='e_0_energy']")[-1]
      energy.append(float(ee.text)/num_atoms) 
 
     
    # Extract virial pressures (i.e., total - ideal gas)
    for pvir_elem in tree.findall("./calculation/varray[@name='stress']"): 
      pvir = 0 
      for n, v in enumerate(pvir_elem): 
        pvir += float(v.text.split()[n]) 
      pvir /= 3.0
      pvir /= 10.0 # convert kbar to GPa 
      pressure_vir.append(pvir) 
    
    # resize the length of all data lists
    list_len +=  len(tree.findall("./calculation/energy/i[@name='total']"))
    position     = position[0:list_len]
    force        = force[0:list_len]
    energy      = energy[0:list_len]
    pressure_vir = pressure_vir[0:list_len]
  if raw_files: 
    _make_raw_files(num_atoms, box_row_vecs, basis, position, force, energy, pressure_vir)



  return {'box_row_vecs': box_row_vecs, 'num_atoms': num_atoms, 'volume_atom': volume_atom, 'basis': basis, 'position': position,\
          'force': force, 'energy': energy, 'pressure_vir': pressure_vir, 'timestep': timestep, 'temperature': temperature}


def _is_large_force(f, tol):

  """checks if the magnitude of force on any atom is larger than the tol tolerance.

  """

  is_large_force = False
  for i,v in enumerate(f):
    v_magn = np.linalg.norm(v)
    if v_magn > tol:
      print(' WARNING! Magnitude of force on atom ',(i+1),' is ', v_magn, ' (> tol = ', tol,').')
      is_large_force = True
  return is_large_force


def  _make_raw_files(num_atoms, box_row_vecs, basis, position, force, energy, pressure_vir):
  """generate raw data files from vasprun.xml. 

  """
  
  file_poscar_eq    = open('poscar_eq.dat', 'w')
  file_posfor       = open('posfor.dat', 'w')
  file_energy       = open('energy.dat', 'w')
  file_pressure_vir = open('pressure_vir.dat', 'w')
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
  for i in range(len(pressure_vir)):
    print('%12.8f' % (pressure_vir[i]) , file=file_pressure_vir)
 
  # close all files 
  file_poscar_eq.close()
  file_posfor.close()
  file_energy.close()
  file_pressure_vir.close()
 
