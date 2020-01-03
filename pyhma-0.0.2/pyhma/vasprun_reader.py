import numpy as np
import sys, os
import lxml.etree 

def read(vasprun_files, force_tol=0.001, raw_files=False, verbose=False):
  box_row_vecs  = [] # box edge (row) vectors  
  basis          = [] 
  positions      = [] 
  forces         = [] 
  energies       = [] 
  pressures_vir  = [] 
  n_files        = len(vasprun_files)
  parser         = lxml.etree.XMLParser(recover=True) 

  for i, vasprun_file_i in enumerate(vasprun_files): 
    tree = lxml.etree.parse(vasprun_file_i, parser) 
    if verbose: 
      if i == 0: 
        print('Reading' , *vasprun_files, 'files') 
        print('==============================================') 
        print('  first configuration data from', vasprun_files[0]) 
        print('  ----------------------------------------------') 
      else: 
        print('  Reading' , vasprun_file_i ,'data (', i+1,'out of' , n_files, ') ... ')

    # Extract simulation parameters
    if i == n_files-1: 
      timestep = float(tree.find("./incar/i[@name='POTIM']").text) 
      temperature  = float(tree.find("./incar/i[@name='TEBEG']").text) 

    # Extract first step information
    if i == 0: 
      num_atoms = int(tree.find("./atominfo/atoms").text) # number of atoms
      volume  = float(tree.find("./structure/crystal/i[@name='volume']").text) # volume
      # box edge (row) vectors
      for v in tree.find("./structure/crystal/varray[@name='basis']"): 
        box_row_vecs.append([float(x) for x in v.text.split()]) 
      # initial positions (must be equilibrium/lattice positions)
      for v in tree.find("./structure[@name='initialpos']/varray[@name='positions']"): 
        basis.append([float(x) for x in v.text.split()]) 

      # extract initial forces (must be smaller than the user-defined force tolerance, force_tol)
      forces_0 = [] 
      for v in tree.find("./calculation/varray[@name='forces']"): 
        forces_0.append([float(x) for x in v.text.split()]) 
      # check for equilubrium
      if _is_large_force(forces_0, force_tol) == True:
        print(' MUST START FROM MINIMIZED CONFIGURATION (ZERO FORCES).')
        print(' EXITING pyHMA!\n')
        return

    # Extract positions
    for pos_elem in tree.findall("./calculation/structure/varray[@name='positions']"): 
      r = [] 
      for v in pos_elem: 
        r.append([float(x) for x in v.text.split()]) 
      positions.append(r) 
   
    # Extract forces
    for for_elem in tree.findall("./calculation/varray[@name='forces']"): 
      f = [] 
      for v in for_elem: 
        f.append([float(x) for x in v.text.split()]) 
      forces.append(f) 
    
    # Extract energies (E0)
    for e in tree.findall("./calculation/energy/i[@name='e_wo_entrp']"): 
      energies.append(float(e.text)) 
 
     
    # Extract virial pressures (i.e., total pressure - ideal gas)
    for pvir_elem in tree.findall("./calculation/varray[@name='stress']"): 
      pvir = 0 
      for n, v in enumerate(pvir_elem): 
        pvir += float(v.text.split()[n]) 
      pvir /= 3.0
      pvir /= 10.0 # convert kbar to GPa 
      pressures_vir.append(pvir) 
    
    # resize the length of all data lists
    min_len =  len(tree.findall("./calculation/energy/i[@name='total']"))
    positions     = positions[0:min_len]
    forces        = forces[0:min_len]
    energies      = energies[0:min_len]
    pressures_vir = pressures_vir[0:min_len]

  if raw_files: 
    _make_raw_files(num_atoms, box_row_vecs, basis, positions, forces, energies, pressures_vir)
 
  return {'box_row_vecs': box_row_vecs, 'num_atoms': num_atoms, 'volume': volume, 'basis': basis, 'positions': positions,\
          'forces': forces, 'energies': energies, 'pressures_vir': pressures_vir, 'timestep': timestep, 'temperature': temperature}


def _is_large_force(f, tol):
  is_large_force = False
  for i,v in enumerate(f):
    v_magn = np.linalg.norm(v)
    if v_magn > tol:
      print(' WARNING! Magnitude of force on atom ',(i+1),' = ', v_magn, ' (> tol = ', tol,')')
      is_large_force = True
  return is_large_force


def  _make_raw_files(num_atoms, box_row_vecs, basis, positions, forces, energies, pressures_vir):
    file_poscar_eq     = open('poscar_eq.dat', 'w')
    file_posfor        = open('posfor.dat', 'w')
    file_energies      = open('energies.dat', 'w')
    file_pressures_vir = open('pressures_vir.dat', 'w')
    # generate poscar_eq file 
    print('Lattice vectors\n1.0 scaling factor', file = file_poscar_eq)
    np.savetxt(file_poscar_eq, box_row_vecs, fmt='%12.8f')
    print(num_atoms,' atoms (total)\nDirect',file=file_poscar_eq)
    np.savetxt(file_poscar_eq, basis, fmt='%12.8f')
    # generate posfor file
    for t in range(len(positions)):
      print(t, file=file_posfor)
      np.savetxt(file_posfor, np.concatenate((positions[t] ,forces[t]),axis=1), fmt='%12.8f')
    # generate energies file
    np.savetxt(file_energies, energies, fmt='%12.8f')
    # generate virial pressure file
    np.savetxt(file_pressures_vir, pressures_vir, fmt='%12.8f')
 
    file_poscar_eq.close()
    file_posfor.close()
    file_energies.close()
    file_pressures_vir.close()
 
