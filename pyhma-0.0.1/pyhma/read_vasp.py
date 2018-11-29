##############################################################################
# pymbar: A Python Library for pyhma
#
# Copyright 2018 University at Buffalo
#
# See the MIT License for more details. 
##############################################################################

"""A module to read VASP OUTCAR data for pyhma

"""


import numpy as np
import sys

class ReadVASP:

  """

  Notes
  -------

  After ReadVASP class


  """

  
  def __init__(self, *filenames):  # tuple of filenames
    """ Initialize ReadVASP ...

    """
    self.filenames = filenames
    self.force_tol = 1e-3

  def read(self):

    """ Read data from OUTCAR file

    Parameters:
    ------------

    """
    # create: energy (u.dat), pressure (p.dat), positions-forces (posfor.dat), and lattice (lattice.dat)
    file_u         = open('u.dat', 'w')
    file_p         = open('p.dat', 'w')
    file_posfor    = open('posfor.dat', 'w')
    file_lattice = open('lattice.dat', 'w')
    with open('u.dat','w') as file_u, open('p.dat','w') as file_p , open('posfor.dat','w') as file_posfor, open('lattice.dat','w') as file_lattice:  
      # loop over OUTCARs
      is_first_step = True
      is_first_lat = False
      is_lat = True
      is_basis = False
      is_posfor = False
      i_atom = 0
      for i in range(len(self.filenames)):
        with open(self.filenames[i]) as  file_outcar:
          for line in file_outcar:
            # Create lattice.dat file
            if is_first_step:
              bits = line.split()
              if ' SYSTEM =' in line:
                print('Lattice vectors', file = file_lattice)
                print(' 1.0 scaling factor' , file = file_lattice)
                print('',self.n_atoms , 'atoms (total)')
                print(' Lattice vectors (Angs)')
                print(' ----------------------')
              if 'number of ions' in line:
                self.n_atoms = int(bits[len(bits)-1])

              # lattice
              if is_first_lat:
                if line.strip():
                  print('% .12f % .12f % .12f' % (float(bits[0]), float(bits[1]), float(bits[2])))
                  print('% .12f % .12f % .12f' % (float(bits[0]), float(bits[1]), float(bits[2])),file = file_lattice)
                  self.lat_vecs[i_lat_vecs] = [float(bits[0]), float(bits[1]), float(bits[2])]
                  i_lat_vecs+=1
                else:
                  print('', self.n_atoms,' atoms (total)',file=file_lattice)
                  print('Basis (Angst)', file=file_lattice)
                  is_lat = False
                  is_first_lat = False
              if ('direct lattice vectors' in line) & is_lat:
                is_first_lat = True
                self.lat_vecs = np.zeros((3,3))
                i_lat_vecs = 0

        
              # basis
              if is_basis:
                if line.strip():
                  print('% .12f % .12f % .12f' % (float(bits[0]), float(bits[1]), float(bits[2]))\
                      ,file = file_lattice)
                else:
                  is_basis = False
              if 'position of ions in cartesian coordinates' in line:
                is_basis = True

            # end of is_first_step case


            # extract pressure (without IG part)  (GPa)
            if '  in kB  ' in line:
              bits=line.split()
              p_vir = (float(bits[2])+float(bits[3])+float(bits[4]))/30.0
              if is_first_step:
                self.p_lat = p_vir
              print('%.10f' % p_vir , file=file_p) 

            # extract positions and forces
            if is_posfor:
              if ('----------' not in line) and ('total drift:' not in line):
                bits = line.split()
                if is_first_step:
                  if i_atom == 0:
                    print('\n','Atom    initial xyz coordinates (Angs)                xyz forces (eV/Ang)   ')
                    print(' --------------------------------------------------------------------------------')
                    initial_forces = np.zeros((self.n_atoms,3))
                  print(' % d     % .6f  % .6f  % .6f        % .7f  % .7f  % .7f ' % (i_atom+1, float(bits[0]),\
                      float(bits[1]),float(bits[2]),float(bits[3]),float(bits[4]),float(bits[5]))) 
                  initial_forces[i_atom] = [float(bits[3]),float(bits[4]),float(bits[5])]
                  i_atom += 1
                  if i_atom == self.n_atoms:
                    if self._get_max_forces(initial_forces) == True:
                      print('\n First configuration must be minimized (i.e., forces = 0)')
                      print(' EXIT pyhma!\n')
                      sys.exit()
                print('% .6f  % .6f  % .6f        % .7f  % .7f  % .7f ' % (float(bits[0]),float(bits[1]),float(bits[2])\
                    ,float(bits[3]),float(bits[4]),float(bits[5])), file=file_posfor) 
              if 'total drift:' in line:
                is_posfor = False

            if 'TOTAL-FORCE' in line:
              print(line , file=file_posfor , end='')
              is_posfor = True

            # extract energy (eV)
            if 'energy  without entropy=' in line:
              bits=line.split()
              print('% .10f' % (float(bits[6])), file=file_u)
              if is_first_step:
                self.u_lat = float(bits[6])/self.n_atoms
                is_first_step = False




    print('\n','Lattice energy (eV/atom): ', self.u_lat)
    print(' Lattice pressure   (GPa): ', self.p_lat)


  def _get_max_forces(self, forces):

    """ Check if some atoms has nonzero forces

    """
    
    print()
    is_large_force = False
    for i,v in enumerate(forces):
      v_magn = np.linalg.norm(v)
      if v_magn > self.force_tol:
        print(' * WARNING: magnitude of force on atom ',(i+1),' = ', v_magn, ' (> force_tol = ',self.force_tol,')')
        is_large_force = True
    return is_large_force




