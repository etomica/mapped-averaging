"""

Before ReadVASP class

"""

class ReadVASP:

  """

  After ReadVASP class


  """

  
  def __init__(self, *filenames):  # tuple of filenames
   self.filenames = filenames
    
  def read(self):
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
      for i in range(len(self.filenames)):
        with open(self.filenames[i]) as  file_outcar:
          for line in file_outcar:
            # Create lattice.dat file
            if is_first_step:
              bits = line.split()
              if ' SYSTEM =' in line:
                print('Lattice vectors', file = file_lattice)
                print(' 1.0 scaling factor' , file = file_lattice)
                print(' ',int(n_atoms) , 'atoms (total)')
                print(' Lattice vectors (Angs):')
              if 'number of ions' in line:
                n_atoms = bits[len(bits)-1]

              # lattice
              if is_first_lat:
                if line.strip():
                  print('  % .12f % .12f % .12f' % (float(bits[0]), float(bits[1]), float(bits[2])))
                  print('% .12f % .12f % .12f' % (float(bits[0]), float(bits[1]), float(bits[2]))\
                      ,file = file_lattice)
                else:
                  print('',n_atoms,' atoms (total)',file=file_lattice)
                  print('Basis (Angst)', file=file_lattice)
                  is_lat = False
                  is_first_lat = False
              if ('direct lattice vectors' in line) & is_lat:
                is_first_lat = True 
        
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
                print(' -------------------------------------------------')
                print(' Lattice pressure (GPa): ', p_vir)
              print('%.10f' % p_vir , file=file_p) 

            # extract positions and forces
            if is_posfor:
              if ('----------' not in line) and ('total drift:' not in line):
                bits=line.split()
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
                print(' Lattice energy (eV): ', bits[6])
                print(' -------------------------------------------------')
                is_first_step = False


