import sys


class Read_VASP:
  def __init__(self, *filenames):
    self.file_outcars = []
    self.file_e = open('e.out', 'w')
    self.file_p = open('p.out', 'w')
    self.file_posfor = open('posfor.out', 'w')
    self.file_poscarlat = open('POSCAR_lat', 'w')
    self.filenames = filenames
    
  def read(self):
    for i in range(len(self.filenames)):
      self.file_outcars.append(open(self.filenames[i], 'r'))
      start     = -1
      is_first  = True
      is_lat    = False
      is_basis  = False
      is_posfor = False

      for line in self.file_outcars[i]:
        # create POSCAR_lat file
        if i == 0 and is_first:
          bits = line.split()
    
          if ' NIONS =' in line: N = bits[len(bits)-1]
          if ' SYSTEM =' in line:
            print(bits[2], file = self.file_poscarlat)
            print(' 1.0', file = self.file_poscarlat)
    
          # lattice
          if is_lat:
            if line.strip():
              print('% .12f % .12f % .12f' % (float(bits[0]), float(bits[1]), float(bits[2]))\
                  ,file = self.file_poscarlat)
            else:
              print('',N,file=self.file_poscarlat)
              print('cartesian coordinates  (Angst)', file=self.file_poscarlat)
              is_lat = False
          if 'direct lattice vectors' in line: is_lat = True 
    
          # basis
          if is_basis:
            if line.strip():
              print('% .12f % .12f % .12f' % (float(bits[0]), float(bits[1]), float(bits[2]))\
                  ,file = self.file_poscarlat)
            else:
              is_basis = False
              is_first = False
          if 'position of ions in cartesian coordinates' in line: is_basis = True
    
        # extract E0 (eV)
        if 'y=' in line:
          bits=line.split()
          print('% .10f' % (float(bits[6])), file=self.file_e)

        # extract P (GPa)
        if 'external' in line:
          bits=line.split()
          print('%.5f' % (float(bits[3])/10), file=self.file_p) 

        # extract positions and forces
        if is_posfor:
          if ('---' not in line) and ('total drift:' not in line):
            bits=line.split()
            print('% .6f  % .6f  % .6f  % .7f  % .7f  % .7f ' % (float(bits[0]),float(bits[1]),float(bits[2])\
                ,float(bits[3]),float(bits[4]),float(bits[5])), file=self.file_posfor) 
          if 'total drift:' in line:
            is_posfor = False

        if 'TOTAL-FORCE' in line:
          print(line , file=self.file_posfor , end='')
          is_posfor = True
          
      if i == len(self.filenames)-1:
        self.close_all_files()


  def close_all_files(self):
    self.file_poscarlat.close()
    self.file_posfor.close()
    self.file_e.close()
    self.file_p.close()
    for j in range(len(self.filenames)):
      self.file_outcars[j].close()
