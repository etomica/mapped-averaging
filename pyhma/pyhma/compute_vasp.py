import math
import sys , os
import numpy as np
from pyhma.nearest_image import Nearest_Image

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('( (( ((( (((( ((((( (((((( pyhma )))))) ))))) )))) ))) )) )')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

class Compute_VASP:
  kB = 0.000086173306373383  # Boltzmann's constant (eV/K)
  def __init__(self):
    files_list = ['POSCAR_lat','posfor.out','e.out','p.out']
    is_files_found = self.check_files(files_list)
    if not is_files_found:
      sys.exit(1)
    self.file_poscarlat = open('POSCAR_lat', 'r')
    self.file_e = open('e.out', 'r')
    self.file_posfor = open('posfor.out', 'r')
    self.file_energies = open('energies.out', 'w')
    self.read_lattice()
    self.nearest_image = Nearest_Image(self.lat_vecs)
    self.read_pyhma_in()

  # read lattice vecors and basis
  def read_lattice(self):
    self.lat_vecs = np.zeros((3,3))
    for l,line in enumerate(self.file_poscarlat):
      bits = line.split()
      if l >= 2 and l <= 4: 
        self.lat_vecs[l-2,:] = [float(bits[0]),float(bits[1]),float(bits[2])]
      if l == 5:
        self.N = 0
        for i in range(0,len(bits)):
          self.N+=int(bits[i])
        self.r_lat=np.zeros((self.N,3))

      if l >= 7 and l <= 6+self.N:
        n = l-7
        self.r_lat[n] = [float(bits[0]), float(bits[1]), float(bits[2])]

    # close file
    self.file_poscarlat.close()


  # write Conv and HMA
  def compute(self):
    kT_eV = Compute_VASP.kB*self.temperature
    n=-1
    timestep=1
    line_posfor = self.file_posfor.readline()
    while line_posfor:
      if 'POSITION' in line_posfor:
        line_e = self.file_e.readline()
        e = float(line_e)
        if n == -1:
          e_lat = e
        de = e - e_lat 
        n = 0
        fdr = 0
      else:
        bits=line_posfor.split()
        r  = np.array([float(bits[0]) , float(bits[1]) , float(bits[2])])
        dr = r - self.r_lat[n]
        if n == 0:
          dr1 = np.copy(dr)
        dr -= dr1 # reference assigment
        self.nearest_image.get_nearest_image(dr)
        f  = np.array([float(bits[3]) , float(bits[4]) , float(bits[5])])
        fdr = fdr + f.dot(dr) 
        if n == self.N-1:
          uah_conv = (de-1.5*kT_eV*(self.N-1))/self.N
          uah_hma  = (de+0.5*fdr)/self.N
          print('% d  % 0.15f  % 0.15f' % (timestep, uah_conv  , uah_hma) , file=self.file_energies)
          if timestep == 1:
            self.energies=np.zeros((1,3))
            self.energies[0] = [timestep, uah_conv  , uah_hma]
          else:
            self.energies = np.append(self.energies , [[timestep, uah_conv  , uah_hma]] , axis = 0)
          timestep+=1
     
        n += 1
      line_posfor = self.file_posfor.readline()
    self.close_out_files()


  # compute statistics: averages, uncertainties, and autocorrelations
  def statistics(self):
    avg = np.average(self.energies[1])
    sd = np.std(self.energies[1],ddof=1)
    return {'uah_conv': avg, 'uah_sd': sd}
 

  # read pyhma.in file
  def read_pyhma_in(self):
    file_pyhma_in = open('pyhma.in', 'r')
    for line in file_pyhma_in:
      bits=line.split()
      if 'TEMP' in bits:
        self.temperature = float(bits[1])
      if 'DPRESS' in bits:
        self.d_pressure_est = float(bits[1])
      if 'NEQ' in bits:
        self.n_eq = int(bits[1])
      if 'BLOCK' in bits:
        self.blocksize = int(bits[1])


  # check if these files exist (obtained by 'pyhma --vasp --read OUTCAR'):
  # POSCAR_lat , posfor.out , e.out , and p.out
  def check_files(self , files_list):
    is_files_found = True
    for file_i in files_list:
      if not os.path.exists(file_i):
        print('ERROR: ' , file_i ,'file not found.')
        is_files_found = False
    if is_files_found:
      return True


  #close all files
  def close_out_files(self):
    self.file_e.close()
    self.file_posfor.close()
    self.file_energies.close()
   
