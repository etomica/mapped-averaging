import math
import sys, os
import numpy as np
from nearest_image import NearestImage

kB = 0.000086173306373383  # Boltzmann's constant (eV/K)
eV2J = 1.60217733e-19
data_size = 2  # u & p

"""
Simulation is a class to compute ensemble averages from outputs of step #1.

"""

class Simulation:

  def __init__(self):
    files_list = ['lattice.dat','posfor.dat','u.dat','p.dat','pyhma.in']
    self.check_files_found(files_list)
    # read files
    self.read_poscar_lat()
    self.read_pyhma_in()

  def read_poscar_lat(self):
    self.lat_vecs = np.zeros((3,3))
    with open('lattice.dat') as file_poscar_lat:
      for l,line in enumerate(file_poscar_lat):
        bits = line.split()
        if l >= 2 and l <= 4: 
          self.lat_vecs[l-2,:] = [float(bits[0]),float(bits[1]),float(bits[2])]
        if l == 5:
          self.n_atoms = int(bits[0])
          self.r_lat=np.zeros((self.n_atoms,3))
        if l >= 7 and l <= 6+self.n_atoms:
          n = l-7
          self.r_lat[n] = [float(bits[0]), float(bits[1]), float(bits[2])]
  
    self.nearest_image = NearestImage(self.lat_vecs)


# read pyhma.in file
  def read_pyhma_in(self):
    with open('pyhma.in') as file_pyhma: 
      for line in file_pyhma:
        bits=line.split()
        if 'TEMP' in line:
          self.temperature = float(bits[1]) # K
        if 'PQH' in line:
          self.p_qh = float(bits[1]) # GPa
        if 'NEQ' in line:
          self.n_eq = int(bits[1])
        if 'BLOCK' in line:
          self.blocksize = int(bits[1])




  # write Conv and HMA
  def run(self):
    kBT_eV = kB*self.temperature # eV
    kBT_J  = kBT_eV*eV2J         # J
    tmp_vec = np.cross(self.lat_vecs[0],self.lat_vecs[1])
    self.volume = np.dot(tmp_vec , self.lat_vecs[2])/self.n_atoms # Ang^3/atom
    self.density = 1.0/self.volume   # 1/m^3
    self.p_ig = self.density*1e30*kBT_J*1e-9  # GPa

    f_v = (self.p_qh*1e9/kBT_J-self.density*1e30)/(3*(self.n_atoms-1))*1e-9  # GPa/J
    is_first_step = True

    print(' Run ....')
    with open('posfor.dat') as file_posfor, open('u.dat') as file_u, open('p.dat') as file_p,\
        open('uah.out','w') as file_uah, open('pah.out','w') as file_pah:
      for l,line in enumerate(file_posfor):
        n = (l % int(self.n_atoms+1))-1  # n=-1,0,1,..N-1
        if n == -1: # beginning of each snap
          u = float(file_u.readline())/self.n_atoms
          p_vir = float(file_p.readline())  # p_vir = p_conv - p_ig
          if is_first_step:
            self.u_lat = u
            self.p_lat = p_vir
            is_first_step = False
            self.md_steps = 1
          fdr = 0
        else: # read each snap data
          bits = line.split()
          r  = np.array([float(bits[0]) , float(bits[1]) , float(bits[2])])
          dr = r - self.r_lat[n]
          if n == 0:
            dr1 = np.copy(dr)
          dr -= dr1 # reference assigment
          self.nearest_image.get_nearest_image(dr)
          f  = np.array([float(bits[3]) , float(bits[4]) , float(bits[5])])
          fdr = fdr + f.dot(dr) 
          if n == self.n_atoms-1: # last atom
            # anharmonic energy
            uah_conv = u-self.u_lat - 1.5*kBT_eV*(self.n_atoms-1)/self.n_atoms
            uah_hma  = u-self.u_lat + 0.5*fdr/self.n_atoms
            print('% d  % 0.15f  % 0.15f' % (self.md_steps, uah_conv  , uah_hma) , file=file_uah)
            
            # pressure (off p_lat)
            pah_conv = p_vir + self.p_ig - self.p_lat - self.p_qh
            pah_hma  = p_vir + f_v*fdr*eV2J - self.p_lat
            print('% d  % 0.15f  % 0.15f' % (self.md_steps, pah_conv, pah_hma) , file=file_pah)
            if self.md_steps == 1:
              self.data = np.zeros((1,data_size*2))
              self.data[0] = [uah_conv, uah_hma, pah_conv, pah_hma]
            else:
              self.data = np.append(self.data , [[uah_conv  , uah_hma, pah_conv, pah_hma]] , axis = 0)
            self.md_steps+=1

    print(' Done. Found', (self.md_steps-1-self.n_eq) , 'steps after', self.n_eq ,'steps of equilibaration\n')

    print(' Temperature           (K):' , self.temperature)
    print(' Volume       (Ang^3/atom):', self.volume)
    print(' Lattice energy  (eV/atom):', self.u_lat)
    print(' Harmonic energy (eV/atom):', 1.5*kBT_eV*(self.n_atoms-1)/self.n_atoms)
    print(' Lattice pressure    (GPa):', self.p_lat)
    print(' Harmonic pressure   (GPa):', self.p_qh,'\n')



  # compute statistics: averages, uncertainties, and autocorrelations
  def get_statistics(self):
    print(' Block averaging statistics')
    print(' ==========================')
    print(' Using ',int((self.md_steps-1-self.n_eq)/self.blocksize), ' blocks of size ' ,  self.blocksize,' MD steps','\n')
    data_prod = np.copy(self.data[self.n_eq:,:])
    data_prod_block = self.block_data(data_prod) 
    # Averages
    data_avg = np.average(data_prod, axis=0) # used raw data because average should be independent on blocks
    # Uncertainty
    data_err = np.std(data_prod_block, ddof = 1, axis=0)/np.sqrt(len(data_prod_block))
    # Correlation
    data_ac = self.get_ac(data_prod_block)

    return {'uahc': {'avg': data_avg[0] , 'err': data_err[0] , 'corr': data_ac[0]},\
            'uahm': {'avg': data_avg[1] , 'err': data_err[1] , 'corr': data_ac[1]}, \
            'pahc': {'avg': data_avg[2] , 'err': data_err[2] , 'corr': data_ac[2]}, \
            'pahm': {'avg': data_avg[3] , 'err': data_err[3] , 'corr': data_ac[3]}}


  def block_data(self, data):
    sum = np.zeros((1,len(data[0])))
    n = 1  # block number
    for i in range(len(data)):
      sum += data[i]
      if (i+1) % self.blocksize == 0:
        block_avg =  sum/self.blocksize
        if n == 1: 
          data_block = block_avg
        else:
          data_block = np.append(data_block, block_avg , axis = 0)
        sum = np.zeros((1,len(data[0])))
        n+=1
    return data_block


  def get_ac(self, data):  # without the N-1 correction in both numerators and denimrators
    data_avg = np.average(data, axis=0)
    data_var = np.var(data, ddof = 0, axis=0) # using N, not N-1
    sum = np.zeros(len(data[0]))
    for i in range(len(data)-1): # 0,1,... n-2
      sum += (data[i]-data_avg)*(data[i+1]-data_avg)
    sum *= 1/(len(data)-1)
    ac = sum/data_var
    return ac



  # check if these files exis: lattice.dat, posfor.dat , u.dat , and p.dat
  def check_files_found(self , files_list):
    is_file_found = True
    for file_i in files_list:
      if not os.path.exists(file_i):
        print(' ERROR: ' , file_i ,'not found.')
        is_file_found = False
    if not is_file_found:
      sys.exit(1)
