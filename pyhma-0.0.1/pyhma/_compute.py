import math
import sys, os
import numpy as np
from nearest_image import NearestImage

kB = 0.000086173306373383  # Boltzmann's constant (eV/K)
eV2J = 1.60217733e-19
data_size = 2  # u & p

"""

COMPUTE AVERAGES

"""

class Compute:
  """

  COMPUTE AVERAGES

  """

  def __init__(self):
    files_list = ['POSCAR_LAT','posfor.dat','u.dat','p.dat','pyhma.in']
    self.check_files_found(files_list)
    # read files
    self.read_poscar_lat()
    self.read_pyhma_in()
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('( (( ((( (((( ((((( (((((( PYHMA )))))) ))))) )))) ))) )) )')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


  def read_poscar_lat(self):
    self.lat_vecs = np.zeros((3,3))
    with open('POSCAR_LAT') as file_poscar_lat:
      for l,line in enumerate(file_poscar_lat):
        bits = line.split()
        if l >= 2 and l <= 4: 
          self.lat_vecs[l-2,:] = [float(bits[0]),float(bits[1]),float(bits[2])]
        if l == 5:
          self.n_atoms = 0
          for i in range(0,len(bits)):
            self.n_atoms+=int(bits[i])
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
        if 'BETADP' in line:
          self.dp_T_est = float(bits[1]) # GPa
        if 'NEQ' in line:
          self.n_eq = int(bits[1])
        if 'BLOCK' in line:
          self.blocksize = int(bits[1])



  # write Conv and HMA
  def run(self):
    kBT_eV = kB*self.temperature # eV
    kBT_J  = kBT_eV*eV2J         # J
    tmp_vec = np.cross(self.lat_vecs[0],self.lat_vecs[1])
    volume = np.dot(tmp_vec , self.lat_vecs[2]) # Ang^3
    density = self.n_atoms/volume*1e30    # 1/m^3
    p_ig = density*kBT_J*1e-9  # GPa
    dp_est = self.dp_T_est*self.temperature
    f_v = (dp_est*1e9/kBT_J-density)/(3*(self.n_atoms-1))*1e-9  # GPa/J
    is_first_step = True
    with open('posfor.dat') as file_posfor, open('u.dat') as file_u, open('p.dat') as file_p, open('uah.out','w') as file_uah, open('dp.out','w') as file_dp:
      for l,line in enumerate(file_posfor):
        n = (l % int(self.n_atoms+1))-1  # n=-1,0,1,..N-1
        if n == -1: # beginning of each snap
          e = float(file_u.readline())
          p_vir = float(file_p.readline())  # p_vir = p_conv - p_ig
          if is_first_step:
            e_lat = e
            p_lat = p_vir
            is_first_step = False
            timestep = 1
          fdr = 0
        else: # read each snap data
          bits = line.split()
          r  = np.array([float(bits[0]) , float(bits[1]) , float(bits[2])])
          dr = r - self.r_lat[n]
          if n == 0:
            dr1 = np.copy(dr)
          dr -= dr1 # reference assigment
          self.nearest_image.get_nearest_image(dr)
          r = self.r_lat[n] + dr
          f  = np.array([float(bits[3]) , float(bits[4]) , float(bits[5])])
          #fdr += f.dot(dr)
          fdr += f.dot(r)
          #fdr += f.dot(self.r_lat[n])
          #fdr += f[0]
          if n == self.n_atoms-1: # last atom
            # anharmonic energy
            uah_conv = (e-e_lat - 1.5*kBT_eV*(self.n_atoms-1))/self.n_atoms
            #uah_hma  = (e-e_lat + 0.5*fdr)/self.n_atoms

            uah_hma  =  fdr
            print('% d  % 0.15f  % 0.15f' % (timestep, uah_conv  , uah_hma) , file=file_uah)
            
            # pressure (off p_lat)
            dp_conv = p_vir + p_ig - p_lat
            dp_hma = dp_est + p_vir + f_v*fdr*eV2J - p_lat
            print('% d  % 0.15f  % 0.15f' % (timestep, dp_conv, dp_hma) , file=file_dp)
            if timestep == 1:
              self.data = np.zeros((1,data_size*2))
              self.data[0] = [uah_conv, uah_hma, dp_conv, dp_hma]
            else:
              self.data = np.append(self.data , [[uah_conv  , uah_hma, dp_conv, dp_hma]] , axis = 0)
            timestep+=1


  # compute statistics: averages, uncertainties, and autocorrelations
  def get_statistics(self):
    data_prod = np.copy(self.data[self.n_eq:,:])
    data_prod_block = self.block_data(data_prod) 
    # Averages
    data_avg = np.average(data_prod, axis=0) # used raw data because average should be independent on blocks
    # Uncertainty
    data_err = np.std(data_prod_block, ddof = 1, axis=0)/np.sqrt(len(data_prod_block))
    # Correlation
    data_ac = self.get_ac(data_prod_block)

    return {'uahc': [data_avg[0], data_err[0],data_ac[0]], 'uahm': [data_avg[1], data_err[1],data_ac[1]], 
            'dpc': [data_avg[2], data_err[2],data_ac[2]], 'dpm':  [data_avg[3], data_err[3],data_ac[3]]}




  
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




  # check if these files exist (obtained by 'pyhma --vasp --read OUTCAR'):
  # POSCAR_LAT , posfor.dat , u.dat , and p.dat
  def check_files_found(self , files_list):
    is_file_found = True
    for file_i in files_list:
      if not os.path.exists(file_i):
        print(' ERROR: ' , file_i ,'not found.')
        is_file_found = False
    if not is_file_found:
      sys.exit(1)
