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
pyHMA module 

*  ss s s s s s

"""

import numpy as np
import pyhma
from pyhma.nearest_image import NearestImage

class Processor:
  """
  This Processor class is a class to compute ensemble averages.

  Parameters
  ----------
  A : a 
   www
  B : b
   www2


  """
  def __init__(self, data, pressure_qh, meV=False):
    """Retrieve MD simulation data ('data' dict of the read() method) and user-defined parameters ('params' dict).
    """
    self.temperature   = data['temperature']            # set emperature (K)
    self.timestep      = data['timestep']               # MD timestep size (fs)
    self.num_atoms     = data['num_atoms']              # total number of atoms
    self.volume_atom   = data['volume_atom']            # specific volume (A^3/atom)
    self.box_row_vecs  = np.array(data['box_row_vecs']) # box edge (raw) vectors (A)
    self.basis         = np.array(data['basis'])        # atomic positions of initial configuration (fractional)
    self.position     = np.array(data['position'])    # positions at each atom at each MD step (fractional)
    self.force        = np.array(data['force'])       # forces at each atom at each MD step (eV/A)
    self.energy      = np.array(data['energy'])     # energy of each configuration (eV/atom)
    self.pressure_vir = data['pressure_vir']          # virial pressure (i.e., without ideal gas) (GPa)
    self.pressure_qh   = pressure_qh                    # quasiharmonic pressure HMA parameter (GPa)
    self.out_data      = np.empty((0,4))                # anharmonic data array ([e_ah_conv, e_ah_hma, p_ah_conv, p_ah_hma])
    self.meV           = meV
    
  """ compute instantinious properties 
  """
  def process(self, steps_tot=None, verbose=False):
    """compute Conv and HMA anharmonic energy and pressure.
    """
    if steps_tot == None:
      self.steps_tot  = len(self.energy)
    else:
      self.steps_tot  = steps_tot
      if steps_tot > len(self.energy):
        print(' WARNING! User-set steps_tot (', steps_tot,') can not be larger than MD simulation steps (', len(self.energy),').')
        print('          Reduce steps_tot and try again.')
        raise RuntimeError('Illegal total number of steps.')

    kB = 0.0000861733063733830                            # Boltzmann's constant (eV/K)
    eV2J = 1.602176634e-19                                # eV to Joules conversion factor
    kBT_eV = kB*self.temperature                          # eV
    kBT_J  = kBT_eV*eV2J                                  # J
    self.density = 1/self.volume_atom                     # atom/A^3
    self.pressure_ig = self.density*1e30*kBT_J*1e-9       # ideal gas pressure (GPa)
    f_v = (self.pressure_qh*1e9/kBT_J-self.density*1e30)\
         /(3*(self.num_atoms-1))*1e-9                     # f_v variable HMA pressure (GPa/J)
    energy_lat   = self.energy[0]
    pressure_lat = self.pressure_vir[0]
    if verbose:
      print('\nSimulation data')
      print('===============')
      print(' Set temperature       (K): %10.5f' % self.temperature)
      print(' Volume         (A^3/atom): %10.5f' % self.volume_atom)
      print(' MD timestep          (fs): %10.5f' % self.timestep)
      print(' Lattice energy  (eV/atom): %10.5f' % energy_lat)
      print(' Harmonic energy (eV/atom): %10.5f' % (1.5*kBT_eV*(self.num_atoms-1)/self.num_atoms))
      print(' Lattice pressure    (GPa): %10.5f' % pressure_lat)
      print(' Harmonic pressure   (GPa): %10.5f' % self.pressure_qh)
      print('\n Found', len(self.energy) ,' total MD steps')
      print(' Using', self.steps_tot, ' user-set MD steps')
      print('\n Computing instantaneous properties ...')
 
    basis_cart = Processor._direct_to_cart(self.basis, self.box_row_vecs)
    e_fac = 1.0
    if self.meV:
      e_fac = 1.0e3
    with open('energy_ah.out','w') as file_energy_ah, open('pressure_ah.out','w') as file_pressure_ah:
      self._nearest_image = NearestImage(self.box_row_vecs) 
      for step in range(self.steps_tot): # snaps
        sim_time = step*self.timestep
        # Conv
        e_ah_conv = e_fac*(self.energy[step] - energy_lat - 1.5*kBT_eV*(self.num_atoms-1)/self.num_atoms)
        p_ah_conv = self.pressure_vir[step] + self.pressure_ig - pressure_lat - self.pressure_qh
        # HMA
        r_cart  = Processor._direct_to_cart(self.position[step], self.box_row_vecs)
        fdr = 0
        for atom in range(self.num_atoms): # atoms
          f  = self.force[step][atom]
          dr = r_cart[atom] - basis_cart[atom]
          if atom == 0:
            dr1 = np.copy(dr)
          dr -= dr1 # reference assigment
          self._nearest_image.get_nearest_image(dr)
          fdr = fdr + f.dot(dr)

        e_ah_hma  = e_fac*(self.energy[step] + 0.5*fdr/self.num_atoms - energy_lat)
        p_ah_hma  = self.pressure_vir[step] + f_v*fdr*eV2J - pressure_lat
        print('%10.1f  %10.5f  %10.5f' % (sim_time, e_ah_conv, e_ah_hma) , file=file_energy_ah)
        print('%10.1f  %10.5f  %10.5f' % (sim_time, p_ah_conv, p_ah_hma) , file=file_pressure_ah)
        self.out_data = np.append(self.out_data , [[e_ah_conv, e_ah_hma, p_ah_conv, p_ah_hma]] , axis=0)

  # Compute statistics: average (avg), stochastic uncertainty (err), and correlation (cor)
  def get_stats(self, steps_eq, blocksize, verbose=False):
    """get_stats() function!!!
    """
    if self.steps_tot < steps_eq:
      print('WARNING! Number of equilibaration steps (', steps_eq,') can not be larger than total steps (', self.steps_tot,').')
      print('         Reduce steps_eq and try again.')
      raise RuntimeError('Illegal equilibaration steps.')
    elif int((self.steps_tot - steps_eq)/blocksize) < 2:
      print('WARNING! Number of blocks ((steps_tot-steps_eq)/blocksize) must be at least two.')
      print('         Reduce blocksize to get finite number of blocks and try again.')
      raise RuntimeError('Illegal block size.') 

    n_prod = self.steps_tot - steps_eq
    data_prod = self.out_data[steps_eq:self.steps_tot,:]

    if verbose:
      print('\nBlock averaging statistics')
      print('==========================')
      print('', n_prod, 'production steps (after', steps_eq ,'equilibration steps)')
      print('',int((self.steps_tot-steps_eq)/blocksize), 'blocks (blocksize =' ,  blocksize,' steps)\n')
      print(' Computing statistics ...')

    data_prod_block = Processor._block_data(data_prod, blocksize) 
    # Averages
    data_avg = np.average(data_prod, axis=0) # used raw data because average should be independent on blocks
    # Uncertainty
    data_err = np.std(data_prod_block, ddof = 1, axis=0)/np.sqrt(len(data_prod_block))
    # Correlation
    data_cor = Processor._get_cor(data_prod_block)

    self.out_stats = {'e_ah_conv': {'avg': data_avg[0] , 'err': data_err[0] , 'cor': data_cor[0]},\
                      'e_ah_hma' : {'avg': data_avg[1] , 'err': data_err[1] , 'cor': data_cor[1]}, \
                      'p_ah_conv': {'avg': data_avg[2] , 'err': data_err[2] , 'cor': data_cor[2]}, \
                      'p_ah_hma' : {'avg': data_avg[3] , 'err': data_err[3] , 'cor': data_cor[3]}}
    return self.out_stats


  """Print statistics (avg, err, cor)
  """
  def print_stats(self, out_stats):
    """ PRINT out_stats
    """
    e_units='eV'
    if self.meV:
      e_units='meV'

    print('\n e_ah_conv (%s/atom): %10.5f +/- %5.1e    cor: %4.2f' % (e_units, out_stats['e_ah_conv']['avg'],\
          out_stats['e_ah_conv']['err'], out_stats['e_ah_conv']['cor']))
    print(' e_ah_hma  (%s/atom): %10.5f +/- %5.1e    cor: %4.2f' % (e_units, out_stats['e_ah_hma']['avg'],\
          out_stats['e_ah_hma']['err'], out_stats['e_ah_hma']['cor']))
    print(' p_ah_conv      (GPa): %10.5f +/- %5.1e    cor: %4.2f' % (out_stats['p_ah_conv']['avg'],\
           out_stats['p_ah_conv']['err'], out_stats['p_ah_conv']['cor']))
    print(' p_ah_hma       (GPa): %10.5f +/- %5.1e    cor: %4.2f\n' % (out_stats['p_ah_hma']['avg'],\
           out_stats['p_ah_hma']['err'], out_stats['p_ah_hma']['cor']),end='')





  """Block data to blocks of size blocksize
  """
  @staticmethod
  def _block_data(data, blocksize):
    """ block data!!!
    """
    sum = np.zeros((1,len(data[0])))
    n = 1  # block number
    for i in range(len(data)):
      sum += data[i]
      if (i+1) % blocksize == 0:
        block_avg =  sum/blocksize
        if n == 1: 
          data_block = block_avg
        else:
          data_block = np.append(data_block, block_avg , axis = 0)
        sum = np.zeros((1,len(data[0])))
        n+=1
    return data_block


  """Compute autocorrelation between adjacent data points
  """
  @staticmethod
  def _get_cor(data):  # without the N-1 correction in both numerators and denominators
    """ AC of data!!!
    """
    data_avg = np.average(data, axis=0)
    data_var = np.var(data, ddof = 0, axis=0) # using N, not N-1
    sum = np.zeros(len(data[0]))
    for i in range(len(data)-1): # 0,1,... n-2
      sum += (data[i]-data_avg)*(data[i+1]-data_avg)
    sum *= 1/(len(data)-1)
    cor = sum/data_var
    return cor


  """Convert direct (fractional) to Cartesian coordinates, for a given box vectors a.
  """
  @staticmethod
  def _direct_to_cart(x, a):
    r = np.empty((0,3))
    for i in range(len(x)): # loop over atoms
      r_i = np.transpose(a).dot(x[i])
      r = np.append(r, [r_i], axis=0)
    return r  
 
