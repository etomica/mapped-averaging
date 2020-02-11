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
This module returns the nearest image of a displacement vector for a given box edge (row) vectors.

"""


import numpy as np

class NearestImage:
  """
  A class to get the nearest image of a displacement (dr), given box edge (raw) vectors (box_row_vecs).

  Parameters
  -----------
  box_row_vecs : numpy.ndarray
    Box edge (row) vectors in Å.

  """  

  def __init__(self, box_row_vecs):
    self.D = 3
    self.halfTol = 0.50000001
    self.tV2 = {}
    self.box_row_vecs = box_row_vecs 
    self.transform_vecs = box_row_vecs
    self.get_transform_vecs(box_row_vecs)

  # Nearest Image 
  def get_nearest_image(self, dr):
    """
    Modify a given displacement (dr) to its nearest image.

    Parameters
    ----------
    dr : numpy.ndarray
      Atomic displacement from lattice site. Will be modified to the nearest image after calling this method. 

    """
    # pretend that we last transformed the last+1 vector
    # this has the effect we want
    lastTransform = len(self.transform_vecs)
    len_transform_vecs = lastTransform
    i = 0
    while i != lastTransform:
      dot = self.transform_vecs[i].dot(dr)/self.tV2[i]
      if abs(dot) > self.halfTol:
        dot = round(dot)
        dr -= dot*self.transform_vecs[i] # MUST do -= to keep the reference/pointer! If you do dr = dr + x, yu only change value and get new reference
        # remember that we transformed with this vector
        # when we finish the vectors, we'll restart back at 0 and continue until we reach this one again.
        lastTransform = i
      if (i == len_transform_vecs-1) and (lastTransform != len_transform_vecs):
        # we finished a pass, but transformed a vector.
        # make another pass (we'll stop at the vector we transformed)
        i = -1
      i += 1


  # Transform Vectors
  def get_transform_vecs(self, box_row_vecs):
    """
    Compute the transform vectors.

    Parameters
    ----------
    box_row_vecs : numpy.ndarray
      Box edge (row) vectors in Å.
 
    """

    tmp1  = 0
    #test edge pairs
    for i in range(self.D-1):
      for k in range(i+1,self.D):
        dot = box_row_vecs[i].dot(box_row_vecs[k])
        if abs(dot) < 1e-10:
          continue
        if dot   < -1e-10:
          tmp1 = box_row_vecs[i] + box_row_vecs[k]
        elif dot >  1e-10:
          tmp1 = box_row_vecs[i] - box_row_vecs[k]

        self.test_transform_vector(tmp1)

    # test edge triplets 
    dot01 = box_row_vecs[0].dot(box_row_vecs[1])
    dot02 = box_row_vecs[0].dot(box_row_vecs[2])
    dot12 = box_row_vecs[1].dot(box_row_vecs[2])

    sum = 0.0

    sum = dot01 + dot02 + dot12 # +0 +1 +2
    if sum < -1e-10:
      tmp1 = box_row_vecs[0] + box_row_vecs[1]
      tmp1 += box_row_vecs[2]
      self.test_transform_vector(tmp1)
            
    sum = dot01 - dot02 - dot12 #+0 +1 -2
    if sum < -1e-10:
      tmp1 = box_row_vecs[0] + box_row_vecs[1]
      tmp1 -= box_row_vecs[2]
      self.test_transform_vector(tmp1)
            
    sum = -dot01 + dot02 - dot12 #+0 -1 +2
    if sum < -1e-10:
      tmp1 = box_row_vecs[0] - box_row_vecs[1]
      tmp1 += box_row_vecs[2]
      self.test_transform_vector(tmp1)

    sum = -dot01 - dot02 + dot12 # +0 -1 -2
    if sum < -1e-10:
      tmp1 = box_row_vecs[0] - box_row_vecs[1]
      tmp1 -= box_row_vecs[2]
      self.test_transform_vector(tmp1)

    for i in range(len(self.transform_vecs)):
      self.tV2[i] = self.transform_vecs[i].dot(self.transform_vecs[i])

  # TEST
  def test_transform_vector(self, v):
    """
    Test transform vectors.

    Parameters
    -----------
    v : numpy.ndarray
      A lattice vector to be tested (Å)

    """

    tmp2 = 0.5*v
    for i in range(len(self.transform_vecs)):
      dot = abs(tmp2.dot(self.transform_vecs[i]))/(self.transform_vecs[i].dot(self.transform_vecs[i]))
      if dot > 1-self.halfTol:
        break
    self.transform_vecs=np.append(self.transform_vecs , [v], axis=0)


