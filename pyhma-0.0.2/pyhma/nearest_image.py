###############################################################################
# pyhma: A Python library for HMA method 
# 
# Copyright (c) 2019
# 
# Authors: Sabry Moustafa, Andrew Schultz, and David Kofke 
# 
# pyhma is a free software
###############################################################################

"""
This module returns the nearest image of a displacement vector for a given box edge vector.
"""


import numpy as np

class NearestImage:
  def __init__(self, latVecs):
    self.D = 3
    self.halfTol = 0.50000001
    self.tV2 = {}
    self.latVecs = latVecs 
    self.transform_vecs = latVecs
    self.get_transform_vecs(latVecs)

  # Nearest Image 
  def get_nearest_image(self, dr):
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
  def get_transform_vecs(self, latVecs):
    tmp1  = 0
    #test edge pairs
    for i in range(self.D-1):
      for k in range(i+1,self.D):
        dot = latVecs[i].dot(latVecs[k])
        if abs(dot) < 1e-10:
          continue
        if dot   < -1e-10:
          tmp1 = latVecs[i] + latVecs[k]
        elif dot >  1e-10:
          tmp1 = latVecs[i] - latVecs[k]

        self.test_transform_vector(tmp1)

    # test edge triplets 
    dot01 = latVecs[0].dot(latVecs[1])
    dot02 = latVecs[0].dot(latVecs[2])
    dot12 = latVecs[1].dot(latVecs[2])

    sum = 0.0

    sum = dot01 + dot02 + dot12 # +0 +1 +2
    if sum < -1e-10:
      tmp1 = latVecs[0] + latVecs[1]
      tmp1 += latVecs[2]
      self.test_transform_vector(tmp1)
            
    sum = dot01 - dot02 - dot12 #+0 +1 -2
    if sum < -1e-10:
      tmp1 = latVecs[0] + latVecs[1]
      tmp1 -= latVecs[2]
      self.test_transform_vector(tmp1)
            
    sum = -dot01 + dot02 - dot12 #+0 -1 +2
    if sum < -1e-10:
      tmp1 = latVecs[0] - latVecs[1]
      tmp1 += latVecs[2]
      self.test_transform_vector(tmp1)

    sum = -dot01 - dot02 + dot12 # +0 -1 -2
    if sum < -1e-10:
      tmp1 = latVecs[0] - latVecs[1]
      tmp1 -= latVecs[2]
      self.test_transform_vector(tmp1)

    for i in range(len(self.transform_vecs)):
      self.tV2[i] = self.transform_vecs[i].dot(self.transform_vecs[i])

  # TEST
  def test_transform_vector(self, v):
    tmp2 = 0.5*v
    for i in range(len(self.transform_vecs)):
      dot = abs(tmp2.dot(self.transform_vecs[i]))/(self.transform_vecs[i].dot(self.transform_vecs[i]))
      if dot > 1-self.halfTol:
        break
    self.transform_vecs=np.append(self.transform_vecs , [v], axis=0)


