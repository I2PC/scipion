# **************************************************************************
# *
# * Authors:     Roberto Marabini
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This module returns the matrices related with the different
point symmetries. Code based on chiomera file sym.py
"""
from constants import SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL, \
    SYM_TETRAHEDRAL, SYM_I222, SYM_I222r, SYM_In25, SYM_In25r
from math import sin, cos, pi
from numpy import array, zeros, float, append
from numpy.linalg import inv as matrix_inverse
from numpy import dot as matrix_multiply


def __translation_matrix(shift):

  tf = array(((1.0, 0, 0, shift[0]),
              (0, 1.0, 0, shift[1]),
              (0, 0, 1.0, shift[2])))
  return tf

def __identity_matrix():

  return ((1.0,0,0,0), (0,1.0,0,0), (0,0,1.0,0))

def __invert_matrix(tf):

    tf = array(tf)
    r = tf[:,:3]
    t = tf[:,3]
    tfinv = zeros((3,4), float)
    rinv = tfinv[:,:3]
    tinv = tfinv[:,3]
    rinv[:,:] = matrix_inverse(r)
    tinv[:] = matrix_multiply(rinv, -t)
    return tfinv

def __multiply_matrices(*mlist):

  if len(mlist) == 2:
    m1, m2 = mlist
    p = [[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0]]
    for r in range(3):
      for c in range(4):
        p[r][c] = m1[r][0]*m2[0][c] + m1[r][1]*m2[1][c] + m1[r][2]*m2[2][c]
        p[r][3] += m1[r][3]
    p = tuple(map(tuple, p))
  else:
    p = __multiply_matrices(*mlist[1:])
    p = __multiply_matrices(mlist[0], p)
  return p

def __matrix_products(mlist1, mlist2):
  plist = []
  for m1 in mlist1:
    for m2 in mlist2:
      m1xm2 = __multiply_matrices(m1, m2)
      plist.append(m1xm2)
  return plist

def __coordinate_transform_list(tflist, ctf):

  ctfinv = invert_matrix(ctf)
  return [__multiply_matrices(ctfinv, tf, ctf) for tf in tflist]


def __recenter_symmetries(tflist, center):

    if center == (0,0,0):
      return tflist
    ctf = __translation_matrix([-x for x in center])
    return __coordinate_transform_list(tflist, ctf)

def getSymmetryMatrices(sym=SYM_CYCLIC, n=1, center = (0,0,0)):
    """ interface between scipion and chimera code
        chimera code uses tuples of tuples as matrices
        but scipion uses np.arrays (lists of lists)
        so let us convert them here
    """
    if sym == SYM_CYCLIC:
        matrices = __cyclicSymmetrySatrices(n, center)
    elif sym == SYM_DIHEDRAL:
        matrices = __dihedralSymmetryMatrices(n, center)
    elif sym == SYM_OCTAHEDRAL:
        pass
    elif sym == SYM_TETRAHEDRAL:
        pass
    elif sym == SYM_I222:
        pass
    elif sym == SYM_I222r:
        pass
    elif sym == SYM_In25:
        pass
    elif sym == SYM_In25r:
        pass

    # convert from 4x 3 to 4x4 matrix, Scipion standard
    extraRow = (0., 0., 0., 1.)
    for i in range(len(matrices)):
        matrices[i] +=  (extraRow,)
    # convert from sets to lists Scipion standard
    return array(matrices)

def __cyclicSymmetrySatrices(n, center = (0, 0, 0)):
    """ Rotation about z axis."""
    tflist = []
    for k in range(n):
        a = 2*pi * float(k) / n
        c = cos(a)
        s = sin(a)
        tf = ((c, -s, 0, 0),
              (s, c, 0, 0),
              (0,0,1,0))
        tflist.append(tf)
    tflist = __recenter_symmetries(tflist, center)
    return tflist

def __dihedralSymmetryMatrices(n, center = (0, 0, 0)):
    """ Rotation about z axis, reflection about x axis. """
    clist = __cyclicSymmetrySatrices(n)
    reflect = ((1,0,0,0),(0,-1,0,0),(0,0,-1,0))
    tflist = __matrix_products([__identity_matrix(), reflect], clist)
    tflist = __recenter_symmetries(tflist, center)
    return tflist


