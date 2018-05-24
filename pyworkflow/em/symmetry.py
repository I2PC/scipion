# coding=utf-8
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
    SYM_TETRAHEDRAL, SYM_TETRAHEDRAL_Z3, SYM_I222, SYM_I222r, \
    SYM_In25, SYM_In25r
from math import sin, cos, pi, acos, sqrt
from numpy import array, zeros, float, append
from numpy.linalg import inv as matrix_inverse
from numpy import dot as matrix_multiply


def __length(v):
  d = sqrt(sum([e*e for e in v]))
  return d

def __normalize_vector(v):
  d = __length(v)
  if d == 0:
    d = 1
  return tuple([e/d for e in v])

def __rotation_transform(axis, angle, center = (0,0,0)):
    """ Angle is in degrees. """
    axis = __normalize_vector(axis)

    arad = angle*pi/180.0
    sa = sin(arad)
    ca = cos(arad)
    k = 1 - ca
    ax, ay, az = axis
    tf = ((1 + k*(ax*ax-1), -az*sa+k*ax*ay, ay*sa+k*ax*az, 0),
          (az*sa+k*ax*ay, 1 + k*(ay*ay-1), -ax*sa+k*ay*az, 0),
          (-ay*sa+k*ax*az, ax*sa+k*ay*az, 1 + k*(az*az-1), 0))
    cx, cy, cz = center
    c_tf = ((1,0,0,cx), (0,1,0,cy), (0,0,1,cz))
    inv_c_tf = ((1,0,0,-cx), (0,1,0,-cy), (0,0,1,-cz))
    rtf = __multiply_matrices(c_tf, tf, inv_c_tf)
    return rtf

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

  ctfinv = __invert_matrix(ctf)
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
        matrices = __octahedral_symmetry_matrices(center)
    elif sym == SYM_TETRAHEDRAL or sym == SYM_TETRAHEDRAL_Z3:
        matrices = __tetrahedralSymmetryMatrices(sym, center)
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
    """ Rotation about z axis.
    This is a local method. do not access directly to it
    """
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

def __octahedral_symmetry_matrices(center = (0,0,0)):
    """ 4-folds along x, y, z axes. """
    c4 = (((0,0,1),0), ((0,0,1),90), ((0,0,1),180), ((0,0,1),270))
    cube = (((1,0,0),0), ((1,0,0),90), ((1,0,0),180), ((1,0,0),270),
            ((0,1,0),90), ((0,1,0),270))
    c4syms = [__rotation_transform(axis, angle) for axis, angle in c4]
    cubesyms = [__rotation_transform(axis, angle) for axis, angle in cube]
    syms = __matrix_products(cubesyms, c4syms)
    syms = __recenter_symmetries(syms, center)
    return syms

def __dihedralSymmetryMatrices(n, center = (0, 0, 0)):
    """ Rotation about z axis, reflection about x axis. """
    clist = __cyclicSymmetrySatrices(n)
    reflect = ((1,0,0,0),(0,-1,0,0),(0,0,-1,0))
    tflist = __matrix_products([__identity_matrix(), reflect], clist)
    tflist = __recenter_symmetries(tflist, center)
    return tflist


def __tetrahedralSymmetryMatrices(orientation = SYM_TETRAHEDRAL, center = (0,0,0)):
    """
    identity
    4 × rotation by 120 clockwise (seen from a vertex): (234), (143), (412), (321)
    4 × rotation by 120 counterclockwise (ditto)
    3 × rotation by 180
    """
    aa = (((0,0,1),0), ((1,0,0),180), ((0,1,0),180), ((0,0,1),180),
          ((1,1,1),120), ((1,1,1),240), ((-1,-1,1),120), ((-1,-1,1),240),
          ((-1,1,-1),120), ((-1,1,-1),240), ((1,-1,-1),120), ((1,-1,-1),240))
    syms = [__rotation_transform(axis, angle) for axis, angle in aa]

    if orientation == SYM_TETRAHEDRAL_Z3:
        # EMAN convention, 3-fold on z, 3-fold in yz plane along neg y.
        tf = __multiply_matrices(
            __rotation_transform((0,0,1), -45.0),
            __rotation_transform((1,0,0), -acos(1/sqrt(3))*180/pi))
        syms = __coordinate_transform_list(syms, tf)

    syms = __recenter_symmetries(syms, center)
    return syms

def __icosahedralSymmetryMatrices(orientation = SYM_I222, center = (0,0,0)):
    pass