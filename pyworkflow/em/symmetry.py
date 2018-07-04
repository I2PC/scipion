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
from math import sin, cos, pi, acos, sqrt, asin
from numpy import array, zeros, float, append, transpose, add
from numpy.linalg import inv as matrix_inverse
from numpy import dot as matrix_multiply
import operator

def _applyMatrix(tf, points):
    """
    Args:multiply point by a matrice list """

    tf = array(tf)
    r = matrix_multiply(points, transpose(tf[:, :3]))
    add(r, tf[:, 3], r)
    return r

def __length(v):
  d = sqrt(sum([e*e for e in v]))
  return d

def _normalizeVector(v):
  d = __length(v)
  if d == 0:
    d = 1
  return tuple([e/d for e in v])

def _rotationTransform(axis, angle, center = (0, 0, 0)):
    """ Angle is in degrees. """
    axis = _normalizeVector(axis)

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
    rtf = _multiplyMatrices(c_tf, tf, inv_c_tf)
    return rtf

def _translationMatrix(shift):

  tf = array(((1.0, 0, 0, shift[0]),
              (0, 1.0, 0, shift[1]),
              (0, 0, 1.0, shift[2])))
  return tf

def _identityMatrix():

  return ((1.0,0,0,0), (0,1.0,0,0), (0,0,1.0,0))

def _invertMatrix(tf):

    tf = array(tf)
    r = tf[:,:3]
    t = tf[:,3]
    tfinv = zeros((3,4), float)
    rinv = tfinv[:,:3]
    tinv = tfinv[:,3]
    rinv[:,:] = matrix_inverse(r)
    tinv[:] = matrix_multiply(rinv, -t)
    return tfinv

def _multiplyMatrices(*mlist):

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
    p = _multiplyMatrices(*mlist[1:])
    p = _multiplyMatrices(mlist[0], p)
  return p

def _matrixProducts(mlist1, mlist2):
  plist = []
  for m1 in mlist1:
    for m2 in mlist2:
      m1xm2 = _multiplyMatrices(m1, m2)
      plist.append(m1xm2)
  return plist

def _coordinateTransformList(tflist, ctf):

  ctfinv = _invertMatrix(ctf)
  return [_multiplyMatrices(ctfinv, tf, ctf) for tf in tflist]


def _recenterSymmetries(tflist, center):

    if center == (0,0,0):
      return tflist
    ctf = _translationMatrix([-x for x in center])
    return _coordinateTransformList(tflist, ctf)

def _transposeMatrix(tf):

  return ((tf[0][0], tf[1][0], tf[2][0], tf[0][3]),
          (tf[0][1], tf[1][1], tf[2][1], tf[1][3]),
          (tf[0][2], tf[1][2], tf[2][2], tf[2][3]))

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
        matrices = __octahedralSymmetryMatrices(center)
    elif sym == SYM_TETRAHEDRAL or sym == SYM_TETRAHEDRAL_Z3:
        matrices = __tetrahedralSymmetryMatrices(sym, center)
    elif sym == SYM_I222 or sym == SYM_I222r or \
        sym == SYM_In25 or sym == SYM_In25r:
        matrices = __icosahedralSymmetryMatrices(sym, center)

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
    tflist = _recenterSymmetries(tflist, center)
    return tflist

def __octahedralSymmetryMatrices(center = (0, 0, 0)):
    """ 4-folds along x, y, z axes. """
    c4 = (((0,0,1),0), ((0,0,1),90), ((0,0,1),180), ((0,0,1),270))
    cube = (((1,0,0),0), ((1,0,0),90), ((1,0,0),180), ((1,0,0),270),
            ((0,1,0),90), ((0,1,0),270))
    c4syms = [_rotationTransform(axis, angle) for axis, angle in c4]
    cubesyms = [_rotationTransform(axis, angle) for axis, angle in cube]
    syms = _matrixProducts(cubesyms, c4syms)
    syms = _recenterSymmetries(syms, center)
    return syms

def __dihedralSymmetryMatrices(n, center = (0, 0, 0)):
    """ Rotation about z axis, reflection about x axis. """
    clist = __cyclicSymmetrySatrices(n)
    reflect = ((1,0,0,0),(0,-1,0,0),(0,0,-1,0))
    tflist = _matrixProducts([_identityMatrix(), reflect], clist)
    tflist = _recenterSymmetries(tflist, center)
    return tflist


def __tetrahedralSymmetryMatrices(orientation = SYM_TETRAHEDRAL,
                                  center = (0,0,0)):
    """
    identity
    4 * rotation by 120 clockwise (seen from a vertex):
                           (234), (143), (412), (321)
    4 * rotation by 120 counterclockwise (ditto)
    3 * rotation by 180
    """
    aa = (((0,0,1),0), ((1,0,0),180), ((0,1,0),180), ((0,0,1),180),
          ((1,1,1),120), ((1,1,1),240), ((-1,-1,1),120), ((-1,-1,1),240),
          ((-1,1,-1),120), ((-1,1,-1),240), ((1,-1,-1),120),
          ((1,-1,-1),240))
    syms = [_rotationTransform(axis, angle) for axis, angle in aa]

    if orientation == SYM_TETRAHEDRAL_Z3:
        # EMAN convention, 3-fold on z, 3-fold in yz plane along neg y.
        tf = _multiplyMatrices(
            _rotationTransform((0, 0, 1), -45.0),
            _rotationTransform((1, 0, 0), -acos(1 / sqrt(3)) * 180 / pi))
        syms = _coordinateTransformList(syms, tf)

    syms = _recenterSymmetries(syms, center)
    return syms


def __icosahedralSymmetryMatrices(orientation = SYM_I222, center = (0,0,0)):
    if orientation == SYM_I222:
        sym = '222'
    elif orientation == SYM_I222r:
        sym = '222r'
    elif orientation == SYM_In25:
        sym = 'n25'
    elif orientation == SYM_In25r:
        sym = 'n25r'

    i = Icosahedron(orientation=sym, center=center)
    return list(i.icosahedralSymmetryMatrices())

icos_matrices = {}  # Maps orientation name to 60 matrices.
class Icosahedron(object):

    def __init__(self, circumscribed_radius = 1, orientation = '222', center=(
            0, 0, 0)):
        """point = np.array([0, 1, PHI]"""
        self.circumscribed_radius = circumscribed_radius
        self.orientation = orientation
        self.center=center
        # Triangle edge length of unit icosahedron.
        self.e= sqrt(2 - 2 / sqrt(5))
        self.vertices = None  # icosahedron vertices
        self.triangles = None # icosahedron faces
        self._3foldAxis = [] #
        self._2foldAxis = [] #
        self.edges = None # icosahedron edges
        # ---------------------------------------------------------------------------------------------
        # Coordinates systems.
        # '222'         2-fold symmetry along x, y, and z axes.
        # '222r'        '222' with 90 degree rotation around z.
        # '2n5'         2-fold symmetry along x and 5-fold along z.
        # '2n5r'        '2n5' with 180 degree rotation about y.
        # 'n25'         2-fold symmetry along y and 5-fold along z.
        # 'n25r'        'n25' with 180 degree rotation about x.
        # '2n3'         2-fold symmetry along x and 3-fold along z.
        # '2n3r'        '2n3' with 180 degree rotation about y.
        #
        self.coordinate_system_names = ('222', '222r', '2n5',
                                        '2n5r', 'n25', 'n25r', '2n3', '2n3r')
        self.icosahedronEdgeLength = \
            self.icosahedronEdgeLength(circumscribed_radius)
        self.angle23, self.angle25, self.angle35 = self.icosahedronAngles()
        self.icosahedronGeometry()

    # ---------------------------------------------------------------------------------------------
    # 60 icosahedral symmetry matrices.
    #
    def icosahedralSymmetryMatrices(self):
        t = self.icosahedralMatrixTable()
        tflist = _recenterSymmetries(t.get(self.orientation, None),
                                     self.center)
        return tflist


    # -------------------------------------------------------------------------
    # Edge length of icosahedron with a certain radio of the circumscribed
    # sphere:
    # According to Radio of the circumscribed sphere = 1 (vertices),
    #  the icosahedronEdgeLength should be 1.0515

    def icosahedronEdgeLength(self, circumscribed_radius):
        return 4 * circumscribed_radius / sqrt(10 + 2 * sqrt(5))

    # -------------------------------------------------------------------------
    # Matrices for mapping between different icosahedron coordinate frames.
    #
    def coordinateSystemTransform(self, from_cs, to_cs):

        self.cst = {}
        if self.cst:
            return self.cst[(from_cs, to_cs)]

        transform = self.cst

        s25 = self.e / 2  # Sin/Cos for angle between 2-fold and 5-fold axis
        c25 = sqrt(1 - s25 * s25)
        # Sin/Cos for angle between 3-fold and 5-fold axis
        s35 = self.e / sqrt(3)
        c35 = sqrt(1 - s35 * s35)

        transform[('2n5', '222')] = ((1, 0, 0, 0),
                                     (0, c25, -s25, 0),
                                     (0, s25, c25, 0))
        transform[('2n5', '2n3')] = ((1, 0, 0, 0),
                                     (0, c35, s35, 0),
                                     (0, -s35, c35, 0))

        # Axes permutations.
        # 90 degree rotation about z
        transform[('222', '222r')] = ((0, 1, 0, 0),
                                      (-1, 0, 0, 0),
                                      (0, 0, 1, 0))
        # 180 degree rotation about y
        transform[('2n3', '2n3r')] = \
            transform[('2n5', '2n5r')] = ((-1, 0, 0, 0),
                                          (0, 1, 0, 0),
                                          (0, 0, -1, 0))
        # 180 degree rotation about x
        transform[('n25', 'n25r')] = ((1, 0, 0, 0),
                                      (0, -1, 0, 0),
                                      (0, 0, -1, 0))
        # x <-> y and z -> -z
        transform[('n25', '2n5')] = ((0, 1, 0, 0),
                                     (1, 0, 0, 0),
                                     (0, 0, -1, 0))

        # Extend to all pairs of transforms.
        tlist = []
        while len(transform) > len(tlist):

            tlist = transform.keys()

            # Add inverse transforms
            for f, t in tlist:
                if not (t, f) in transform:
                    transform[(t, f)] = _transposeMatrix(transform[(f, t)])

            # Use transitivity
            for f1, t1 in tlist:
                for f2, t2 in tlist:
                    if f2 == t1 and f1 != t2 and not (f1, t2) in transform:
                        transform[(f1, t2)] = _multiplyMatrices(transform[(
                            f2, t2)], transform[(f1, t1)])

        i = _identityMatrix()
        for s in self.coordinate_system_names:
            transform[(s, s)] = i

        return transform[(from_cs, to_cs)]

    # ---------------------------------------------------------------------------------------
    # Compute icosahedral transformation matrices for different
    # coordinate systems.
    #

    def icosahedralMatrixTable(self):
        global icos_matrices
        if icos_matrices:
            return icos_matrices

        c = cos(2 * pi / 5)  # .309016994
        c2 = cos(4 * pi / 5)  # -.809016994

        icos_matrices['222'] = (

            ((1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0)),

            ((c2, -0.5, c, 0.0),
             (-0.5, c, c2, 0.0),
             (c, c2, -0.5, 0.0)),

            ((0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0)),

            ((-c2, -0.5, -c, 0.0),
             (-0.5, -c, c2, 0.0),
             (c, -c2, -0.5, 0.0)),

            ((0.5, c, c2, 0.0),
             (-c, c2, -0.5, 0.0),
             (c2, 0.5, -c, 0.0)),

            ((-c, c2, -0.5, 0.0),
             (c2, 0.5, -c, 0.0),
             (0.5, c, c2, 0.0)),

            ((c2, 0.5, -c, 0.0),
             (0.5, c, c2, 0.0),
             (-c, c2, -0.5, 0.0)),

            ((c2, -0.5, -c, 0.0),
             (0.5, -c, c2, 0.0),
             (c, c2, 0.5, 0.0)),

            ((-c, -c2, -0.5, 0.0),
             (c2, -0.5, -c, 0.0),
             (-0.5, c, -c2, 0.0)),

            ((0.5, -c, c2, 0.0),
             (-c, -c2, -0.5, 0.0),
             (-c2, 0.5, c, 0.0)),

            ((0.0, 0.0, -1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0)),

            ((-0.5, -c, c2, 0.0),
             (c, -c2, -0.5, 0.0),
             (-c2, -0.5, -c, 0.0)),

            ((-0.5, c, c2, 0.0),
             (c, c2, -0.5, 0.0),
             (c2, -0.5, c, 0.0)),

            ((-c, c2, -0.5, 0.0),
             (-c2, -0.5, c, 0.0),
             (-0.5, -c, -c2, 0.0)),

            ((c2, 0.5, -c, 0.0),
             (-0.5, -c, -c2, 0.0),
             (c, -c2, 0.5, 0.0)),

            ((0.5, c, c2, 0.0),
             (c, -c2, 0.5, 0.0),
             (-c2, -0.5, c, 0.0)),

            ((-0.5, c, c2, 0.0),
             (-c, -c2, 0.5, 0.0),
             (-c2, 0.5, -c, 0.0)),

            ((0.0, 0.0, -1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0)),

            ((-0.5, -c, c2, 0.0),
             (-c, c2, 0.5, 0.0),
             (c2, 0.5, c, 0.0)),

            ((0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0)),

            ((c2, 0.5, c, 0.0),
             (0.5, c, -c2, 0.0),
             (c, -c2, -0.5, 0.0)),

            ((-c2, 0.5, -c, 0.0),
             (0.5, -c, -c2, 0.0),
             (c, c2, -0.5, 0.0)),

            ((-c, -c2, -0.5, 0.0),
             (-c2, 0.5, c, 0.0),
             (0.5, -c, c2, 0.0)),

            ((0.5, -c, c2, 0.0),
             (c, c2, 0.5, 0.0),
             (c2, -0.5, -c, 0.0)),

            ((c2, -0.5, -c, 0.0),
             (-0.5, c, -c2, 0.0),
             (-c, -c2, -0.5, 0.0)),

            ((-c, c2, 0.5, 0.0),
             (c2, 0.5, c, 0.0),
             (-0.5, -c, c2, 0.0)),

            ((-c, -c2, 0.5, 0.0),
             (-c2, 0.5, -c, 0.0),
             (-0.5, c, c2, 0.0)),

            ((1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0)),

            ((c, -c2, -0.5, 0.0),
             (-c2, -0.5, -c, 0.0),
             (-0.5, -c, c2, 0.0)),

            ((c, c2, -0.5, 0.0),
             (c2, -0.5, c, 0.0),
             (-0.5, c, c2, 0.0)),

            ((-1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0)),

            ((-c2, 0.5, -c, 0.0),
             (-0.5, c, c2, 0.0),
             (-c, -c2, 0.5, 0.0)),

            ((0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0)),

            ((c2, 0.5, c, 0.0),
             (-0.5, -c, c2, 0.0),
             (-c, c2, 0.5, 0.0)),

            ((-0.5, -c, -c2, 0.0),
             (-c, c2, -0.5, 0.0),
             (-c2, -0.5, c, 0.0)),

            ((c, -c2, 0.5, 0.0),
             (c2, 0.5, -c, 0.0),
             (-0.5, -c, -c2, 0.0)),

            ((-c2, -0.5, c, 0.0),
             (0.5, c, c2, 0.0),
             (c, -c2, 0.5, 0.0)),

            ((-c2, 0.5, c, 0.0),
             (0.5, -c, c2, 0.0),
             (-c, -c2, -0.5, 0.0)),

            ((c, c2, 0.5, 0.0),
             (c2, -0.5, -c, 0.0),
             (0.5, -c, c2, 0.0)),

            ((-0.5, c, -c2, 0.0),
             (-c, -c2, -0.5, 0.0),
             (c2, -0.5, -c, 0.0)),

            ((0.0, 0.0, 1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0)),

            ((0.5, c, -c2, 0.0),
             (c, -c2, -0.5, 0.0),
             (c2, 0.5, c, 0.0)),

            ((0.5, -c, -c2, 0.0),
             (c, c2, -0.5, 0.0),
             (-c2, 0.5, -c, 0.0)),

            ((c, -c2, 0.5, 0.0),
             (-c2, -0.5, c, 0.0),
             (0.5, c, c2, 0.0)),

            ((-c2, -0.5, c, 0.0),
             (-0.5, -c, -c2, 0.0),
             (-c, c2, -0.5, 0.0)),

            ((-0.5, -c, -c2, 0.0),
             (c, -c2, 0.5, 0.0),
             (c2, 0.5, -c, 0.0)),

            ((0.5, -c, -c2, 0.0),
             (-c, -c2, 0.5, 0.0),
             (c2, -0.5, c, 0.0)),

            ((0.0, 0.0, 1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0)),

            ((0.5, c, -c2, 0.0),
             (-c, c2, 0.5, 0.0),
             (-c2, -0.5, -c, 0.0)),

            ((0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0)),

            ((-c2, -0.5, -c, 0.0),
             (0.5, c, -c2, 0.0),
             (-c, c2, 0.5, 0.0)),

            ((c2, -0.5, c, 0.0),
             (0.5, -c, -c2, 0.0),
             (-c, -c2, 0.5, 0.0)),

            ((c, c2, 0.5, 0.0),
             (-c2, 0.5, c, 0.0),
             (-0.5, c, -c2, 0.0)),

            ((-0.5, c, -c2, 0.0),
             (c, c2, 0.5, 0.0),
             (-c2, 0.5, c, 0.0)),

            ((-c2, 0.5, c, 0.0),
             (-0.5, c, -c2, 0.0),
             (c, c2, 0.5, 0.0)),

            ((c, -c2, -0.5, 0.0),
             (c2, 0.5, c, 0.0),
             (0.5, c, -c2, 0.0)),

            ((c, c2, -0.5, 0.0),
             (-c2, 0.5, -c, 0.0),
             (0.5, -c, -c2, 0.0)),

            ((-1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0)),

            ((-c, c2, 0.5, 0.0),
             (-c2, -0.5, -c, 0.0),
             (0.5, c, -c2, 0.0)),

            ((-c, -c2, 0.5, 0.0),
             (c2, -0.5, c, 0.0),
             (0.5, -c, -c2, 0.0)),

        )

        for cs in self.coordinate_system_names:
            if cs != '222':
                t = self.coordinateSystemTransform(cs, '222')
                tinv = self.coordinateSystemTransform('222', cs)
                icos_matrices[cs] = [_multiplyMatrices(tinv, m, t)
                                     for m in icos_matrices['222']]
        return icos_matrices

    # -----------------------------------------------------------------------------

    def icosahedronAngles(self):
        # Angles between 2-fold, 3-fold and 5-fold
        # symmetry axes of an icosahedron.
        angle25 = asin(self.e/2)
        angle35 = asin(self.e/sqrt(3))
        angle23 = asin(self.e/(2*sqrt(3)))
        return angle23, angle25, angle35

    def icosahedronGeometry(self):

        a = 2 * self.angle25  # Angle spanned by edge from center
        # 5-fold symmetry axis along z
        c5 = cos(2 * pi / 5)
        s5 = sin(2 * pi / 5)
        tf5 = ((c5, -s5, 0, 0),
               (s5, c5, 0, 0),
               (0, 0, 1, 0))

        # 2-fold symmetry axis along x
        tf2 = ((1, 0, 0, 0),
               (0, -1, 0, 0),
               (0, 0, -1, 0))

        p = tuple(map(operator.add, self.center,
                      (0, 0, self.circumscribed_radius)))
        p50 = tuple(map(operator.add, self.center,
                        (0, self.circumscribed_radius * sin(a),
                         self.circumscribed_radius * cos(a))))
        p51 = _applyMatrix(tf5, p50)
        p52 = _applyMatrix(tf5, p51)
        p53 = _applyMatrix(tf5, p52)
        p54 = _applyMatrix(tf5, p53)
        self.vertices = [p, p50, p51, p52, p53, p54]
        self.vertices.extend([_applyMatrix(tf2, q) for q in self.vertices])
        if self.orientation != '2n5':
            self.tf = self.coordinateSystemTransform('2n5', self.orientation)
            self.vertices = [_applyMatrix(self.tf, p) for p in self.vertices]

            #
            # Vertex numbering
            #
            #  Top   1          Bottom
            #      2   5        9   10
            #        0            6
            #      3   4        8   11
            #                     7
        # 20 triangles composing icosahedron.
        #
        self.triangles = ((0,1,2), (0,2,3), (0,3,4), (0,4,5), (0,5,1),
                         (6,7,8), (6,8,9), (6,9,10), (6,10,11), (6,11,7),
                         (1,9,2), (2,9,8), (2,8,3), (3,8,7), (3,7,4),
                         (4,7,11), (4,11,5), (5,11,10), (5,10,1), (1,10,9))

        for triangle in self.triangles:
            self._3foldAxis.append( [(item1 + item2 + item3) /3. \
                for item1,  item2,  item3 in zip(self.vertices[triangle[0]],
                                                 self.vertices[triangle[1]],
                                                 self.vertices[triangle[2]])])

        self.edges = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
                      (1, 2), (2, 3), (3, 4), (4, 5), (5, 1),
                      (6, 7), (6, 8), (6, 9), (6, 10), (6, 11),
                      (1, 9), (2, 8), (2, 9), (3, 7), (3, 8),
                      (7, 8), (8, 9), (9, 10), (10, 11), (11, 7),
                      (4, 7), (4, 11), (5, 10), (5, 11), (1, 10))

        for edge in self.edges:
            self._2foldAxis.append([(item1 + item2 ) / 2. \
                for item1, item2 in zip(self.vertices[edge[0]],
                                        self.vertices[edge[1]])])

    def getVertices(self):
        return self.vertices

    def getTriangles(self):
        return self.triangles

    def getEdges(self):
        return self.edges

    def get3foldAxis(self):
        return self._3foldAxis

    def get2foldAxis(self):
        return self._2foldAxis
