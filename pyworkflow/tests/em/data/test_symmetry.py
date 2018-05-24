# **************************************************************************
# *
# * Authors:     roberto Marabini (roberto@cnb.csic.es)
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
#
from pyworkflow.tests import *
from pyworkflow.em.constants import SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL, \
    SYM_TETRAHEDRAL, SYM_I222, SYM_I222r, SYM_In25, SYM_In25r
from pyworkflow.em.symmetry import getSymmetryMatrices
from pyworkflow.em.transformations import identity_matrix
import numpy


class TestSymmetry(unittest.TestCase):

    def assertArrayAlmostEqual(self, a1, a2):
        try:
            numpy.testing.assert_array_almost_equal(a1, a2, decimal=3)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

    def testSymmetryCyclic(self):
        matrices = getSymmetryMatrices(SYM_CYCLIC, n=3)
        m1 = matrices[0]
        m2 = matrices[1]
        m3 = matrices[2]

        #m1 is identity
        r = identity_matrix()
        self.assertArrayAlmostEqual(r, m1)

        #m2 rotated 120 degrees
        c = -0.5
        s = 0.86602540378
        r = [[c, -s, 0, 0],
             [s, c, 0, 0],
             [0, 0, 1.0, 0],
             [0, 0, 0.0, 1.0]]

        self.assertArrayAlmostEqual(r, m2)
        #m3 rotated -120 degrees
        c = -0.5
        s = -0.86602540378
        r = [[c, -s, 0, 0],
             [s, c, 0, 0],
             [0, 0, 1.0, 0],
             [0, 0, 0.0, 1.0]]

        self.assertArrayAlmostEqual(r, m3)

    def testSymmetryDihedral(self):
        matrices = getSymmetryMatrices(SYM_DIHEDRAL, n=3)
        # skip m1, m2 and m3 since there are identical to CYCLIC
        m4 = matrices[3]
        m5 = matrices[4]
        m6 = matrices[5]

        #m4 mirrored
        r = identity_matrix()
        r[1][1] *= -1.
        r[2][2] *= -1.
        self.assertArrayAlmostEqual(r, m4)

        #m5 rotated 120 degrees and mirrored
        c = -0.5
        s = 0.86602540378
        r = [[c, -s, 0, 0],
             [-s, -c, 0, 0],
             [0, 0, -1.0, 0],
             [0, 0, 0.0, 1.0]]
        self.assertArrayAlmostEqual(r, m5)

        # m6 rotated -120 degrees and mirrored
        c = -0.5
        s = -0.86602540378
        r = [[c, -s, 0, 0],
             [-s, -c, 0, 0],
             [0, 0, -1.0, 0],
             [0, 0, 0.0, 1.0]]
        self.assertArrayAlmostEqual(r, m6)

    @classmethod
    def tearDownClass(cls):
        pass

    @classmethod
    def setUpClass(cls):
        pass