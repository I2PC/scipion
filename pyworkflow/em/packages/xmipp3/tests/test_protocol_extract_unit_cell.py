# ***************************************************************************
# * Authors:    Roberto Marabini (roberto@cnb.csic.es)
# *             Marta Martinez (mmmtnez@cnb.csic.es)
# *
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
# ***************************************************************************/

# general protocol to test the appropriate extraction of the unit cell from
# symmetries cyclic,dihedral, tetrahedral, octahedral, and icosahedral,
# that are available to select in the form for the extraction of the unit
# cell from a specific volume

import os
import math
from tempfile import mkstemp

from pyworkflow.em.symmetry import Icosahedron
from pyworkflow.em.constants import SYM_I222r, SYM_I222, SCIPION_SYM_NAME, \
    SYM_In25, SYM_In25r, SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL, \
    SYM_OCTAHEDRAL
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Transform
from pyworkflow.em.convert_header.CCP4.convert import Ccp4Header
from pyworkflow.em.packages.xmipp3 import getEnviron
from pyworkflow.em.packages.xmipp3.constants import XMIPP_SYM_NAME
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms \
    import XmippProtConvertToPseudoAtoms
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base \
    import NMA_MASK_THRE
from pyworkflow.em.packages.xmipp3.protocol_extract_unit_cell \
    import XmippProtExtractUnit
from pyworkflow.em.protocol import ProtImportVolumes
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import runJob

OFFSET = 22.5

# function to write the coordinates of a phantom (3D map) for icosahedral
# symmetry in a text file


def generate_ico(sym, mode, f):
    icosahedron = Icosahedron(orientation=sym)
    center = []
    x = 0.
    y = 0.
    z = 0.
    if mode == 'xmipp':
        f.write("sph  + 1. %.3f %.3f %.3f .65\n" % (x, y, z))
    else:
        f.write('.sphere %.3f %.3f %.3f .65\n' % (x, y, z))
    # print 5fold points
    for vertice in icosahedron.getVertices():
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n" %
                    (vertice[0], vertice[1], vertice[2]))
        else:
            f.write('.sphere %.3f %.3f %.3f .15\n' %
                    (vertice[0], vertice[1], vertice[2]))

    # print 3fold points
    if mode == 'xmipp':
        pass
    else:
        f.write(""".comment 3fold
    .color yellow
    """)

    for _3fold in icosahedron.get3foldAxis():
        x, y, z = _3fold
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .10\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .10\n' % (x, y, z))

    # print 2fold points
    if mode == 'xmipp':
        pass
    else:
        f.write(""".comment 2fold
    .color green
    """)
    for _2fold in icosahedron.get2foldAxis():
        x, y, z = _2fold

        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .09\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .09\n' % (x, y, z))

# function to write the coordinates of a phantom (3D map) for cyclic symmetry
# (order 8) in a text file


def generate_cyclic(order, offset, mode, f, ):
    center = []
    z_value = [-.45, 0]
    if mode == 'xmipp':
        f.write("cyl  +  1. 0 0 +.45 1.0 1.0 .2 0 0 0\n")
        f.write("cyl  + -1. 0 0 +.45 0.8 0.8 .2 0 0 0\n")
        f.write("cyl  +  1. 0 0 +.45 0.5 0.5 .2 0 0 0\n")
        f.write("cyl  + -1. 0 0 +.45 0.5 0.5 .2 0 0 0\n")
    else:
        f.write('.cylinder 0 0 +.46 0 0 +.44 1 open\n')
        f.write('.cylinder 0 0 +.46 0 0 +.44 .5 open\n')
    for z in z_value:
        x = 0.
        y = 0.
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .15\n' % (x, y, z))
        for point in range(order):
            x = math.cos(2*point*math.pi/order + offset)
            y = math.sin(2*point*math.pi/order + offset)
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x, y, z))
            else:
                f.write('.sphere %.3f %.3f %.3f .15\n' % (x, y, z))
        for point in range(order):
            x = math.cos(2*point*math.pi/order + offset)/2
            y = math.sin(2*point*math.pi/order + offset)/2
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .10\n" % (x, y, z))
            else:
                f.write('.sphere %.3f %.3f %.3f .10\n' % (x, y, z))
        for point in range(order):
            x = math.cos(2*point*math.pi/order + offset)/4
            y = math.sin(2*point*math.pi/order + offset)/4
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .05\n" % (x, y, z))
            else:
                f.write('.sphere %.3f %.3f %.3f .05\n' % (x, y, z))

# function to write the coordinates of a phantom (3D map) for dihedral symmetry
# (order 8) in a text file


def generate_dihedral(order, offset, mode, f, ):
    center = []
    z_value = [-.45, 0., .45]
    for z in z_value:
        x = 0.
        y = 0.
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .15\n' % (x, y, z))
        for point in range(order):
            x = math.cos(2*point*math.pi/order + offset)
            y = math.sin(2*point*math.pi/order + offset)
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x, y, z))
            else:
                f.write('.sphere %.3f %.3f %.3f .15\n' % (x, y, z))
        for point in range(order):
            x = math.cos(2*point*math.pi/order + offset)/2
            y = math.sin(2*point*math.pi/order + offset)/2
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .10\n" % (x, y, z))
            else:
                f.write('.sphere %.3f %.3f %.3f .10\n' % (x, y, z))
        for point in range(order):
            x = math.cos(2*point*math.pi/order + offset)/4
            y = math.sin(2*point*math.pi/order + offset)/4
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .05\n" % (x, y, z))
            else:
                f.write('.sphere %.3f %.3f %.3f .05\n' % (x, y, z))

# function to write the coordinates of a phantom (3D map) for tetrahedral
# symmetry in a text file


def generate_tetrahedral(mode, f, ):
    centroid = [0., 0., 0.]
    rmax = 1.
    _3f = [0., 0., 1.]
    _3fp = [0., 0.94281, - 0.33333]
    _3fpp = [0.81650, - 0.47140, - 0.33333]
    _3fppp = [-0.81650, - 0.47140, - 0.33333]
    points = [centroid, _3f, _3fp, _3fpp, _3fppp]
    vertices = []
    for point in points:
        x = point[0] * rmax
        y = point[1] * rmax
        z = point[2] * rmax
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .15\n' % (x, y, z))
        if (x, y, z) != (0., 0., 0.):
            vertices.append([x, y, z])
    edges = []
    for i in range(len(vertices)):
        for j in range(len(vertices)):
            if i != j:
                x = (0.5) * (vertices[i][0] + vertices[j][0])
                y = (0.5) * (vertices[i][1] + vertices[j][1])
                z = (0.5) * (vertices[i][2] + vertices[j][2])
                edge = [x, y, z]
                if (x, y, z) != (0., 0., 0.):
                    if edge not in edges:
                        edges.append(edge)
                        if mode == 'xmipp':
                            f.write("sph  + 1. %.3f %.3f %.3f .05\n" %
                                    (edge[0], edge[1], edge[2]))
                        else:
                            f.write('.sphere %.3f %.3f %.3f .05\n' %
                                    (edge[0], edge[1], edge[2]))

# function to write the coordinates of a phantom (3D map) for octahedral
# symmetry in a text file


def generate_octahedral(mode, f, ):
    centre = [0., 0., 0.]
    rmax = 1.
    _4f = [1., 0., 0.]
    _4fp = [0., 1., 0.]
    _4fpp = [0., 0., 1.]
    _4_2f = [-1., 0., 0.]
    _4_2fp = [0., -1., 0.]
    _4_2fpp = [0., 0., -1.]
    points = [centre, _4f, _4fp, _4fpp, _4_2f, _4_2fp, _4_2fpp]
    vertices = []
    for point in points:
        x = point[0] * rmax
        y = point[1] * rmax
        z = point[2] * rmax
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .15\n' % (x, y, z))
        if (x, y, z) != (0., 0., 0.):
            vertices.append([x, y, z])
    edges = []
    for i in range(len(vertices)):
        for j in range(len(vertices)):
            if i != j:
                x = (0.5) * (vertices[i][0] + vertices[j][0])
                y = (0.5) * (vertices[i][1] + vertices[j][1])
                z = (0.5) * (vertices[i][2] + vertices[j][2])
                edge = [x, y, z]
                if (x, y, z) != (0., 0., 0.):
                    if edge not in edges:
                        edges.append(edge)
                        if mode == 'xmipp':
                            f.write("sph  + 1. %.3f %.3f %.3f .10\n" %
                                    (edge[0], edge[1], edge[2]))
                        else:
                            f.write('.sphere %.3f %.3f %.3f .10\n' %
                                    (edge[0], edge[1], edge[2]))
    baricentres = [[(1./3.)*rmax, (1./3.)*rmax, (1./3.)*rmax],
                   [(1./3.)*rmax, (1./3.)*rmax, (-1./3.)*rmax],
                   [(1./3.)*rmax, (-1./3.)*rmax, (1./3.)*rmax],
                   [(1./3.)*rmax, (-1./3.)*rmax, (-1./3.)*rmax],
                   [(-1./3.)*rmax, (1./3.)*rmax, (1./3.)*rmax],
                   [(-1./3.)*rmax, (1./3.)*rmax, (-1./3.)*rmax],
                   [(-1./3.)*rmax, (-1./3.)*rmax, (1./3.)*rmax],
                   [(-1./3.)*rmax, (-1./3.)*rmax, (-1./3.)*rmax]]
    for baricentre in baricentres:
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .05\n" %
                    (baricentre[0], baricentre[1], baricentre[2]))
        else:
            f.write('.sphere %.3f %.3f %.3f .05\n' %
                    (baricentre[0], baricentre[1], baricentre[2]))

# general function to generate the coordinates files for phantoms (3D map)
# for every symmetry


def generate(sym='I2n5', mode='xmipp', suffix="_i2", offset=0.):
    offset = math.radians(offset)
    symPreffix = sym[:1]
    symSuffix = sym[1:]
    if mode == 'xmipp':
        suffix += ".feat"
    else:
        suffix += ".bild"

    (fd, filename) = mkstemp(suffix=suffix)
    f = os.fdopen(fd, "w")

    if mode == 'xmipp':
        f.write("""# Phantom description file, (generated with phantom help)
# General Volume Parameters:
#      Xdim      Ydim      Zdim   Background_Density Scale
       3 3 3 0 60
# Feature Parameters:
#Type  +/=  Density X_Center Y_Center Z_Center
""")
    else:
        f.write(""".scale 60
.comment system of coordinates
.color 1 0 0
.arrow 0 0 0 2 0 0 0.025
.color 0 1 0
.arrow 0 0 0 0 2 0 0.025
.color 0 0 1
.arrow 0 0 0 0 0 2 0.025
.comment vertices
.color red
""")
    if symPreffix == 'I':
        generate_ico(symSuffix, mode, f)
    elif symPreffix == 'C':
        generate_cyclic(int(symSuffix), offset,  mode, f)
    elif symPreffix == 'D':
        generate_dihedral(int(symSuffix), offset,  mode, f)
    elif symPreffix == 'T':
        generate_tetrahedral(mode, f)
    elif symPreffix == 'O':
        generate_octahedral(mode, f)

    f.close()
    return filename

# general class to test the extraction of the unit cell of cyclic, dihedral,
# tetrahedral, octahedral, and isosahedral symmetries


class TestProtModelBuilding(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.mode = 'xmipp'
        # Cyclic
        cls.symOrder = 8  # phantom symmetry order for Cn and similars

        # general
        cls.innerRadius = None
        cls.outerRadius = None
        setupTestProject(cls)
        cls.filename = {}
        cls.box = {}

    # function to extract the unit cell of cyclic symmetry
    def test_extractunitCellcyclic(self):
        self.innerRadius = 0.
        self.outerRadius = 40.

        # C8
        self.filename[SYM_CYCLIC] = generate(
            SCIPION_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder),
            'xmipp', XMIPP_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder))
        self.box[SYM_CYCLIC] = (50, 45, 81)
        self.extractunitCell(SYM_CYCLIC)

        # C8 + offset
        self.filename[SYM_CYCLIC] = generate(
            SCIPION_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder),
            'xmipp', XMIPP_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder), OFFSET)
        self.box[SYM_CYCLIC] = (46, 48, 81)
        self.extractunitCell(SYM_CYCLIC, OFFSET)

        # C1
        self.filename[SYM_CYCLIC] = generate(
            SCIPION_SYM_NAME[SYM_CYCLIC][:1] + str(self.symOrder),
            'xmipp', XMIPP_SYM_NAME[SYM_CYCLIC][:1] + str(self.symOrder))
        self.box[SYM_CYCLIC] = (81, 81, 81)
        self.symOrder = 1
        self.extractunitCell(SYM_CYCLIC)

        # C2
        self.symOrder = 8  # set to 8 so a pretty 8 fold phantom is created
        self.filename[SYM_CYCLIC] = generate(
            SCIPION_SYM_NAME[SYM_CYCLIC][:1] + str(self.symOrder),
            'xmipp', XMIPP_SYM_NAME[SYM_CYCLIC][:1] + str(self.symOrder))
        self.box[SYM_CYCLIC] = (81, 81, 81)
        self.symOrder = 2
        self.extractunitCell(SYM_CYCLIC)

    # function to extract the unit cell of dihedral symmetry
    def test_extractunitCelldihedral(self):
        self.innerRadius = 0.
        self.outerRadius = 40.

        # D8
        self.filename[SYM_DIHEDRAL] = generate(
            SCIPION_SYM_NAME[SYM_DIHEDRAL][:1] + str(self.symOrder),
            'xmipp', XMIPP_SYM_NAME[SYM_DIHEDRAL][:1] + str(self.symOrder))
        self.box[SYM_DIHEDRAL] = (50, 45, 81)
        self.extractunitCell(SYM_DIHEDRAL)

        # D8 + offset
        self.filename[SYM_DIHEDRAL] = generate(
            SCIPION_SYM_NAME[SYM_DIHEDRAL][:1] + str(self.symOrder), 'xmipp',
            XMIPP_SYM_NAME[SYM_DIHEDRAL][:1] + str(self.symOrder), OFFSET)
        self.box[SYM_DIHEDRAL] = (46, 48, 81)
        self.extractunitCell(SYM_DIHEDRAL, OFFSET)

    # function to extract the unit cell of tetrahedral symmetry
    def test_extractunitCelltetrahedral(self):
        self.innerRadius = 0.
        self.outerRadius = 90.
        self.filename[SYM_TETRAHEDRAL] = generate(
            SCIPION_SYM_NAME[SYM_TETRAHEDRAL],
            'xmipp', XMIPP_SYM_NAME[SYM_TETRAHEDRAL])
        self.box[SYM_TETRAHEDRAL] = (80, 140, 181)
        self.extractunitCell(SYM_TETRAHEDRAL)

    # function to extract the unit cell of octahedral symmetry
    def test_extractunitCelloctahedral(self):
        self.innerRadius = 0.
        self.outerRadius = 90.
        self.filename[SYM_OCTAHEDRAL] = generate(
            SCIPION_SYM_NAME[SYM_OCTAHEDRAL],
            'xmipp', XMIPP_SYM_NAME[SYM_OCTAHEDRAL])
        self.box[SYM_OCTAHEDRAL] = (110, 55, 181)
        self.extractunitCell(SYM_OCTAHEDRAL)

    # function to extract the unit cell of icosahedral symmetry
    def test_extractunitCellIco(self):
        self.innerRadius = 37.
        self.outerRadius = 79.
        self.filename[SYM_I222] = generate(SCIPION_SYM_NAME[SYM_I222],
                                           'xmipp', XMIPP_SYM_NAME[SYM_I222])
        self.box[SYM_I222] = (94, 92, 53)
        self.filename[SYM_I222r] = generate(SCIPION_SYM_NAME[SYM_I222r],
                                            'xmipp', XMIPP_SYM_NAME[SYM_I222r])
        self.box[SYM_I222r] = (92, 71, 53)
        self.filename[SYM_In25] = generate(SCIPION_SYM_NAME[SYM_In25],
                                           'xmipp', XMIPP_SYM_NAME[SYM_In25])
        self.box[SYM_In25] = (70, 94, 115)
        self.filename[SYM_In25r] = generate(SCIPION_SYM_NAME[SYM_In25r],
                                            'xmipp', XMIPP_SYM_NAME[SYM_In25r])
        self.box[SYM_In25r] = (70, 71, 55)

        self.extractunitCell(SYM_I222)  # no crowther 222
        self.extractunitCell(SYM_I222r)  # crowther 222
        self.extractunitCell(SYM_In25)
        self.extractunitCell(SYM_In25r)

    def test_extractunitCellHalfIco(self):
        self.innerRadius = 37.
        self.outerRadius = 79.
        self.filename[SYM_I222r] = generate(SCIPION_SYM_NAME[SYM_I222r],
                                            'xmipp', XMIPP_SYM_NAME[SYM_I222r])
        self.box[SYM_I222r] = (91, 70, 53)

        self.extractunitCell(SYM_I222r, cropZ=True)  # crowther 222

    # general function to extract the unit cell
    def extractunitCell(self, sym, offset=0, cropZ=False):

        """ extract unit cell from icosahedral phantom
            using xmipp_i2 symmetry
        """
        # create phantom (3D map)
        _samplingRate = 1.34
        _, outputFile1 = mkstemp(suffix=".mrc")
        command = "xmipp_phantom_create "
        args = " -i %s" % self.filename[sym]
        args += " -o %s" % outputFile1
        runJob(None, command, args, env=getEnviron())
        ccp4header = Ccp4Header(outputFile1, readHeader=True)
        x, y, z = ccp4header.getDims()
        t = Transform()

        if cropZ:
            _, outputFile2 = mkstemp(suffix=".mrc")
            args = "-i %s -o %s" % (outputFile1, outputFile2)
            args += " --corners "
            args += " %d " % (- x / 2.)
            args += " %d " % (- y / 2.)
            args += " %d " % (0.)
            args += " %d " % (+ x / 2.)
            args += " %d " % (+ y / 2.)
            args += " %d " % (+ z / 2.)
            runJob(None, "xmipp_transform_window", args, env=getEnviron())
            t.setShifts(0, 0, 0)
            outputFile = outputFile2
            ccp4header = Ccp4Header(outputFile2, readHeader=True)

        else:
            t.setShifts(0, 0, 0)
            outputFile = outputFile1

        ccp4header.setSampling(_samplingRate)
        ccp4header.setOrigin(t.getShifts())
        ccp4header.writeHeader()

        # import volume
        if cropZ:
            args = {'filesPath': outputFile,
                    'filesPattern': '',
                    'samplingRate': _samplingRate,
                    'copyFiles': True,
                    'setOrigCoord': True,
                    'x': 90. * _samplingRate,
                    'y': 90. * _samplingRate,
                    'z': 0.
                    # x, y, z in Angstroms
                    }
        else:
            args = {'filesPath': outputFile,
                    'filesPattern': '',
                    'samplingRate': _samplingRate,
                    'copyFiles': True,
                    'setDefaultOrigin': False,
                    }
        prot = self.newProtocol(ProtImportVolumes, **args)
        prot.setObjLabel('import volume(%s)' % XMIPP_SYM_NAME[sym])
        self.launchProtocol(prot)
        # execute protocol extract unitCell
        args = {'inputVolumes': prot.outputVolume,
                'symmetryGroup': sym,
                'symmetryOrder': self.symOrder,
                'innerRadius': self.innerRadius,
                'outerRadius': self.outerRadius,
                'expandFactor': .2,
                'offset': offset
                }
        prot = self.newProtocol(XmippProtExtractUnit, **args)
        prot.setObjLabel('extract unit cell')
        self.launchProtocol(prot)

        # check results
        ih = ImageHandler()
        xdim, ydim, zdim, ndim = \
            ih.getDimensions(prot.outputVolume.getFileName())
        self.assertTrue(abs(xdim - self.box[sym][0])<2)
        self.assertTrue(abs(ydim - self.box[sym][1])<2)
        self.assertTrue(abs(zdim - self.box[sym][2])<2)

        # create pdb fileoutput
        args = {'inputStructure': prot.outputVolume,
                'maskMode': NMA_MASK_THRE,
                'maskThreshold': 0.5,
                'pseudoAtomRadius': 1.5
                }
        prot = self.newProtocol(XmippProtConvertToPseudoAtoms, **args)
        prot.setObjLabel('get pdb')
        self.launchProtocol(prot)

        # check results
        filenamePdb = prot._getPath('pseudoatoms.pdb')
        self.assertTrue(os.path.isfile(filenamePdb))
        # delete temporary files
        os.remove(self.filename[sym])
        os.remove(outputFile)
