# ***************************************************************************
# * Authors:    Roberto Marabini (roberto@cnb.csic.es)
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
import os
from tempfile import mkstemp
from pyworkflow.utils import runJob
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.packages.xmipp3 import getEnviron
from pyworkflow.em.protocol import ProtImportVolumes
from pyworkflow.em.packages.xmipp3.protocol_extract_unit_cell import XmippProtExtractUnit
from pyworkflow.em.constants import SYM_I222r
from pyworkflow.em.convert import ImageHandler

class TestProtModelBuilding(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        geometricPhantom  = """# Phantom description file, (generated with phantom help)
# General Volume Parameters: 
#      Xdim      Ydim      Zdim   Background_Density Scale 
       5 5 5 0 60  
# Feature Parameters: 
#Type  +/=  Density X_Center Y_Center Z_Center 
sph  + 1.   1.000 0.000 1.618 .15
sph  + 1.   1.000 0.000 -1.618 .15
sph  + 1.   -1.000 0.000 1.618 .15
sph  + 1.   -1.000 0.000 -1.618 .15
sph  + 1.   0.000 1.618 1.000 .15
sph  + 1.   0.000 -1.618 1.000 .15
sph  + 1.   0.000 1.618 -1.000 .15
sph  + 1.   0.000 -1.618 -1.000 .15
sph  + 1.   1.618 1.000 0.000 .15
sph  + 1.   -1.618 1.000 0.000 .15
sph  + 1.   1.618 -1.000 0.000 .15
sph  + 1.   -1.618 -1.000 0.000 .15
sph  + 1.   0.000 -0.539 1.412 .10
sph  + 1.   0.000 0.539 1.412 .10
sph  + 1.   0.873 -0.873 0.873 .10
sph  + 1.   0.873 0.873 0.873 .10
sph  + 1.   1.412 0.000 0.539 .10
sph  + 1.   0.000 0.539 -1.412 .10
sph  + 1.   0.873 0.873 -0.873 .10
sph  + 1.   0.000 -0.539 -1.412 .10
sph  + 1.   1.412 0.000 -0.539 .10
sph  + 1.   0.873 -0.873 -0.873 .10
sph  + 1.   -0.873 0.873 0.873 .10
sph  + 1.   -1.412 0.000 0.539 .10
sph  + 1.   -0.873 -0.873 0.873 .10
sph  + 1.   -0.873 0.873 -0.873 .10
sph  + 1.   -1.412 0.000 -0.539 .10
sph  + 1.   -0.873 -0.873 -0.873 .10
sph  + 1.   -0.539 1.412 0.000 .10
sph  + 1.   0.539 1.412 0.000 .10
sph  + 1.   -0.539 -1.412 0.000 .10
sph  + 1.   0.539 -1.412 0.000 .10
sph  + 1.   0.000 0.000 1.618 .10
sph  + 1.   -0.500 -0.809 1.309 .10
sph  + 1.   0.500 -0.809 1.309 .10
sph  + 1.   0.500 0.809 1.309 .10
sph  + 1.   -0.500 0.809 1.309 .10
sph  + 1.   0.000 0.000 1.618 .10
sph  + 1.   0.500 -0.809 1.309 .10
sph  + 1.   0.809 -1.309 0.500 .10
sph  + 1.   1.309 -0.500 0.809 .10
sph  + 1.   1.309 0.500 0.809 .10
sph  + 1.   0.809 1.309 0.500 .10
sph  + 1.   0.500 0.809 1.309 .10
sph  + 1.   1.309 -0.500 0.809 .10
sph  + 1.   1.618 0.000 0.000 .10
sph  + 1.   1.309 0.500 0.809 .10
sph  + 1.   0.000 0.000 -1.618 .10
sph  + 1.   -0.500 0.809 -1.309 .10
sph  + 1.   0.500 0.809 -1.309 .10
sph  + 1.   0.500 0.809 -1.309 .10
sph  + 1.   0.809 1.309 -0.500 .10
sph  + 1.   1.309 0.500 -0.809 .10
sph  + 1.   0.500 -0.809 -1.309 .10
sph  + 1.   -0.500 -0.809 -1.309 .10
sph  + 1.   0.000 0.000 -1.618 .10
sph  + 1.   1.309 0.500 -0.809 .10
sph  + 1.   1.618 0.000 0.000 .10
sph  + 1.   1.309 -0.500 -0.809 .10
sph  + 1.   1.309 -0.500 -0.809 .10
sph  + 1.   0.809 -1.309 -0.500 .10
sph  + 1.   0.500 -0.809 -1.309 .10
sph  + 1.   -0.500 0.809 1.309 .10
sph  + 1.   -0.809 1.309 0.500 .10
sph  + 1.   -1.309 0.500 0.809 .10
sph  + 1.   -1.309 0.500 0.809 .10
sph  + 1.   -1.618 0.000 0.000 .10
sph  + 1.   -1.309 -0.500 0.809 .10
sph  + 1.   -1.309 -0.500 0.809 .10
sph  + 1.   -0.809 -1.309 0.500 .10
sph  + 1.   -0.500 -0.809 1.309 .10
sph  + 1.   -0.500 0.809 -1.309 .10
sph  + 1.   -0.809 1.309 -0.500 .10
sph  + 1.   -1.309 0.500 -0.809 .10
sph  + 1.   -1.309 0.500 -0.809 .10
sph  + 1.   -1.618 0.000 0.000 .10
sph  + 1.   -1.309 -0.500 -0.809 .10
sph  + 1.   -1.309 -0.500 -0.809 .10
sph  + 1.   -0.809 -1.309 -0.500 .10
sph  + 1.   -0.500 -0.809 -1.309 .10
sph  + 1.   0.000 1.618 0.000 .10
sph  + 1.   -0.809 1.309 -0.500 .10
sph  + 1.   -0.809 1.309 0.500 .10
sph  + 1.   0.809 1.309 0.500 .10
sph  + 1.   0.809 1.309 -0.500 .10
sph  + 1.   0.000 1.618 0.000 .10
sph  + 1.   0.000 -1.618 0.000 .10
sph  + 1.   -0.809 -1.309 -0.500 .10
sph  + 1.   -0.809 -1.309 0.500 .10
sph  + 1.   0.809 -1.309 0.500 .10
sph  + 1.   0.809 -1.309 -0.500 .10
sph  + 1.   0.000 -1.618 0.000 .10"""
        (fd, cls.filename) = mkstemp(suffix=".feat")
        f = os.fdopen(fd, "w")
        f.write(geometricPhantom)
        f.close()

    def test_extractunitCell(self):
        """ extract unit cell from icosahedral pahntom
            using xmipp_i2 symmetry
        """
        # create phantom (3D map)
        _, outputFile = mkstemp(suffix=".vol")
        command = "xmipp_phantom_create "
        args    = " -i %s"%self.filename
        args += " -o %s"%outputFile
        runJob(None, command, args,env=getEnviron())

        #import volume
        args = {'filesPath': outputFile,
                'filesPattern': '',
                'samplingRate': 1
                }
        prot = self.newProtocol(ProtImportVolumes, **args)
        prot.setObjLabel('import volume')
        self.launchProtocol(prot)

        # execute protocol extract unitCell
        args = {'inputVolumes': prot.outputVolume,
                'symmetryGroup': SYM_I222r,
                'innerRadius': 71.,
                'outerRadius': 141.,
                'expandFactor': .2,
                'offset': 0.
                }
        prot = self.newProtocol(XmippProtExtractUnit, **args)
        prot.setObjLabel('extract unit cell')
        self.launchProtocol(prot)

        #check results
        print "qqqq",prot.outputVolume.getFileName()
        ih = ImageHandler()
        xdim, ydim, zdim, ndim = ih.getDimensions(prot.outputVolume.getFileName())
        self.assertEqual(xdim,125)
        self.assertEqual(ydim,121)
        self.assertEqual(zdim,91)
        #delete temporary file
        os.remove(self.filename)
        os.remove(outputFile)
