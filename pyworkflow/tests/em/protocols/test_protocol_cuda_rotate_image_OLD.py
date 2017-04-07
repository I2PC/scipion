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
SCALE=60
DIMENSION = 64#300
NUMIMAGES = 10#00
SKIP=False

class TestProtModelBuilding(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        geometricPhantom  = """# Phantom description file, (generated with phantom help)
# General Volume Parameters: 
#      Xdim      Ydim      Zdim   Background_Density Scale 
       5 5 5 0 %d
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
sph  + 1.   0.000 -1.618 0.000 .10"""%SCALE
        (fd, cls.filename) = mkstemp(suffix=".feat")
        f = os.fdopen(fd, "w")
        f.write(geometricPhantom)
        f.close()

        geometricParams = """# XMIPP_STAR_1 *
# Projection Parameters
data_block1
# X and Y projection dimensions [Xdim Ydim]
_dimensions2D   '%d %d'
# Rotation range and number of samples [Start Finish NSamples]
_projRotRange    '0 0.00001 1'
# Rotation angle added noise  [noise (bias)]
_projRotNoise   '0'
# Tilt range and number of samples for Tilt
_projTiltRange    '0 0.00001 1'
# Tilt angle added noise
_projTiltNoise   '0'
# Psi range and number of samples
_projPsiRange    '0 360 %d'
# Psi added noise
_projPsiNoise   '0'
# Noise applied to pixels [noise (bias)]
_noisePixelLevel   '0'
# Noise applied to particle center coordenates [noise (bias)]
_noiseCoord   '0'
"""%(DIMENSION,DIMENSION,NUMIMAGES)
        (fd, cls.filenameImages) = mkstemp(suffix=".xmd")
        f = os.fdopen(fd, "w")
        f.write(geometricParams)
        f.close()

    def test_cudaRotateImage(self):
        """ extract unit cell from icosahedral pahntom
            using xmipp_i2 symmetry
        """
        # create phantom (3D map)
        _, volFile = mkstemp(suffix=".vol")
        _, projFile = mkstemp(suffix=".stk")
        command = "xmipp_phantom_create "
        args    = " -i %s"%self.filename
        args += " -o %s"%volFile
        runJob(None, command, args,env=getEnviron())

        #create projections
        command = "xmipp_phantom_project "
        args    = " -i %s"%volFile
        args += " -o %s"%projFile
        args += " --params %s"%self.filenameImages
        runJob(None, command, args,env=getEnviron())

        #os.remove(self.filename)
        #os.remove(outputFile)
