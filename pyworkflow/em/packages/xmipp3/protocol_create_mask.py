# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains protocols for creating 3D masks.
"""


from pyworkflow.em import *  
from constants import *

      
class XmippProtCreateGeo3DMask(ProtCreateMask3D):
    """ Create a 3D mask from a geometrical description (Sphere, Box, Cylinder...) """

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('size', IntParam, label="Mask size", 
                      help='Select the mask dimensions in voxels.')
        form.addParam('geo', EnumParam, label='Mask type', default=MASK3D_SPHERE,
                      choices = ['Sphere', 'Box', 'Crown', 'Cylinder', 
                                 'Gaussian', 'Raised cosine', 'Raised crown'])
        form.addParam('radius', IntParam, default=-1, 
                      condition='geo==%d or geo==%d' % (MASK3D_SPHERE, MASK3D_CYLINDER),
                      label="Radius", help="Mask radius, if -1, the radius will be MaskSize/2")
        form.addParam('boxSize', IntParam, condition='geo==%d' % MASK3D_BOX,
                      label="Box size", help="Mask radius, if -1, the radius will be MaskSize/2")
        radiusCondition = 'geo==%d or geo==%d or geo==%d' % (MASK3D_CROWN, MASK3D_RAISED_COSINE, MASK3D_RAISED_CROWN)
        form.addParam('innerRadius', IntParam, default=0, 
                      condition=radiusCondition,
                      label="Inner radius (pix)", help="Inner radius in pixels")
        form.addParam('outerRadius', IntParam, default=-1, 
                      condition=radiusCondition,
                      label="Outer radius (pix)", help="Outer radius in pixels")
        form.addParam('height', IntParam, condition='geo==%d' % MASK3D_CYLINDER,
                      label="Height (pix)", help="Cylinder height in pixels")
        form.addParam('sigma', IntParam, condition='geo==%d' % MASK3D_GAUSSIAN,
                      label="Sigma (pix)", help="Cylinder height in pixels")                
        form.addParam('pixelWidth', IntParam, condition='geo==%d' % MASK3D_RAISED_CROWN,
                      label="Pixel width (pix)", help="in pixels")        
        
    def _defineSteps(self):
        self.maskFile = self._getPath('mask.vol')
        self._insertFunctionStep('createMask', self.size.get(), self.geo.get())
        self._insertFunctionStep('createOutput')
        
    def createMask(self, size, geo):
        # Create empty volume file with desired dimensions
        import xmipp
        xmipp.createEmptyFile(self.maskFile, size, size, size)
        # Create the mask
        program = "xmipp_transform_mask"
        args = '-i %s --mask ' % self.maskFile
        r = self.radius.get()
        if r == -1:
            r = size / 2
        r = -r
        iR, oR = -self.innerRadius.get(), -self.outerRadius.get()
        
        if geo == MASK3D_SPHERE:
            args += 'circular %d' % r
        elif geo == MASK3D_BOX:
            b = -(self.boxSize.get())
            args += 'rectangular %d %d %d' % (b, b, b)
        elif geo == MASK3D_CROWN:
            args += 'crown %d %d' % (iR, oR)
        elif geo == MASK3D_CYLINDER:
            args += 'cylinder %d %d' % (r, -self.height.get())
        elif geo == MASK3D_GAUSSIAN:
            args += 'gaussian %f' % -self.sigma.get()                        
        elif geo == MASK3D_RAISED_COSINE:
            args += 'raised_cosine %d %d' % (iR, oR)
        elif geo == MASK3D_RAISED_CROWN:
            args += 'raised_crown %d %d %d' % (iR, oR, self.pixelWidth.get())               
        else:
            raise Exception("No valid mask type value %d" % geo)
        args += ' --create_mask %s' % self.maskFile
        
        self.runJob(None, program, args)
    
    def createOutput(self):
        volMask = VolumeMask()
        volMask.setFileName(self.maskFile)
        self._defineOutputs(outputMask=volMask)
        
    def _summary(self):
        summary = []
        return summary
    
