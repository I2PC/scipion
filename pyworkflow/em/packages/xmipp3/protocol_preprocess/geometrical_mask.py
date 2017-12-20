# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
#                J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.protocol.params import (IntParam, EnumParam, FloatParam,
                                        BooleanParam)
from ..constants import *


class XmippGeometricalMask3D:
    """ Basic class for protocols using geometrical masks 3D"""
    
    def defineParams(self, form, isGeometry, addSize):
        # For geometrical sources
        if addSize:
            form.addParam('size', IntParam, condition=isGeometry, 
                          label="Mask size (px)", 
                          help='Select the mask dimensions in voxels. The mask '
                               'will be size x size x size voxels')

        form.addParam('geo', EnumParam, label='Mask type', default=MASK3D_SPHERE, 
                      choices = ['Sphere', 'Box', 'Crown', 'Cylinder',
                                 'Gaussian', 'Raised cosine', 'Raised crown'],
                      condition=isGeometry)
        # TODO add wizard
        form.addParam('radius', IntParam, default=-1, 
                      condition='(geo==%d or geo==%d) and %s'
                                % (MASK3D_SPHERE, MASK3D_CYLINDER, isGeometry),
                      label="Radius (px)",
                      help="Mask radius, if -1, the radius will be MaskSize/2")
        
        form.addParam('shiftCenter', BooleanParam, default=False,
                      label='Shift center of the mask?')
        
        line = form.addLine('Shift Center (px) :', condition='shiftCenter',
                            help='Shift Mask Center to a new origin.')
        line.addParam('centerX', IntParam, default=0, label='X')
        line.addParam('centerY', IntParam, default=0, label='Y')
        line.addParam('centerZ', IntParam, default=0, label='Z')
        
        form.addParam('boxSize', IntParam, default=-1,
                      condition='geo==%d and %s' % (MASK3D_BOX,isGeometry),
                      label="Box size",
                      help="Mask box size, if -1, the box size will be "
                           "MaskSize/2")
        radiusCondition = '(geo==%d or geo==%d or geo==%d) and %s' % (MASK3D_CROWN,
                           MASK3D_RAISED_COSINE, MASK3D_RAISED_CROWN, isGeometry)
        # TODO add wizard
        form.addParam('innerRadius', IntParam, default=0, 
                      condition=radiusCondition,
                      label="Inner radius (px)", help="Inner radius in pixels")
        # TODO add wizard
        form.addParam('outerRadius', IntParam, default=-1, 
                      condition=radiusCondition, label="Outer radius (px)",
                      help="Outer radius in pixels, if -1, the outer radius "
                           "will be MaskSize/2")
        form.addParam('height', IntParam, default=-1,
                      condition='geo==%d and %s' % (MASK3D_CYLINDER,isGeometry), 
                      label="Height (px)", 
                      help="Cylinder height in pixels. If -1, height will be "
                           "MaskSize")
        form.addParam('sigma', FloatParam, default=-1, label="Sigma (px)",
                      condition='geo==%d and %s' % (MASK3D_GAUSSIAN,isGeometry),
                      help="Gaussian sigma in pixels. If -1, sigma will be "
                           "MaskSize/6")
        form.addParam('borderDecay', IntParam, default=0,
                      label="Border decay (px)",
                      condition='geo==%d and %s' % (MASK3D_RAISED_CROWN,isGeometry),
                      help="This is the fall-off of the two borders of the crown")        

    def argsForTransformMask(self,size):
        # Create empty volume file with desired dimensions
        geo = self.geo.get()
        
        # Create the mask
        args = ' --mask '
        r = self.radius.get()
        if r == -1:
            r = size / 2
        r = -r
        iR, oR = self.innerRadius.get(), self.outerRadius.get()
        if oR == -1:
            oR = size / 2
        iR*=-1
        oR*=-1
        
        if geo == MASK3D_SPHERE:
            args += 'circular %d' % r
        elif geo == MASK3D_BOX:
            b = self.boxSize.get()
            if b==-1:
                b = size/2
            args += 'rectangular %d %d %d' % (-b, -b, -b)
        elif geo == MASK3D_CROWN:
            args += 'crown %d %d' % (iR, oR)
        elif geo == MASK3D_CYLINDER:
            height = self.height.get()
            if height==-1:
                height = size
            args += 'cylinder %d %d' % (r, -height)
        elif geo == MASK3D_GAUSSIAN:
            sigma = self.sigma.get()
            if sigma==-1:
                sigma = size/6.0
            args += 'gaussian %f' % -sigma                        
        elif geo == MASK3D_RAISED_COSINE:
            args += 'raised_cosine %d %d' % (iR, oR)
        elif geo == MASK3D_RAISED_CROWN:
            args += 'raised_crown %d %d %d' % (iR, oR, self.borderDecay.get())               
        else:
            raise Exception("No valid mask type value %d" % geo)
        if self.shiftCenter.get():
            args += ' --center %d %d %d'%(self.centerX.get(), self.centerY.get(),
                                          self.centerZ.get())
        return args
    
    def summary(self):
        messages = []      
        geo = self.geo.get()
        if geo == MASK3D_SPHERE:
            messages.append("   Sphere of radius %d"%self.radius.get())
        elif geo == MASK3D_BOX:
            messages.append("   Box of size %d"%self.boxSize.get())
        elif geo == MASK3D_CROWN:
            messages.append("   Crown between %d and %d"
                            % (self.innerRadius.get(), self.outerRadius.get()))
        elif geo == MASK3D_CYLINDER:
            messages.append("   Cylinder of radius %f and height %f"
                            % (self.radius.get(), self.height.get()))
        elif geo == MASK3D_GAUSSIAN:
            messages.append("   Gaussian of sigma %f"%(self.sigma.get()))
        elif geo == MASK3D_RAISED_COSINE:
            messages.append("   Raised cosine between %f and %f"
                            % (self.innerRadius.get(), self.outerRadius.get()))
        elif geo == MASK3D_RAISED_CROWN:
            messages.append("   Raised crown between %f and %f (decay=%f)"
                            % (self.innerRadius.get(), self.outerRadius.get(),
                               self.borderDecay.get()))
        return messages

    def methods(self):
        messages = []      
        geo = self.geo.get()
        if geo == MASK3D_SPHERE:
            messages.append("The mask represented a sphere of radius %d. "
                             % self.radius.get())
        elif geo == MASK3D_BOX:
            messages.append("The mask represented a box of size %d. "
                             % self.boxSize.get())
        elif geo == MASK3D_CROWN:
            messages.append("The mask represented a crown between %d and %d. "
                             % (self.innerRadius.get(), self.outerRadius.get()))
        elif geo == MASK3D_CYLINDER:
            messages.append("The mask represented a cylinder of radius %f "
                            "and height %f. " % (self.radius.get(),
                                                 self.height.get()))
        elif geo == MASK3D_GAUSSIAN:
            messages.append("The mask represented a Gaussian of sigma %f. "
                             % (self.sigma.get()))
        elif geo == MASK3D_RAISED_COSINE:
            messages.append("The mask represented a raised cosine between %f "
                            "and %f. " % (self.innerRadius.get(),
                                          self.outerRadius.get()))
        elif geo == MASK3D_RAISED_CROWN:
            messages.append("The mask represented a raised crown between %f "
                            "and %f (decay=%f)" % (self.innerRadius.get(),
                                                   self.outerRadius.get(),
                                                   self.borderDecay.get()))
        return messages


class XmippGeometricalMask2D:
    """ Basic class for protocols using geometrical masks2D """
    
    def defineParams(self, form, isGeometry, addSize):
        # For geometrical sources
        if addSize:
            form.addParam('size', IntParam, condition=isGeometry,
                          label="Mask size (px)", 
                          help='Select the mask dimensions in pixels. The mask '
                               'will be size x size pixels')

        form.addParam('geo', EnumParam, label='Mask type',
                      default=MASK2D_CIRCULAR, condition=isGeometry, 
                      choices = ['Circular', 'Box', 'Crown', 'Gaussian', 
                                 'Raised cosine', 'Raised crown'])
        # TODO add wizard
        form.addParam('radius', IntParam, default=-1, label="Radius (px)",
                      condition='geo==%d and %s' %(MASK2D_CIRCULAR, isGeometry),
                      help="Mask radius, if -1, the radius will be MaskSize/2")
        form.addParam('boxSize', IntParam, default=-1, label="Box size", 
                      condition='geo==%d and %s' % (MASK2D_BOX, isGeometry),
                      help="Mask box size, if -1, the box size will be MaskSize/2")
        radiusCondition = '(geo==%d or geo==%d or geo==%d) and %s' % (MASK2D_CROWN, 
                              MASK2D_RAISED_COSINE, MASK2D_RAISED_CROWN, isGeometry)
        # TODO add wizard
        form.addParam('innerRadius', IntParam, default=0, 
                      condition=radiusCondition,
                      label="Inner radius (px)", help="Inner radius in pixels")
        # TODO add wizard
        form.addParam('outerRadius', IntParam, default=-1, 
                      condition=radiusCondition, label="Outer radius (px)", 
                      help='Outer radius in pixels, if -1, '
                           'the outer radius will be MaskSize/2')
        form.addParam('sigma', FloatParam, default=-1, label="Sigma (px)", 
                      condition='geo==%d and %s' % (MASK2D_GAUSSIAN,isGeometry),
                      help="Gaussian sigma in pixels. If -1, sigma will be MaskSize/6")                
        form.addParam('borderDecay', IntParam, default=0, label="Border decay (px)", 
                      condition='geo==%d and %s' % (MASK2D_RAISED_CROWN,isGeometry),
                      help="This is the fall-off of the two borders of the crown")
        form.addParam('shiftCenter', BooleanParam, default=False, 
                      label='Shift Center',
                      help='Shift Mask Center to a new origin.')
        form.addParam('centerX', IntParam, default=0, label='X center offset',
                      condition='shiftCenter', help='New x center coordinate')
        form.addParam('centerY', IntParam, default=0, label='Y center offset', 
                      condition='shiftCenter', help='New y center coordinate')

    def argsForTransformMask(self,size):
        # Create empty volume file with desired dimensions
        geo = self.geo.get()
        
        # Create the mask
        args = ' --mask '
        r = self.radius.get()
        if r == -1:
            r = size / 2
        r = -r
        iR, oR = self.innerRadius.get(), self.outerRadius.get()
        if oR == -1:
            oR = size / 2
        iR*=-1
        oR*=-1
        
        if geo == MASK2D_CIRCULAR:
            args += 'circular %d' % r
        elif geo == MASK2D_BOX:
            b = self.boxSize.get()
            if b==-1:
                b = size/2
            args += 'rectangular %d %d' % (-b, -b)
        elif geo == MASK2D_CROWN:
            args += 'crown %d %d' % (iR, oR)
        elif geo == MASK2D_GAUSSIAN:
            sigma = self.sigma.get()
            if sigma==-1:
                sigma = size/6.0
            args += 'gaussian %f' % -sigma                        
        elif geo == MASK2D_RAISED_COSINE:
            args += 'raised_cosine %d %d' % (iR, oR)
        elif geo == MASK2D_RAISED_CROWN:
            args += 'raised_crown %d %d %d' % (iR, oR, self.borderDecay.get())               
        else:
            raise Exception("No valid mask type value %d" % geo)
        if self.shiftCenter.get():
            args += ' --center %d %d' % (self.centerX.get(), self.centerY.get()) 
        return args
    
    def summary(self):
        messages = []      
        geo = self.geo.get()
        if geo == MASK2D_CIRCULAR:
            messages.append("   Circle of radius %d"%self.radius.get())
        elif geo == MASK2D_BOX:
            messages.append("   Box of size %d"%self.boxSize.get())
        elif geo == MASK2D_CROWN:
            messages.append("   Crown between %d and %d" 
                             % (self.innerRadius.get(), self.outerRadius.get()))
        elif geo == MASK2D_GAUSSIAN:
            messages.append("   Gaussian of sigma %f" % (self.sigma.get()))
        elif geo == MASK2D_RAISED_COSINE:
            messages.append("   Raised cosine between %f and %f" 
                             % (self.innerRadius.get(), self.outerRadius.get()))
        elif geo == MASK2D_RAISED_CROWN:
            messages.append("   Raised crown between %f and %f (decay=%f)"
                             % (self.innerRadius.get(), self.outerRadius.get(),
                                self.borderDecay.get()))
        return messages

    def methods(self):
        messages = []      
        geo = self.geo.get()
        if geo == MASK2D_CIRCULAR:
            messages.append("The mask represented a sphere of radius %d. "
                             % self.radius.get())
        elif geo == MASK2D_BOX:
            messages.append("The mask represented a box of size %d. "
                             % self.boxSize.get())
        elif geo == MASK2D_CROWN:
            messages.append("The mask represented a crown between %d and %d. "
                             % (self.innerRadius.get(), self.outerRadius.get()))
        elif geo == MASK2D_GAUSSIAN:
            messages.append("The mask represented a Gaussian of sigma %.2f. "
                             % (self.sigma.get()))
        elif geo == MASK2D_RAISED_COSINE:
            messages.append("The mask represented a raised cosine between %.2f and %.2f. "
                             % (self.innerRadius.get(), self.outerRadius.get()))
        elif geo == MASK2D_RAISED_CROWN:
            messages.append("The mask represented a raised crown between %.2f and %.2f (decay=%.2f)"
                             % (self.innerRadius.get(), self.outerRadius.get(),
                                self.borderDecay.get()))
        return messages
    
