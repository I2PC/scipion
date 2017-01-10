# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow
import pyworkflow.object as pwobj
from pyworkflow.em import *  
from xmipp import MetaData, MDL_ANGLE_ROT, MDL_ANGLE_TILT
from pyworkflow.em.packages.xmipp3.convert import readSetOfParticles
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class XmippProtCreateGallery(ProtAnalysis3D):
    """
    Create a gallery of projections from a volume.
    This gallery of projections may help to understand the images
    observed in the microscope.
    """
    _label = 'create gallery'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input volume')
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See'
                           'http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry '
                           'for a description of the symmetry groups format. '
                           'If no symmetry is present, give c1')

        rot = form.addLine('Rotational angle',
                           help='Minimum, maximum and step values for '
                                'rotational angle range, all in degrees.')
        rot.addParam('rot0', FloatParam, default=0, label='Min')
        rot.addParam('rotF', FloatParam, default=360, label='Max')
        rot.addParam('rotStep', FloatParam, default=5, label='Step')

        tilt = form.addLine('Tilt angle',
                            help='In degrees. tilt=0 is a top view, '
                                 'while tilt=90 is a side view"')
        tilt.addParam('tilt0', FloatParam, default=0, label='Min')
        tilt.addParam('tiltF', FloatParam, default=180, label='Max')
        tilt.addParam('tiltStep', FloatParam, default=5, label='Step')

        form.addParam('maxFreq',FloatParam, default=0.25,
                      expertLevel=LEVEL_ADVANCED,
                      label='Maximum frequency', help="Normalized to 0.5")

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('copyInput')
        self._insertFunctionStep('createGallery')
        self._insertFunctionStep('createOutput')
    
    #--------------------------- STEPS functions -------------------------------
    def copyInput(self):
        ImageHandler().convert(self.inputVolume.get(),
                               self._getTmpPath("volume.vol"))
                        
    def createGallery(self):
        xdim = self.inputVolume.get().getXDim()
        rotN = round((self.rotF.get()-self.rot0.get())/self.rotStep.get())
        tiltN = round((self.tiltF.get()-self.tilt0.get())/self.tiltStep.get())

        paramContent ="""# XMIPP_STAR_1 *
data_block1
_dimensions2D   '%d %d'
_projRotRange    '%f %f %d'
_projRotRandomness   even 
_projTiltRange    '%f %f %d'
_projTiltRandomness   even 
_projPsiRange    '0 0 1'
_projPsiRandomness   even 
""" % (xdim, xdim, self.rot0, self.rotF,rotN, self.tilt0, self.tiltF, tiltN)
        fhParam = open(self._getExtraPath("projectionParameters.xmd"), 'w')
        fhParam.write(paramContent)
        fhParam.close()

        self.runJob("xmipp_phantom_project",
                    "-i %s -o %s --params %s --method fourier 2 %f --sym %s" %
                    (self._getTmpPath("volume.vol"),
                     self._getPath("images.stk"),
                     self._getExtraPath("projectionParameters.xmd"),
                     self.maxFreq, self.symmetryGroup))

    def createOutput(self):
        imgSetOut = self._createSetOfAverages()
        imgSetOut.setSamplingRate(self.inputVolume.get().getSamplingRate())
        imgSetOut.setAlignmentProj()
        readSetOfParticles(self._getPath("images.xmd"), imgSetOut)

        self._defineOutputs(outputReprojections=imgSetOut)
        self._defineSourceRelation(self.inputVolume, imgSetOut)

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        messages = []
        messages.append("Rot.angle from %0.2f to %0.2f in steps of %0.2f" %
                        (self.rot0, self.rotF, self.rotStep))
        messages.append("Tilt.angle from %0.2f to %0.2f in steps of %0.2f" %
                        (self.tilt0, self.tiltF, self.tiltStep))
        return messages
