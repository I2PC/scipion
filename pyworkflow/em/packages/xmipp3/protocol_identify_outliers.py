# **************************************************************************
# *
# * Authors:     Javier Vargas and Adrian Quintana (jvargas@cnb.csic.es aquintana@cnb.csic.es)
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
#from pyworkflow.em.packages.xmipp3.convert import locationToXmipp,\
#    writeSetOfVolumes
from pyworkflow.em.packages.xmipp3.convert import locationToXmipp, createXmippInputImages

"""
This sub-package contains the protocol projection outliers
"""
from pyworkflow.em import *  
from xmipp3 import ProjMatcher
from pyworkflow.utils.path import cleanPath

class XmippProtIdentifyOutliers(ProtClassify3D, ProjMatcher):
    """ Protocol for screening a number of classes comparing them to a volume. """
    _label = 'identify outliers'
    
    def __init__(self, **args):
        ProtClassify3D.__init__(self, **args)
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImages', PointerParam, label="Input particles", important=True, 
              pointerClass='SetOfParticles',
              help='Select the input iamges from the project.'
                   'It should be a SetOfParticles class') 
        form.addParam('inputVolume', PointerParam, label="Volume to compare images to", important=True, 
              pointerClass='Volume',
#              pointerClass='SetOfVolumes',
              help='Provide a volume against which images will be compared'
                   'It should be a Volume class')
        form.addParam('volumeIsCTFCorrected', BooleanParam, label="Volume has been corrected by CTF", default=False,
                      help='Volume has been corrected by CTF')
        form.addParam('symmetryGroup', StringParam, default="c1",
              label='Symmetry group', 
              help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                'If no symmetry is present, give c1')
        form.addParam('angularSampling', FloatParam, default=3, expertLevel=LEVEL_EXPERT,
                      label="Angular sampling",
                      help='In degrees. This sampling defines how fine the projection gallery from the volume is explored.')
         
        form.addParallelSection(mpi=8)
        
    def _insertAllSteps(self):
        # Projection matching
        self.fnAngles = self._getTmpPath('angles.xmd')
#        self.images = locationToXmipp(*self.inputImages.get().getLocation())
        self.images = createXmippInputImages(self, self.inputImages.get())
        
        self._insertFunctionStep("projMatchStep", locationToXmipp(*self.inputVolume.get().getLocation()), self.angularSampling.get(), self.symmetryGroup.get(),\
                                 self.images, self.fnAngles, self.Xdim)        

        # Prepare output
        self.fnOutputImages = self._getPath('images.xmd')
        self.visualizeInfoOutput = String(self.fnOutputImages)
        self._insertRunJobStep("xmipp_metadata_utilities","-i %s --set join %s -o %s" %
                                (self.fnAngles, self.images, self.fnOutputImages), numberOfMpi=1)
        
        # Produce difference images
        fnDiff = self._getExtraPath("diff.stk")
        fnAligned = self._getExtraPath("images_aligned.xmd")
        self._insertFunctionStep("produceAlignedImagesStep", self.volumeIsCTFCorrected.get(),\
                                 fnAligned, self.fnOutputImages)
        
        # Evaluate each image
        fnAutoCorrelations = self._getExtraPath("autocorrelations.xmd")
        self._insertRunJobStep("xmipp_image_residuals", "-i %s -o %s --save_metadata_stack %s" %
                              (fnDiff, self._getExtraPath("autocorrelations.stk"), fnAutoCorrelations),
                              [fnAutoCorrelations], numberOfMpi=1)
        self._insertRunJobStep("xmipp_metadata_utilities", "-i %s --set merge %s" % 
                               (fnAligned,fnAutoCorrelations), numberOfMpi=1)
        
        cleanPath(fnAutoCorrelations)

        # Prepare output
        self._insertRunJobStep("xmipp_metadata_utilities", "-i %s --set join %s" % 
                               (self.fnOutputImages, fnAligned), numberOfMpi=1)
        self._insertRunJobStep("xmipp_metadata_utilities", "-i %s --operate sort zScoreResCov desc" % 
                               self.fnOutputImages, numberOfMpi=1)  
            
    def getVisualizeInfo(self):
        return self.visualizeInfoOutput

    def _summary(self):
        summary = []
        summary.append("Set of particles: [%s] " % self.inputImages.get().getNameId())
        summary.append("Volume: [%s] " % self.inputVolume.getNameId())
        summary.append("Symmetry: %s " % self.symmetryGroup.get())
        summary.append("Volume is CTF corrected: %s " % self.volumeIsCTFCorrected.get())
        return summary
        
    def validate(self):
        self.Xdim = self.inputImages.get().getDimensions()[0]
        volumeXdim = self.inputVolume.get().getDim()[0]
         
        errors = []
        if volumeXdim != self.Xdim:
            errors.append("Make sure that the volume and the images have the same size")
        return errors  
        
    def _citations(self):
       return []
            