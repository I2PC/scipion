# **************************************************************************
# *
# * Authors:     Laura del Cano (laura.cano@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
This sub-package contains wrapper around Screen Particles Xmipp program
"""

from pyworkflow.em import *  
import xmipp

from convert import createXmippInputImages, readSetOfParticles

# Automatic Particle rejection enum
REJ_NONE = 0
REJ_MAXZSCORE = 1
REJ_PERCENTAGE =2
        
        
class XmippProtScreenParticles(ProtProcessParticles):
    """ Protocol to screen a set of particles in the project. """
    _label = 'screen particles'

    def _defineProcessParams(self, form):
        
        form.addParam('autoParRejection', EnumParam, choices=['None', 'MaxZscore', 'Percentage'],
                      label="Automatic particle rejection", default=REJ_NONE,
                      display=EnumParam.DISPLAY_COMBO, expertLevel=LEVEL_EXPERT,
                      help='How to automatically reject particles. It can be none (no rejection), '
                      'maxZscore (reject a particle if its Zscore is larger than this value), '
                      'Percentage (reject a given percentage in each one of the screening criteria). ')
        form.addParam('maxZscore', IntParam, default=3, condition='autoParRejection==1',
                      label='Maximum Zscore', expertLevel=LEVEL_EXPERT,
                      help='Maximum Zscore.', validators=[Positive])      
        form.addParam('percentage', IntParam, default=5, condition='autoParRejection==2',
                      label='Percentage (%)', expertLevel=LEVEL_EXPERT,
                      help='Percentage.', validators=[Range(0, 100, error="Percentage must be between 0 and 100.")])        
        
                
    def _defineSteps(self):
        """ Mainly prepare the command line for call cl2d program"""
        
        # Convert input images if necessary
        imgsFn = createXmippInputImages(self, self.inputParticles.get())
        
        self._insertFunctionStep('sortImages', imgsFn) 
        
        self._insertFunctionStep('createOutput')

    def sortImages(self, inputFile):
        from xmipp import MetaDataInfo
        #(Xdim, Ydim, Zdim, Ndim, _) = MetaDataInfo(inputFile)
        args=""
        # copy file to run path
        self.outputMd = String(self._getPath(replaceBaseExt(inputFile, 'xmd')))
        self.outputMd._objDoStore = True
        if inputFile != self.outputMd.get():
            copyFile(inputFile, self.outputMd.get())
        if self.autoParRejection.get()==REJ_MAXZSCORE:
            args+=" --zcut "+str(self.maxZscore.get())
        elif self.autoParRejection.get()==REJ_PERCENTAGE:
            args+=" --percent "+str(self.percentage.get())
        #if Ndim > 0:
        self.runJob(None, "xmipp_image_sort_by_statistics", "-i " + self.outputMd.get() + " --addToInput"+args)

                        
    def createOutput(self):
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        readSetOfParticles(self.outputMd.get(), imgSet, imgSet.hasCTF())
        imgSet.write()

        self._defineOutputs(outputParticles=imgSet)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Input Images: %s" % self.inputParticles.get().getNameId())
            summary.append("Output particles: %s" % self.outputParticles.get())
        return summary
