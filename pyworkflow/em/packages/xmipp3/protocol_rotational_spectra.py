# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains wrapper around rotational spectra Xmipp program
"""

from pyworkflow.em import *  
import xmipp

import xmipp, xmipp3
from protocol_kerdensom import KendersomBaseClassify

        
class XmippProtRotSpectra(KendersomBaseClassify):
    """Protocol to compute the rotational spectrum of the given particles"""
    _label = 'rotational spectra'
    
    def __init__(self, **args):
        KendersomBaseClassify.__init__(self, **args)
        
    def _addParams(self, form):
        form.addParam('howCenter', EnumParam, choices=['Middle of the image', 'Minimize first harmonic'], 
                      default=xmipp3.ROTSPECTRA_CENTER_MIDDLE, important=True, 
                      label='How to find the center of rotation', display=EnumParam.DISPLAY_COMBO, 
                      help='Select how to find the center of rotation.')  
        form.addParam('spectraInnerRadius', IntParam, default=15,
                      label='Inner radius for rotational harmonics (%)',
                      help='A percentage of the image radius', expertLevel=LEVEL_ADVANCED)
        form.addParam('spectraOuterRadius', IntParam, default=80,
                      label='Outer radius for rotational harmonics (%)',
                      help='A percentage of the image radius', expertLevel=LEVEL_ADVANCED)
        form.addParam('spectraLowHarmonic', IntParam, default=1,
                      label='Lowest harmonic to calculate',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('spectraHighHarmonic', IntParam, default=15,
                      label='Highest harmonic to calculate',
                      expertLevel=LEVEL_ADVANCED) 
        form.addSection(label='Kerdensom')
        
    def _defineParamsStep(self):
        KendersomBaseClassify._defineParamsStep(self)
        self._params['extraDir'] = self._getExtraPath()
        self._params['R1'] = self.spectraInnerRadius.get()
        self._params['R2'] = self.spectraOuterRadius.get()
        self._params['spectraLowHarmonic'] = self.spectraLowHarmonic.get()
        self._params['spectraHighHarmonic'] = self.spectraHighHarmonic.get()
        self._params['vectors'] = self._getExtraPath("rotSpectra.xmd")
    
    def _defineProccesStep(self):
        imagesFn = self._params['imgsFn']
        centerFn = self._getExtraPath("center2d_center.xmd")
        # After any of this steps the file "center2d_center.xmd" should be produced
        if self.howCenter == xmipp3.ROTSPECTRA_CENTER_MIDDLE:
            self._insertMiddleStep(imagesFn, centerFn)
        else:
            self._insertFunctionStep('centerFirstHarmonicStep', imagesFn, centerFn)
        # Produce "rotSpectra.xmd" vectors
        self._insertFunctionStep('calculateSpectraStep', imagesFn, centerFn, self._params['vectors'])
        # Call kerdensom for classification
        self._insertKerdensomStep()
    
    def _insertMiddleStep(self, inputImages, outputCenter):
        R2 = self._params['R2']
        
        if R2 + 20 > 100:
            R3 = R2 + (100 - R2) / 2
            R4 = 100
        else:
            R3 = R2 + 10
            R4 = R2 + 20
        self._params['R3'] = R3
        self._params['R4'] = R4
        
        program = 'xmipp_image_find_center'
        args = '-i ' + inputImages
        args += ' --oroot %(extraDir)s/center2d --r1 %(R1)d --r2 %(R2)d --r3 %(R3)d --r4 %(R4)d'
        # Run the command with formatted parameters
        self._insertRunJobStep(program, args % self._params)
            
    def centerFirstHarmonicStep(self, inputImages, outputCenter):
        dims = xmipp.MetaDataInfo(str(inputImages))
        md = xmipp.MetaData()
        objId = md.addObject()
        md.setValue(xmipp.MDL_X, float(dims[0] / 2), objId)
        md.setValue(xmipp.MDL_Y, float(dims[1] / 2), objId)
        md.write(outputCenter)
        return [outputCenter] # this file should exists after the step
            
    def calculateSpectraStep(self, inputImages, inputCenter, outputSpectra):     
        md = xmipp.MetaData(inputCenter)
        objId = md.firstObject()
        self._params['xOffset'] = md.getValue(xmipp.MDL_X, objId)
        self._params['yOffset'] = md.getValue(xmipp.MDL_Y, objId)
        
        program = 'xmipp_image_rotational_spectra'
        args = "-i %s -o %s" % (inputImages, outputSpectra)
        args += ' --x0 %(xOffset)d --y0 %(yOffset)d --r1 %(R1)d --r2 %(R2)d' + \
                     ' --low %(spectraLowHarmonic)d --high %(spectraHighHarmonic)d'
        # Run the command with formatted parameters
        self.runJob(program, args % self._params)
        return [outputSpectra]
         
    def _citations(self):
        return ['Montano2000']
        