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

from os.path import join
import numpy as np

import pyworkflow.em as em  
from pyworkflow.utils.path import makePath
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol.params import EnumParam, IntParam
import xmipp, xmipp3

from protocol_kerdensom import KendersomBaseClassify
from convert import readSetOfClasses2D

        
class XmippProtRotSpectra(KendersomBaseClassify):
    """Protocol to compute the rotational spectrum of the given particles.
    """
    _label = 'rotational spectra'
    
    def __init__(self, **args):
        KendersomBaseClassify.__init__(self, **args)
        
    #--------------------------- DEFINE param functions --------------------------------------------
    def _addParams(self, form):
        form.addSection(label='Spectra')
        form.addParam('howCenter', EnumParam, 
                      choices=['Middle of the image', 'Minimize first harmonic'], 
                      default=xmipp3.ROTSPECTRA_CENTER_FIRST_HARMONIC, 
                      display=EnumParam.DISPLAY_COMBO, 
                      label='How to find the center of rotation', important=True,  
                      help='Select how to find the center of rotation.')
        
        line = form.addLine('Rotational harmonics radius (px)')
        line.addParam('spectraInnerRadius', IntParam, default=5,
                      label='Inner',
                      help='A percentage of the image radius')
        line.addParam('spectraOuterRadius', IntParam, default=43,
                      label='Outer',
                      help='A percentage of the image radius')
        
        line = form.addLine('Harmonics to calculate')
        line.addParam('spectraLowHarmonic', IntParam, default=1,
                      label='Lowest')
        line.addParam('spectraHighHarmonic', IntParam, default=15,
                      label='Highest') 
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _prepareParams(self):
        KendersomBaseClassify._prepareParams(self)
        self._params['extraDir'] = self._getExtraPath()
        # we need to convert R1 and R2 to percentage of the radius
        radius = self.inputParticles.get().getDim()[0] / 2
        percent = 100. / radius
        self._params['R1'] = self.spectraInnerRadius.get() * percent
        self._params['R2'] = self.spectraOuterRadius.get() * percent
        self._params['spectraLowHarmonic'] = self.spectraLowHarmonic.get()
        self._params['spectraHighHarmonic'] = self.spectraHighHarmonic.get()
        self._params['vectors'] = self._getExtraPath("rotSpectra.xmd")
    
    def _insertProccesStep(self):
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
    
    def _insertMiddleStep(self, imagesFn, outputCenter):
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
        args = '-i ' + imagesFn
        args += ' --oroot %(extraDir)s/center2d --r1 %(R1)d --r2 %(R2)d --r3 %(R3)d --r4 %(R4)d'
        # Run the command with formatted parameters
        self._insertRunJobStep(program, args % self._params)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def centerFirstHarmonicStep(self, imagesFn, outputCenter):
        dims = xmipp.MetaDataInfo(str(imagesFn))
        md = xmipp.MetaData()
        objId = md.addObject()
        md.setValue(xmipp.MDL_X, float(dims[0] / 2), objId)
        md.setValue(xmipp.MDL_Y, float(dims[1] / 2), objId)
        md.write(outputCenter)
        return [outputCenter] # this file should exists after the step
            
    def calculateSpectraStep(self, imagesFn, inputCenter, outputSpectra):     
        md = xmipp.MetaData(inputCenter)
        objId = md.firstObject()
        self._params['xOffset'] = md.getValue(xmipp.MDL_X, objId)
        self._params['yOffset'] = md.getValue(xmipp.MDL_Y, objId)
        
        program = 'xmipp_image_rotational_spectra'
        args = "-i %s -o %s" % (imagesFn, outputSpectra)
        args += ' --x0 %(xOffset)s --y0 %(yOffset)s --r1 %(R1)d --r2 %(R2)d' + \
                     ' --low %(spectraLowHarmonic)d --high %(spectraHighHarmonic)d'
        # Run the command with formatted parameters
        self.runJob(program, args % self._params)
        return [outputSpectra]
    
    def createOutputStep(self):
        self.plotsDir = self._getExtraPath('plots')
        makePath(self.plotsDir)
        
        fnClassVectors = self._params['kvectors'].replace('xmd', 'vec')
        f = open(fnClassVectors)
        self.classArray = np.fromfile(f, dtype=np.float32)
        f.close()
        
        fnImgVectors = self._params['vectors'].replace('xmd', 'vec')
        f = open(fnImgVectors)
        self.imgArray = np.fromfile(f, dtype=np.float32)
        f.close()
        self.imgCount = 0
        
        imgSet = self.inputParticles.get()
        classes2DSet = self._createSetOfClasses2D(imgSet)
        readSetOfClasses2D(classes2DSet, self._params['kclasses'], 
                           preprocessClass=self._preprocessClass,
                           postprocessImageRow=self._postprocessImageRow)
        self._defineOutputs(outputClasses=classes2DSet)
        self._defineSourceRelation(self.inputParticles, classes2DSet)
                
    def _preprocessClass(self, classItem, classRow):
        KendersomBaseClassify._preprocessClass(self, classItem, classRow)
        ref = classRow.getValue(xmipp.MDL_REF) # get class number
        classItem.spectraPlot = em.Image()
        classItem.spectraPlot.setFileName(self._createSpectraPlot('class', 
                                                                  self.classArray, 
                                                                  ref))
        
    def _postprocessImageRow(self, img, imgRow):
        self.imgCount += 1
        img.spectraPlot = em.Image()
        img.spectraPlot.setFileName(self._createSpectraPlot('image', 
                                                            self.imgArray, 
                                                            self.imgCount, 
                                                            img.getObjId()))
        
    def _createSpectraPlot(self, label, array, index, objId=None):
        if objId is None:
            objId = index
        # Number of harmonics calculated
        plotter = Plotter()
        a = plotter.createSubPlot('Spectrum for %s %d' % (label, objId), 
                                  xlabel='Harmonics', ylabel='Energy')
        
        n = self.spectraHighHarmonic.get() - self.spectraLowHarmonic.get() + 1
        i1 = (index-1) * n
        i2 = i1 + n
        xdata = range(self.spectraLowHarmonic.get(), self.spectraHighHarmonic.get()+1)
        ydata = array[i1:i2]
        a.plot(xdata, ydata)
        
        plotFile = join(self.plotsDir, 'spectra_%s%06d.png' % (label, objId))
        plotter.savefig(plotFile)
        
        return plotFile
        
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        return KendersomBaseClassify._validate(self)
    
    def _citations(self):
        return ['Pascual2000']
    
    def _summary(self):
        return KendersomBaseClassify._summary(self)
    
    def _methods(self):
        return KendersomBaseClassify._methods(self)
   
