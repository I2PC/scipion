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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains wrapper around reconstruct_significant Xmipp program
"""

from pyworkflow.utils import Timer
from pyworkflow.em import *  
from pyworkflow.em.packages.xmipp3.convert import writeSetOfVolumes, volumeToRow
from pyworkflow.em.packages.xmipp3.xmipp3 import XmippMdRow
from convert import writeSetOfClasses2D, writeSetOfParticles
import pyworkflow.em.metadata as metadata



class XmippProtReconstructSignificant(ProtInitialVolume):
    """ 
    This algorithm addresses the initial volume problem in SPA
    by setting it in a Weighted Least Squares framework and 
    calculating the weights through a statistical approach based on
    the cumulative density function of different image similarity measures. 
    """    
    _label = 'reconstruct significant'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, label="Input classes", important=True, 
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      help='Select the input classes2D from the project.\n'
                           'It should be a SetOfClasses2D class with class representative')
        form.addParam('symmetryGroup', TextParam, default='c1',
                      label="Symmetry group",
                      help='See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')  
        form.addParam('thereisRefVolume', BooleanParam, default=False,
                      label="Is there a reference volume(s)?", 
                       help='You may use a reference volume to initialize the calculations. For instance, '
                            'this is very useful to obtain asymmetric volumes from symmetric references. The symmetric '
                            'reference is provided as starting point, choose no symmetry group (c1), and reconstruct_significant'
                            'will tend to break the symmetry finding a suitable volume. The reference volume can also be useful, '
                            'for instance, when reconstructing a fiber. Provide in this case a cylinder of a suitable size.')
        form.addParam('refVolume', PointerParam, label='Initial 3D reference volumes',
                      pointerClass='SetOfVolumes, Volume', condition="thereisRefVolume")
        form.addParam('Nvolumes', IntParam, label='Number of volumes', help="Number of volumes to reconstruct",
                      default=1,condition="not thereisRefVolume")
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label='Angular sampling',
                      help='Angular sampling in degrees for generating the projection gallery.')
        form.addParam('minTilt', FloatParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label='Minimum tilt (deg)',
                      help='Use the minimum and maximum tilts to limit the angular search. This can be useful, for instance, '
                           'in the reconstruction of fibers from side views. 0 degrees is a top view, while 90 degrees is a side view.')
        form.addParam('maxTilt', FloatParam, default=90, expertLevel=LEVEL_ADVANCED,
                      label='Maximum tilt (deg)',
                      help='Use the minimum and maximum tilts to limit the angular search. This can be useful, for instance, '
                           'in the reconstruction of fibers from side views. 0 degrees is a top view, while 90 degrees is a side view.')
        form.addParam('maximumShift', FloatParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label='Maximum shift (px):', help="Set to -1 for free shift search")
        form.addParam('keepIntermediate', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Keep intermediate volumes',
                      help='Keep all volumes and angular assignments along iterations')

        form.addSection(label='Criteria')
        form.addParam('alpha0', FloatParam, default=80,
                      label='Starting significance',
                      help='80 means 80% of significance. Use larger numbers to relax the starting significance and have a smoother '
                           'landscape of solutions')
        form.addParam('iter', IntParam, default=50,
                      label='Number of iterations',
                      help='Number of iterations to go from the initial significance to the final one')
        form.addParam('alphaF', FloatParam, default=99.5,
                      label='Final significance',
                      help='99.5 means 99.5% of significance. Use smaller numbers to be more strict and have a sharper reconstruction.'
                           'Be aware that if you are too strict, you may end with very few projections and the reconstruction becomes very'
                           'noisy.')
        form.addParam('useImed', BooleanParam, default=True, expertLevel=LEVEL_ADVANCED,
                      label='Use IMED', help='Use IMED for the weighting. IMED is an alternative to correlation that can '
                      'discriminate better among very similar images')
        form.addParam('strictDir', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Strict direction', help='If the direction is strict, then only the most significant experimental images '
                      'can contribute to it. As a consequence, many experimental classes are lost and only the best contribute to the 3D '
                      'reconstruction. Be aware that only the best can be very few depending on the cases.')
        form.addParam('angDistance', IntParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label='Angular neighborhood', help='Images in an angular neighborhood also determines the weight of each image. '
                      'It should be at least twice the angular sampling')
        form.addParam('dontApplyFisher', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Do not apply Fisher', help="Images are preselected using Fisher's confidence interval on the correlation "
                      "coefficient. Check this box if you do not want to make this preselection.")

        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def getSignificantArgs(self, imgsFn):
        """ Return the arguments needed to launch the program. """
        # Prepare arguments to call program: xmipp_classify_CL2D
        self._params = {'imgsFn': imgsFn, 
                        'extraDir': self._getExtraPath(),
                        'symmetryGroup': self.symmetryGroup.get(),
                        'angularSampling': self.angularSampling.get(),
                        'minTilt': self.minTilt.get(),
                        'maxTilt': self.maxTilt.get(),
                        'maximumShift': self.maximumShift.get(),
                        'angDistance': self.angDistance.get()
                        }
        args = '-i %(imgsFn)s --sym %(symmetryGroup)s --angularSampling %(angularSampling)f '\
               '--minTilt %(minTilt)f --maxTilt %(maxTilt)f --maxShift %(maximumShift)f '\
               '--dontReconstruct --angDistance %(angDistance)f' % self._params
        
        if self.useImed:
            args += " --useImed"
        if self.strictDir:
            args += " --strictDirection"
        if self.dontApplyFisher:
            args += " --dontApplyFisher"
            
        return args
        
    def _insertAllSteps(self):
        # Convert input images if necessary
        self.imgsFn = self._getExtraPath('input_classes.xmd')
        self._insertFunctionStep('convertInputStep', self.imgsFn)
        
        args = self.getSignificantArgs(self.imgsFn)    
        n = self.iter.get()
        alpha0 = self.alpha0.get()
        deltaAlpha = (self.alphaF.get() - alpha0) / n
        
        # Insert one step per iteration
        for i in range(n):
            alpha = 1 - (alpha0 + deltaAlpha * i)/100.0
            self._insertFunctionStep('significantStep', i+1, alpha, args)           
            
        self._insertFunctionStep('createOutputStep')        

    #--------------------------- STEPS functions --------------------------------------------   
    def significantStep(self, iterNumber, alpha, args):
        iterDir = self._getTmpPath('iter%03d' % iterNumber)
        makePath(iterDir)
        args += ' --odir %s' % iterDir
        args += ' --alpha0 %f --alphaF %f' % (alpha, alpha)
        prevVolFn = self.getIterVolume(iterNumber-1)
        volFn = self.getIterVolume(iterNumber)
        
        if iterNumber == 1:
            if self.thereisRefVolume:
                args += " --initvolumes " + self._getExtraPath('input_volumes.xmd')
            else:
                args += " --numberOfVolumes %d" % self.Nvolumes
        else:
            args += " --initvolumes %s" % prevVolFn
        
        t = Timer()
        t.tic()
        self.runJob("xmipp_reconstruct_significant", args)
        t.toc('Significant took: ')
        
        anglesFn = self._getExtraPath('angles_iter%03d.xmd' % iterNumber)
        moveFile(os.path.join(iterDir, 'angles_iter001_00.xmd'), anglesFn)
        reconsArgs = ' -i %s' %  anglesFn
        reconsArgs += ' -o %s' % volFn
        reconsArgs += ' --weight -v 0  --sym %s ' % self.symmetryGroup

        print "Number of images for reconstruction: ", metadata.getSize(anglesFn)
        t.tic()
        self.runJob("xmipp_reconstruct_fourier", reconsArgs)
        t.toc('Reconstruct fourier took: ')
        
        xdim = self.inputSet.get().getDimensions()[0]
        maskArgs = "-i %s --mask circular %d -v 0" % (volFn, -xdim/2)
        self.runJob('xmipp_transform_mask', maskArgs, numberOfMpi=1)
        
        if not self.keepIntermediate:
            cleanPath(prevVolFn, iterDir)
            
        
    def convertInputStep(self, classesFn):
        inputSet = self.inputSet.get()
        
        if isinstance(inputSet, SetOfClasses2D):
            writeSetOfClasses2D(inputSet, classesFn, writeParticles=False)
        else:
            writeSetOfParticles(inputSet, classesFn)
            
        if self.thereisRefVolume:
            inputVolume= self.refVolume.get()
            fnVolumes = self._getExtraPath('input_volumes.xmd')
            if isinstance(inputVolume, SetOfVolumes):
                writeSetOfVolumes(inputVolume, fnVolumes)
            else:
                row = XmippMdRow()
                volumeToRow(inputVolume, row, alignType = ALIGN_NONE)
                md = xmipp.MetaData()
                row.writeToMd(md, md.addObject())
                md.write(fnVolumes)
        
    def createOutputStep(self):
        Nvolumes = self.getNumberOfVolumes()
        lastIter = self.getLastIteration(Nvolumes)
        if Nvolumes==1:
            vol = Volume()
            vol.setObjComment('significant volume 1')
            vol.setLocation(self.getIterVolume(lastIter))
            vol.setSamplingRate(self.inputSet.get().getSamplingRate())
            self._defineOutputs(outputVolume=vol)
            output = vol
        else:
            volSet = self._createSetOfVolumes()
            volSet.setSamplingRate(self.inputSet.get().getSamplingRate())
            fnVolumes = glob(self._getExtraPath('volume_iter%03d_*.vol')%lastIter)
            fnVolumes.sort()
            for i, fnVolume in enumerate(fnVolumes):
                vol = Volume()
                vol.setObjComment('significant volume %02d' % (i+1))
                vol.setLocation(fnVolume)
                volSet.append(vol)
            self._defineOutputs(outputVolumes=volSet)
            output = volSet

        self._defineSourceRelation(self.inputSet, output)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.thereisRefVolume:
            if self.refVolume.hasValue():
                refVolume = self.refVolume.get()
                x1, y1, _ = refVolume.getDim()
                x2, y2, _ = self.inputSet.get().getDimensions()
                if x1!=x2 or y1!=y2:
                    errors.append('The input images and the reference volume have different sizes') 
            else:
                errors.append("Please, enter a reference image")
        return errors
        
    def _summary(self):
        summary = []
        summary.append("Input classes: %s" % self.getObjectTag('inputSet'))
        if self.thereisRefVolume:
            summary.append("Starting from: %s" % self.getObjectTag('refVolume'))
        else:
            summary.append("Starting from: %d random volumes" % self.Nvolumes)
        summary.append("Significance from %f%% to %f%% in %d iterations" % (self.alpha0, self.alphaF, self.iter))
        if self.useImed:
            summary.append("IMED used")
        if self.strictDir:
            summary.append("Strict directions")    
        return summary
    
    def _citations(self):
        return ['Sorzano2015']
    
    def _methods(self):
        retval = ""
        if self.inputSet.get() is not None:
            retval = "We used reconstruct significant to produce an initial volume "
            retval += "from the set of classes %s." % self.getObjectTag('inputSet')
            if self.thereisRefVolume:
                retval += " We used %s volume " % self.getObjectTag('refVolume')
                retval += "as a starting point of the reconstruction iterations."
            else:
                retval += " We started the iterations with %d random volumes." % self.Nvolumes
            retval += " %d iterations were run going from a " % self.iter
            retval += "starting significance of %f%% to a final one of %f%%." % (self.alpha0, self.alphaF)
            if self.useImed:
                retval += " IMED weighting was used."
            if self.strictDir:
                retval += " The strict direction criterion was employed." 
        
        if self.getNumberOfVolumes() > 1:
            if self.hasAttribute('outputVolumes'):
                retval += " The set of reconstructed volumes was %s." % self.getObjectTag('outputVolumes')
        else:
            if self.hasAttribute('outputVolume'):
                retval+=" The reconstructed volume was %s." % self.getObjectTag('outputVolume')
        return [retval]


    #--------------------------- UTILS functions --------------------------------------------

    def getIterVolume(self, iterNumber):
        return self._getExtraPath('volume_iter%03d.vol' % iterNumber)
    
    def getIterTmpVolume(self, iterNumber):
        self._getTmpPath('iter%03d' % iterNumber, 'volume_iter001.vol')
        
    def getLastIteration(self,Nvolumes):
        lastIter =-1
        for n in range(1, self.iter.get()+1):
            NvolumesIter=len(glob(self._getExtraPath('volume_iter%03d*.vol' % n)))
            if NvolumesIter==0:
                continue
            elif NvolumesIter==Nvolumes:
                lastIter=n
            else:
                break
        return lastIter
    
    def getNumberOfVolumes(self):
        if self.thereisRefVolume:
            inputVolume= self.refVolume.get()
            if isinstance(inputVolume, SetOfVolumes):
                Nvolumes=inputVolume.getSize()
            else:
                Nvolumes=1
        else:
            Nvolumes=self.Nvolumes.get()
        return Nvolumes