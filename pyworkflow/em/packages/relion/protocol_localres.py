# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

from pyworkflow.protocol.params import (PointerParam, FloatParam, FileParam,
                                        IntParam, LabelParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D, ImageHandler
from convert import isVersion2
from pyworkflow.utils import exists
import pyworkflow.utils.path as putils


class ProtRelionLocalRes(ProtAnalysis3D):
    """
    Relion local-resolution estimation protocol.

    This program basically performs a series of post-processing operations
    with a small soft, spherical mask that is moved over the entire map,
    while using phase-randomisation to estimate the convolution effects
    of that mask.
    """
    _label = 'local resolution'

    @classmethod
    def isDisabled(cls):
        return not isVersion2()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
                 'half1': self._getTmpPath("relion_half1_class001_unfil.mrc"),
                 'half2': self._getTmpPath("relion_half2_class001_unfil.mrc"),
                 'finalMap' : self._getExtraPath('relion_locres_filtered.mrc'),
                 'resolMap' : self._getExtraPath('relion_locres.mrc')
                 }

        self._updateFilenamesDict(myDict)

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        
        form.addSection(label='Input')
        form.addParam('protRefine', PointerParam,
                      pointerClass="ProtRefine3D",
                      label='Select a previous refinement protocol',
                      help='Select any previous refinement protocol to get the '
                           '3D half maps. Note that it is recommended that the '
                           'refinement protocol uses a gold-standard method.')
        form.addParam('mtf', FileParam,
                       label='MTF-curve file',
                       help='User-provided STAR-file with the MTF-curve '
                            'of the detector.'
                            'Relion param: <--mtf>')
        form.addParam('bfactor', FloatParam, default=-250,
                       label='Provide B-factor:',
                       help='Probably, the overall B-factor as was '
                            'estimated in the postprocess is a useful '
                            'value for here. Use negative values for '
                            'sharpening. Be careful: if you over-sharpen '
                            'your map, you may end up interpreting '
                            'noise for signal!')

        form.addSection(label='LocalRes')
        form.addParam('Msg', LabelParam,
                      label='Select Advanced level if you want to adjust the '
                            'parameters')
        form.addParam('locResSamp', IntParam, default=25,
                      label='Sampling rate (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Sampling rate (in Angstroms) with which to '
                           'sample the local-resolution map')
        form.addParam('locResMaskRad', IntParam, default=-1,
                      label='Mask radius (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Radius (in A) of spherical mask for '
                           'local-resolution map (default = 0.5*sampling)')
        form.addParam('locResEdgeWidth', IntParam, default=-1,
                      label='Edge width (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Width of soft edge (in A) on masks for '
                           'local-resolution map (default = sampling)')
        form.addParam('locResRand', FloatParam, default=25.0,
                      label='Randomize phases from (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Randomize phases from this resolution (in A)')
        form.addParam('locResMin', IntParam, default=50,
                      label='Lowest res limit (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Lowest local resolution allowed (in A)')

        form.addParallelSection(threads=0, mpi=1)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        objId = self.protRefine.get().getObjId()
        self._createFilenameTemplates()
        self._defineParamDict()
        self._insertFunctionStep('initializeStep', objId)
        self._insertFunctionStep('postProcessStep', self.paramDict)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def initializeStep(self, protId):
        protRef = self.protRefine.get()
        protClassName = protRef.getClassName()
        if protClassName.startswith('ProtRelionRefine3D'):
            protRef._createFilenameTemplates()
            half1Map = protRef._getFileName("final_half1_volume", ref3d=1)
            half2Map = protRef._getFileName("final_half2_volume", ref3d=1)
            
            putils.copyFile(self._getRelionMapFn(half1Map), self._getTmpPath())
            putils.copyFile(self._getRelionMapFn(half2Map), self._getTmpPath())
            
        elif protClassName.startswith('ProtFrealign'):
            protRef._createFilenameTemplates()
            lastIter = protRef._getLastIter()

            half1Map = protRef._getFileName('iter_vol1', iter=lastIter)
            half2Map = protRef._getFileName('iter_vol2', iter=lastIter)
            
            putils.copyFile(half1Map, self._getFileName("half1"))
            putils.copyFile(half2Map, self._getFileName("half2"))

        elif protClassName.startswith('XmippProtProjMatch'):
            iterN = protRef.getLastIter()
            protRef._initialize()
            half1 = protRef._getFileName('reconstructedFileNamesItersSplit1',
                                         iter=iterN, ref=1)
            half2 = protRef._getFileName('reconstructedFileNamesItersSplit2',
                                         iter=iterN, ref=1)
            ih = ImageHandler()
            ih.convert(half1, self._getFileName("half1"))
            ih.convert(half2, self._getFileName("half2"))
            
        elif protClassName.startswith('EmanProtRefine'):
            protRef._createFilenameTemplates()
            numRun = protRef._getRun()
            protRef._createIterTemplates(numRun)
            iterN = protRef._lastIter()
            print "RUN ITER: ", numRun, iterN, protRef._iterTemplate
            half1 = protRef._getFileName("mapEvenUnmasked", run=numRun)
            half2 = protRef._getFileName("mapOddUnmasked", run=numRun)
            
            ih = ImageHandler()
            ih.convert(half1, self._getFileName("half1"))
            ih.convert(half2, self._getFileName("half2"))
    
    def postProcessStep(self, paramDict):
        params = ' '.join(['%s %s' % (k, str(v))
                           for k, v in self.paramDict.iteritems()])
        program = 'relion_postprocess'
        if self.numberOfMpi.get() > 1:
            program += '_mpi'

        self.runJob(program, params)
    
    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getFileName("finalMap"))
        vol = self.protRefine.get().outputVolume
        pxSize = vol.getSamplingRate()
        volume.setSamplingRate(pxSize)
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(vol, volume)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        errors = []
        mtfFile = self.mtf.get()
        
        if mtfFile and not exists(mtfFile):
            errors.append("Missing MTF-file '%s'" % mtfFile)

        protClassName = self.protRefine.get().getClassName()
        if protClassName.startswith('SpiderProtRefinement'):
            errors.append("Relion local resolution protocol is not "
                          "implemented for Spider - refinement.")
        
        if protClassName.startswith('XmippProtReconstructHighRes'):
            errors.append("Relion local resolution protocol is not "
                          "implemented for Xmipp - highres.")
        return errors
    
    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume not ready yet.")
        else:
            output = self.outputVolume
            summary.append("%s: Output volume was locally filtered "
                           "and sharpened" % self.getObjectTag(output))
        return summary
    
    #--------------------------- UTILS functions -------------------------------
    def _defineParamDict(self):
        """ Define all parameters to run relion_postprocess"""
        volume = self.protRefine.get().outputVolume
        self.paramDict = {'--i': self._getTmpPath("relion"),
                          '--o': self._getExtraPath('relion'),
                          '--angpix': volume.getSamplingRate(),
                          '--adhoc_bfac': self.bfactor.get(),
                          '--locres': '',
                          # Expert options
                          '--locres_sampling': self.locResSamp.get(),
                          '--locres_maskrad': self.locResMaskRad.get(),
                          '--locres_edgwidth': self.locResEdgeWidth.get(),
                          '--locres_randomize_at': self.locResRand.get(),
                          '--locres_minres': self.locResMin.get()
                          }

        mtfFile = self.mtf.get()
        if mtfFile:
            self.paramDict['--mtf'] = mtfFile

    def _getRelionMapFn(self, fn):
        return fn.split(':')[0]
