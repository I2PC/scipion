# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory Of Molecular Biology (MRC-LMB)
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


import re
from glob import glob

from pyworkflow.em import ProtClassify3D, ALIGN_PROJ, Float
from pyworkflow.em.data import SetOfVolumes, SetOfClassesVol
from pyworkflow.em.constants import ALIGN_NONE
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from xmipp3 import getEnviron
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from convert import (readSetOfClassesVol, getImageLocation,
                     writeSetOfVolumes, rowToAlignment)


MISSING_WEDGE_Y = 0
MISSING_WEDGE_X = 1
MISSING_PYRAMID = 2
MISSING_CONE = 3


class XmippProtMLTomo(ProtClassify3D):
    """ Align and classify 3D images with missing data regions in Fourier
    space, e.g. subtomograms or RCT reconstructions, by a 3D
    multi-reference refinement based on a maximum-likelihood (ML) target
    function.

    See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Ml_tomo_v31
    for further documentation
    """
    _label = 'mltomo'

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('results/mltomo_it%(iter)06d')
        myDict = {
                  'data_it': self.extraIter + '_img.xmd',
                  #'data': self._getExtraPath('results/mltomo_img.xmd'),
                  'classes_scipion': self.extraIter + '_classes_scipion.sqlite',
                  'log_it': self.extraIter + '_log.xmd',
                  'ref_it': self.extraIter + '_ref.xmd',
                  'ref': self._getExtraPath('results/mltomo_ref.xmd'),
                  'fsc_it': self.extraIter + '.fsc',
                  #'fsc': self._getExtraPath('results/mltomo.fsc'),
                  'volume': self.extraIter + '_ref%(ref3d)06d.vol'
                  }
        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('ref_it', iter=0).replace('000000', '??????')
        # Iterations will be identify by _itXXXXXX_ where XXXXXX is the
        # iteration number and is restricted to only 6 digits.
        self._iterRegex = re.compile('_it(\d{6,6})_')

    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='General params')
        form.addParam('inputVols', params.PointerParam,
                      pointerClass="SetOfVolumes",
                      label='Set of volumes',
                      help="Set of input volumes")
        form.addParam('copyAlignment', params.BooleanParam, default=True,
                      label='Consider previous alignment?',
                      help='If set to Yes, then alignment information '
                           'from input volumes will be considered.')
        form.addParam('generateRefs', params.BooleanParam, default=True,
                      label='Automatically generate references',
                      help="If set to true, 3D classes will be generated "
                           "automatically. Otherwise you can provide initial "
                           "reference volumes yourself.")
        form.addParam('numberOfReferences', params.IntParam,
                      label='Number of references', default=3,
                      condition="generateRefs",
                      help="Number of references to generate automatically")
        form.addParam('inputRefVols', params.PointerParam,
                      pointerClass="SetOfVolumes, Volume",
                      condition="not generateRefs",
                      label='Input reference volume(s)',
                      help="Provide a set of initial reference volumes")
        form.addParam('numberOfIterations', params.IntParam,
                      label='Number of iterations', default=25,
                      help="Maximum number of iterations to perform")
        form.addParam('symmetry', params.StringParam, default='c1',
                      label='Symmetry group',
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry "
                           "for a description of the symmetry groups format. If no "
                           "symmetry is present, give c1")
        form.addParam('missingDataType', params.EnumParam,
                      choices=['wedge_y', 'wedge_x', 'pyramid', 'cone'],
                      default=0,
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Missing data regions',
                      help="Provide missing data region type:\n\n"
                           "a) wedge_y for a missing wedge where the tilt "
                           "axis is along Y\n"
                           "b) wedge_x for a missing wedge where the tilt "
                           "axis is along X\n"
                           "c) pyramid for a missing pyramid where the tilt "
                           "axes are along Y and X\n"
                           "d) cone for a missing cone (pointing along Z)")
        form.addParam('missingAng', params.StringParam, default='-60 60',
                      label='Angles of missing data',
                      help='Provide angles for missing data area in the '
                           'following format:\n\n'
                           'for wedge_y or wedge_x: -60 60\n'
                           'for pyramid: -60 60 -60 60 (for y and x, respectively)\n'
                           'for cone: 45')
        form.addParam('maxCC', params.BooleanParam, default=False,
                      label='Use CC instead of ML',
                      help='Use constrained cross-correlation and weighted '
                           'averaging instead of ML')

        form.addSection(label='Sampling')
        form.addParam('angSampling', params.FloatParam, default=10.0,
                      label='Angular sampling (deg)',
                      help="Angular sampling rate (in degrees)")
        form.addParam('angSearch', params.FloatParam, default=-1.0,
                      label='Angular search range (deg)',
                      help="Angular search range around orientations of "
                           "input particles (by default [-1.0], exhaustive "
                           "searches are performed)")
        form.addParam('globalPsi', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Exhaustive psi search',
                      help="Exhaustive psi searches (only for c1 symmetry)")
        form.addParam('limitTrans', params.FloatParam, default=-1.0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Max shift (px)',
                      help="Maximum allowed shifts "
                           "(negative value means no restriction)")

        line = form.addLine('Tilt angle limits (deg)',
                            expertLevel=LEVEL_ADVANCED,
                            help='Limits for tilt angle search (in degrees)')
        line.addParam('tiltMin', params.FloatParam, default=0.0, label='Min')
        line.addParam('tiltMax', params.FloatParam, default=180.0, label='Max')
        form.addParam('psiSampling', params.FloatParam, default=-1.0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Psi angle sampling (deg)',
                      help="Angular sampling rate for the in-plane "
                           "rotations (in degrees)")

        form.addSection(label='Restrictions')
        form.addParam('dim', params.IntParam, default=-1,
                      label='Downscale input to (px)',
                      help="Use downscaled (in fourier space) images of this "
                           "size (in pixels)")
        form.addParam('maxRes', params.FloatParam, default=0.5,
                      label='Maximum resolution (px^-1)',
                      help="Maximum resolution (in pixel^-1) to use")
        form.addParam('doPerturb', params.BooleanParam, default=False,
                      label='Perturb',
                      help="Apply random perturbations to angular sampling "
                           "in each iteration")
        form.addParam('dontRotate', params.BooleanParam, default=False,
                      label='Do not rotate',
                      help="Keep orientations fixed, only translate and classify")
        form.addParam('dontAlign', params.BooleanParam, default=False,
                      label='Do not align',
                      help="Keep angles and shifts fixed "
                           "(otherwise start from random)")
        form.addParam('maskFile', params.PointerParam,
                      pointerClass='VolumeMask',
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Mask', condition='dontAlign',
                      help='Mask input volumes; only valid in combination '
                           'with --dont_align')
        form.addParam('onlyAvg', params.BooleanParam, default=False,
                      label='Only average',
                      help="Keep orientations and classes, "
                           "only output weighted averages")
        form.addParam('dontImpute', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Do not impute',
                      help='Use weighted averaging, rather than imputation')
        form.addParam('noImpThresh', params.FloatParam, default=1.0,
                      expertLevel=LEVEL_ADVANCED,
                      condition='dontImpute',
                      label='Threshold for averaging',
                      help='Threshold to avoid division by zero '
                           'for weighted averaging')
        form.addParam('fixSigmaNoise', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Fix sigma noise',
                      help='Do not re-estimate the standard deviation in '
                           'the pixel noise')
        form.addParam('fixSigmaOffset', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Fix sigma offsets',
                      help='Do not re-estimate the standard deviation in the '
                           'origin offsets')
        form.addParam('fixFrac', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Fix model fractions',
                      help='Do not re-estimate the model fractions. '
                           'Calculations start with even distribution.')

        form.addSection(label='Advanced')
        group = form.addGroup('Regularization parameters')
        line = group.addLine('Regularization',
                             help='Regularization parameters (in N/K^2)')

        line.addParam('regIni', params.FloatParam, default=0.0,
                      label='Initial')
        line.addParam('regFinal', params.FloatParam, default=0.0,
                      label='Final')
        group.addParam('regSteps', params.IntParam, default=5,
                       label='Steps',
                       help='Number of iterations in which the regularization '
                            'is changed from reg0 to regF')

        form.addParam('numberOfImpIterations', params.IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      label='Iterations in inner imputation loop', default=1,
                      help="Number of iterations for inner imputation loop")
        # FIXME: next param is to continue from iter X,
        #form.addParam('iterStart', params.IntParam,
         #             expertLevel=LEVEL_ADVANCED,
         #             label='Initial iteration', default=1,
         #             help="Number of initial iteration")
        form.addParam('eps', params.FloatParam, default=5e-5,
                      expertLevel=LEVEL_ADVANCED,
                      label='Stopping criterium',
                      help="Stopping criterium")
        form.addParam('stdNoise', params.FloatParam, default=1.0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Expected noise std',
                      help="Expected standard deviation for pixel noise")
        form.addParam('stdOrig', params.FloatParam, default=3.0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Expected origin offset std (px)',
                      help="Expected standard deviation for origin offset "
                           "(in pixels)")

        form.addParallelSection(threads=1, mpi=3)

    #--------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('convertInputs')
        self._insertFunctionStep('runMLTomo')
        self._insertFunctionStep('createOutput')
    
    #--------------------------- STEPS functions ------------------------------
    def convertInputs(self):
        inputVols = self.inputVols.get()
        self.createInputMd(inputVols)
        if not self.generateRefs:
            refVols = self.inputRefVols.get()
            self.createRefMd(refVols)

    def createInputMd(self, vols):
        fnVols = self._getExtraPath('input_volumes.xmd')
        if self.copyAlignment:
            alignType = vols.getAlignment()
        else:
            alignType = ALIGN_NONE

        writeSetOfVolumes(vols, fnVols,
                          postprocessImageRow=self._postprocessVolumeRow,
                          alignType=alignType)
        if not vols.hasAlignment() or not self.copyAlignment:
            mdFn = md.MetaData(fnVols)
            mdFn.fillConstant(md.MDL_ANGLE_ROT, 0.)
            mdFn.fillConstant(md.MDL_ANGLE_TILT, 0.)
            mdFn.fillConstant(md.MDL_ANGLE_PSI, 0.)
            mdFn.fillConstant(md.MDL_SHIFT_X, 0.)
            mdFn.fillConstant(md.MDL_SHIFT_Y, 0.)
            mdFn.fillConstant(md.MDL_SHIFT_Z, 0.)
            mdFn.write(fnVols, md.MD_OVERWRITE)

        # set missing angles
        missType = ['wedge_y', 'wedge_x', 'pyramid', 'cone']
        missNum = self.missingDataType.get()
        missAng = self.missingAng.get()
        missDataFn = self._getExtraPath('wedges.xmd')
        missAngValues = str(missAng).strip().split()
        thetaY0, thetaYF, thetaX0, thetaXF = 0, 0, 0, 0
        mdFn = md.MetaData()

        if missNum == MISSING_WEDGE_X:
            thetaX0, thetaXF = missAngValues
        elif missNum == MISSING_WEDGE_Y:
            thetaY0, thetaYF = missAngValues
        elif missNum == MISSING_PYRAMID:
            thetaY0, thetaYF, thetaX0, thetaXF = missAngValues
        else:  # MISSING_CONE
            thetaY0 = missAngValues[0]

        for i in range(1, vols.getSize() + 1):
            row = md.Row()
            row.setValue(md.MDL_MISSINGREGION_NR, missNum + 1)
            row.setValue(md.MDL_MISSINGREGION_TYPE, missType[missNum])
            row.setValue(md.MDL_MISSINGREGION_THX0, float(thetaX0))
            row.setValue(md.MDL_MISSINGREGION_THXF, float(thetaXF))
            row.setValue(md.MDL_MISSINGREGION_THY0, float(thetaY0))
            row.setValue(md.MDL_MISSINGREGION_THYF, float(thetaYF))
            row.addToMd(mdFn)

        mdFn.write(missDataFn, md.MD_APPEND)

    def createRefMd(self, vols):
        refVols = self._getExtraPath('ref_volumes.xmd')
        mdFn = md.MetaData()

        if self.isSetOfVolumes():
            for vol in vols:
                imgId = vol.getObjId()
                row = md.Row()
                row.setValue(md.MDL_ITEM_ID, long(imgId))
                row.setValue(md.MDL_IMAGE, getImageLocation(vol))
                row.setValue(md.MDL_ENABLED, 1)
                row.addToMd(mdFn)
        else:
            imgId = vols.getObjId()
            row = md.Row()
            row.setValue(md.MDL_ITEM_ID, long(imgId))
            row.setValue(md.MDL_IMAGE, getImageLocation(vols))
            row.setValue(md.MDL_ENABLED, 1)
            row.addToMd(mdFn)

        mdFn.write(refVols, md.MD_APPEND)

    def runMLTomo(self):
        fnVols = self._getExtraPath('input_volumes.xmd')
        refVols = self._getExtraPath('ref_volumes.xmd')
        missDataFn = self._getExtraPath('wedges.xmd')
        outDir = self._getExtraPath("results")
        pwutils.makePath(outDir)

        params = ' -i %s' % fnVols
        params += ' --oroot %s' % (outDir + '/mltomo')
        params += ' --iter %d' % self.numberOfIterations.get()
        params += ' --sym %s' % self.symmetry.get()
        params += ' --missing %s' % missDataFn
        params += ' --maxres %0.2f' % self.maxRes.get()
        params += ' --dim %d' % self.dim.get()
        params += ' --ang %0.1f' % self.angSampling.get()
        params += ' --ang_search %0.1f' % self.angSearch.get()
        params += ' --limit_trans %0.1f' % self.limitTrans.get()
        params += ' --tilt0 %0.1f --tiltF %0.1f' % (self.tiltMin.get(),
                                                    self.tiltMax.get())
        params += ' --psi_sampling %0.1f' % self.psiSampling.get()
        params += ' --reg0 %0.1f --regF %0.1f --reg_steps %d' % \
                  (self.regIni.get(), self.regFinal.get(), self.regSteps.get())
        params += ' --impute_iter %d' % self.numberOfImpIterations.get()
        #params += ' --istart %d' % self.iterStart.get()
        params += ' --eps %0.2f' % self.eps.get()
        params += ' --pixel_size %0.2f' % self.inputVols.get().getSamplingRate()
        params += ' --noise %0.1f --offset %0.1f' % (self.stdNoise.get(),
                                                     self.stdOrig.get())
        params += ' --thr %d' % self.numberOfThreads.get()

        if self.generateRefs:
            params += ' --nref %d' % self.numberOfReferences.get()
        else:
            params += ' --ref %s' % refVols
        if self.copyAlignment:
            params += ' --keep_angles'
            # when considering alignment, use it,
            # otherwise starts from random angles
        if self.doPerturb:
            params += ' --perturb'
        if self.dontRotate:
            params += ' --dont_rotate'
        if self.dontAlign:
            params += ' --dont_align'
        if self.onlyAvg:
            params += ' --only_average'
        if self.dontImpute:
            params += ' --dont_impute'
            params += ' --noimp_threshold %0.1f' % self.noImpThresh.get()
        if self.fixSigmaNoise:
            params += ' --fix_sigma_noise'
        if self.fixSigmaOffset:
            params += ' --fix_sigma_offset'
        if self.fixFrac:
            params += ' --fix_fractions'
        if self.globalPsi:
            params += ' --dont_limit_psirange'
        if self.maxCC:
            params += ' --maxCC'
        if self.dontAlign and self.maskFile:
            params += ' --mask %s' % self.maskFile.get().getFileName()

        self.runJob('xmipp_ml_tomo', '%s' % params,
                    env=self.getMLTomoEnviron(),
                    numberOfMpi=self.numberOfMpi.get(),
                    numberOfThreads=self.numberOfThreads.get())
    
    def createOutput(self):
        # output files:
        #   mltomo_ref.xmd contains all info for output 3D classes
        #   mltomo_refXXXXXX.vol output volume - 3D class
        #   mltomo_img.xmd contains alignment metadata for all vols
        #   mltomo.fsc
        outputGlobalMdFn = self._getFileName('ref')
        setOfClasses = self._createSetOfClassesVol()
        setOfClasses.setImages(self.inputVols.get())
        readSetOfClassesVol(setOfClasses, outputGlobalMdFn)
        self._defineOutputs(outputClasses=setOfClasses)
        self._defineSourceRelation(self.inputVols, self.outputClasses)

    #--------------------------- INFO functions -------------------------------
    def _summary(self):
        messages = []
        if not hasattr(self, 'outputClasses'):
            messages.append('Output is not ready')
        else:
            messages.append('Number of input volumes: %d' %
                            self.inputVols.get().getSize())
            if self.generateRefs:
                messages.append('References were auto-generated')
            else:
                messages.append('References were provided by user')
            messages.append('Number of output classes: %d' %
                            self.outputClasses.getSize())

        return messages

    def _validate(self):
        errors = []
        missNum = self.missingDataType.get()
        angString = self.missingAng.get()
        angs = str(angString).split()
        if (missNum == 0 or missNum == 1) and len(angs) != 2:
            errors.append('Wrong angles of missing data! Provide two values '
                          'for a missing wedge')
        elif missNum == 2 and len(angs) != 4:
            errors.append('Wrong angles of missing data! Provide four values '
                          'for a missing pyramid')
        elif missNum == 3 and len(angs) != 1:
            errors.append('Wrong angles of missing data! Provide one value '
                          'for a missing cone')

        if not self.copyAlignment and self.angSearch != -1:
            errors.append('You cannot do local searches when '
                          'ignoring input alignments.')

        return errors

    def _citations(self):
        return ['Scheres2009c']

    #--------------------------- UTILS functions ------------------------------
    def getMLTomoEnviron(self):
        env = getEnviron()
        return env

    def isSetOfVolumes(self):
        return isinstance(self.inputRefVols.get(), SetOfVolumes)

    def _postprocessVolumeRow(self, img, imgRow):
        # explicitly set this from protocol input
        # to avoid conflict with input metadata
        missNum = self.missingDataType.get()
        imgRow.setValue(md.MDL_MISSINGREGION_NR, missNum + 1)

        if not imgRow.hasLabel(md.MDL_REF):
            imgRow.setValue(md.MDL_REF, 1)
        if not imgRow.hasLabel(md.MDL_LL):
            imgRow.setValue(md.MDL_LL, 1.)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 6 digits iteration number
        return result

    def _getIterClasses(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration mltomo_it??????_img.xmd file.
        """
        data_classes = self._getFileName('classes_scipion', iter=it)

        if clean:
            pwutils.cleanPath(data_classes)

        if not pwutils.exists(data_classes):
            clsSet = SetOfClassesVol(filename=data_classes)
            clsSet.setImages(self.inputVols.get())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes

    def _loadClassesInfo(self, iteration):
        """ Read some information about the produced MLtomo 3D classes
        from the *ref.xmd file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
        modelFn = md.MetaData('classes@' +
                              self._getFileName('ref_it', iter=iteration))

        for classNumber, row in enumerate(md.iterRows(modelFn)):
            fn = str(row.getValue(md.MDL_IMAGE))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            self._classesInfo[classNumber + 1] = (fn, row.clone())

    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClassesVol from a given iteration. """
        self._loadClassesInfo(iteration)
        dataFn = self._getFileName('data_it', iter=iteration)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows(dataFn,
                                                          sortByLabel=md.MDL_ITEM_ID))

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.MDL_REF))
        item.setTransform(rowToAlignment(row, ALIGN_PROJ))
        item._xmippLogLikeliContribution = Float(row.getValue(md.MDL_LL))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            fn, row = self._classesInfo[classId]
            item.setAlignmentProj()
            item.getRepresentative().setLocation(1, fn)
            item._xmippWeight = Float(row.getValue(md.MDL_WEIGHT))
            item._xmippSignalChange = Float(row.getValue(md.MDL_SIGNALCHANGE))
