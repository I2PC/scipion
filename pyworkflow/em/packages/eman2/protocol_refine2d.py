# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

import os
import re
from os.path import exists
from glob import glob

import pyworkflow.em as em
from pyworkflow import VERSION_1_2
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        EnumParam, StringParam,
                                        BooleanParam, LabelParam)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import makePath, cleanPath


from convert import (rowToAlignment, createEmanProcess,
                     writeSetOfParticles, convertReferences)
from constants import *
from eman2 import getEmanProgram, validateVersion, isNewVersion, SCRATCHDIR



class EmanProtRefine2D(em.ProtClassify2D):
    """
    This Protocol wraps *e2refine2d.py* Eman2 program.

    This program is used to produce reference-free class averages
    from a population of mixed, unaligned particle images. These averages
    can be used to generate initial models or assess the structural
    variability of the data. They are not normally themselves used as part
    of the single particle reconstruction refinement process, which
    uses the raw particles in a reference-based classification
    approach. However, with a good structure, projections of the
    final 3-D model should be consistent with the results of
    this reference-free analysis.

    This program uses a fully automated iterative alignment/MSA approach.
    You should normally target a minimum of 10-20 particles per
    class-average, though more is fine.

    Default parameters should give a good start, but are likely not
    optimal for any given system.

    Note that it does have the --parallel option, but a few steps of the
    iterative process are not parallellised, so don't be surprised if
    multiple cores are not always active.
"""
    _label = 'refine 2D'
    _lastUpdateVersion = VERSION_1_2

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {
            'partFlipSet': 'sets/inputSet__ctf_flip.lst',
            'initialAvgSet': self._getExtraPath('initial_averages.hdf'),
            'classes_scipion': self._getExtraPath('classes_scipion_it%(iter)02d.sqlite'),
            'classes': 'r2d_01/classes_%(iter)02d.hdf',
            'cls': 'r2d_01/classmx_%(iter)02d.hdf',
            'results': self._getExtraPath('results_it%(iter)02d.txt'),
            'allrefs': self._getExtraPath('r2d_01/allrefs_%(iter)02d.hdf'),
            'alirefs': self._getExtraPath('r2d_01/aliref_%(iter)02d.hdf'),
            'basis': self._getExtraPath('r2d_01/basis_%(iter)02d.hdf')
        }
        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        clsFn = self._getExtraPath(self._getFileName('classes', iter=1))
        self._iterTemplate = clsFn.replace('classes_01', 'classes_??')
        # Iterations will be identify by classes_XX_ where XX is the iteration
        #  number and is restricted to only 2 digits.
        self._iterRegex = re.compile('classes_(\d{2,2})')

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      label="Input particles",
                      important=True, pointerClass='SetOfParticles',
                      help='Select the input particles.')
        form.addParam('inputClassAvg', PointerParam,
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label="Input class averages",
                      pointerClass='SetOfAverages',
                      help='Select starting class averages. If not provided, '
                           'they will be generated automatically.')
        form.addParam('numberOfClassAvg', IntParam, default=32,
                      label='Number of class-averages',
                      help='Number of class-averages to generate. Normally you '
                           'would want a minimum of ~10-20 particles per class on '
                           'average, but it is fine to have 100-200 for a large data '
                           'set. If you plan on making a large number (>100) of '
                           'classes, you should use the *Fast seed* option. Note '
                           'that these averages are not used for final 3-D '
                           'refinement, so generating a very large number is not '
                           'useful in most situations.')
        form.addParam('numberOfIterations', IntParam, default=3,
                      label='Number of iterations',
                      help='Number of iterations of the overall 2-D refinement '
                           'process to run. For high contrast data, 4-5 iterations '
                           'may be more than enough, but for low contrast data '
                           'it could take 10-12 iterations to converge well.')
        form.addParam('nbasisfp', IntParam, default=8,
                      label='Number of MSA vectors to use',
                      help='Number of MSa basis vectors to use when '
                           'classifying particles.')
        form.addParam('numberOfAlignRef', IntParam, default=5,
                      label='Number of alignment references',
                      help='The number of alignment references to use in each '
                           'iteration. You can look at this as the number of '
                           'different highly distinct views your particle has. '
                           'With something like GroEL with mostly side views '
                           'and top views, 3-4 is sufficient. With something like '
                           'a ribosome something more like 10-15 would be '
                           'appropriate.')
        form.addParam('doNormProj', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Normalize projected vectors?',
                      help='Normalizes each projected vector into the MSA '
                           'subspace. Note that this is different from normalizing '
                           'the input images since the subspace is not expected to '
                           'fully span the image')
        form.addParam('doFastSeed', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Use fast seed?',
                      help='Will seed the k-means loop quickly, but may '
                           'produce less consistent results. Always use this '
                           'when generating >~ 100 classes.')
        if isNewVersion():
            form.addParam('doAutomask', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label='Automask class-averages?',
                          help='This will perform a 2-D automask on class-averages '
                               'to help with centering. May be useful for negative '
                               'stain data particularly.')

        line = form.addLine('Centering: ',
                            help="If the default centering algorithm "
                                 "(xform.center) doesn't work well, "
                                 "you can specify one of the others "
                                 "here (e2help.py processor center)")
        line.addParam('centerType', EnumParam,
                      choices=['nocenter', 'xform.center',
                               'xform.centeracf', 'xform.centerofmass', 'None'],
                      label="type", default=XFORM_CENTER,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('centerParams', StringParam, default='',
                      label='params')

        form.addSection(label='Similarity matrix')
        form.addParam('paramsMsg', LabelParam, default=True,
                      label='These parameters are for advanced users only!\n',
                      help='For help please address to EMAN2 %s or run:\n'
                            '*scipion run e2help.py cmp -v 2* or\n'
                            '*scipion run e2help.py aligners -v 2*' % WIKI_URL)
        form.addParam('shrink', IntParam, default=1,
                      label='Shrink particles',
                      help='Optionally shrink the input particles by an integer '
                           'amount prior to computing similarity scores. '
                           'For speed purposes.')
        line = form.addLine('simcmp: ',
                            help='The name of a cmp to be used in comparing '
                                 'the aligned images (default=ccc)')
        line.addParam('simcmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simcmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('First stage aligner')
        line = group.addLine('simalign: ')
        line.addParam('simalignType', EnumParam,
                      choices=['frm2d', 'rotate_flip',
                               'rotate_flip_iterative', 'rotate_precenter',
                               'rotate_trans_flip_scale',
                               'rotate_trans_flip_scale_iter',
                               'rotate_trans_scale_iter',
                               'rotate_translate', 'rotate_translate_flip',
                               'rotate_translate_flip_iterative',
                               'rotate_translate_flip_resample',
                               'rotate_translate_iterative',
                               'rotate_translate_resample',
                               'rotate_translate_scale', 'rotate_translate_tree',
                               'rotational', 'rotational_iterative', 'rtf_exhaustive',
                               'rtf_slow_exhaustive', 'scale', 'symalign',
                               'symalignquat', 'translational', 'None'],
                      label='type', default=ALN_ROTATE_TRANSLATE_TREE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simalignParams', StringParam,
                      default='', label='params')
        line = group.addLine('simaligncmp: ')
        line.addParam('simaligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simaligncmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('Second stage aligner')
        line = group.addLine('simralign: ')
        line.addParam('simralignType', EnumParam,
                      choices=['None', 'refine'],
                      label='type', default=RALN_NONE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simralignParams', StringParam,
                      default='', label='params')
        line = group.addLine('simraligncmp: ')
        line.addParam('simraligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_DOT,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simraligncmpParams', StringParam,
                      default='', label='params')

        form.addSection(label='Class averaging')
        form.addParam('paramsMsg2', LabelParam, default=True,
                      label='These parameters are for advanced users only!\n',
                      help='For help please address to EMAN2 %s or run:\n'
                           '*scipion run e2help.py cmp -v 2* or\n'
                           '*scipion run e2help.py averagers -v 2*' % WIKI_URL)
        form.addParam('classIter', IntParam, default=4,
                      label='Number of iterations for class-averages',
                      help='Number of iterations to use when making '
                           'class-averages (default=5)')
        form.addParam('classKeep', FloatParam, default=0.8,
                      label='Fraction of particles to keep',
                      help='The fraction of particles to keep in each class, '
                      'based on the similarity score generated by cmp '
                      '(default=0.8)')
        form.addParam('classKeepSig', BooleanParam, default=False,
                      label='Keep particles based on sigma?',
                      help='Change the *keep* criterion from fraction-based '
                      'to sigma-based')
        form.addParam('classRefSf', BooleanParam, default=False,
                      label='Filter class-averages?',
                      expertLevel=LEVEL_ADVANCED,
                      help='Use setsfref option in class averaging to '
                           'produce better class averages')

        line = form.addLine('classnormproc: ',
                            help='Normalization applied during class-averaging')
        line.addParam('classnormprocType', EnumParam,
                      choices=['normalize', 'normalize.bymass',
                               'normalize.circlemean', 'normalize.edgemean',
                               'normalize.local', 'normalize.lredge',
                               'normalize.mask', 'normalize.maxmin',
                               'normalize.ramp.normvar', 'normalize.rows',
                               'normalize.toimage', 'normalize.unitlen',
                               'normalize.unitsum', 'None'],
                      label='type',
                      default=PROC_NORMALIZE_EDGEMEAN,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classnormprocParams', StringParam,
                      default='', label='params')

        form.addParam('classAveragerType', EnumParam,
                      choices=['absmaxmin', 'ctf.auto', 'ctf.weight',
                               'ctf.weight.autofilt', 'ctfw.auto', 'iteration',
                               'localweight', 'mean', 'mean.tomo',
                               'minmax', 'sigma', 'weightedfourier'],
                      label='Class averager type: ',
                      default=AVG_CTF_WEIGHT_AUTOFILT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='The averager used to generated class-averages')

        line = form.addLine('classnormproc: ',
                            help='Normalization applied during class-averaging')
        line.addParam('classnormprocType', EnumParam,
                      choices=['normalize', 'normalize.bymass',
                               'normalize.circlemean', 'normalize.edgemean',
                               'normalize.local', 'normalize.lredge',
                               'normalize.mask', 'normalize.maxmin',
                               'normalize.ramp.normvar', 'normalize.rows',
                               'normalize.toimage', 'normalize.unitlen',
                               'normalize.unitsum', 'None'],
                      label='type',
                      default=PROC_NORMALIZE_EDGEMEAN,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classnormprocParams', StringParam,
                      default='', label='params')

        line = form.addLine('classcmp: ')
        line.addParam('classcmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classcmpParams', StringParam,
                      default='', label='params',
                      help='The name of a cmp to be used in class averaging '
                           '(default=ccc)')

        group = form.addGroup('First stage aligner (clsavg)')
        line = group.addLine('classalign: ')
        line.addParam('classalignType', EnumParam,
                      choices=['frm2d', 'rotate_flip',
                               'rotate_flip_iterative', 'rotate_precenter',
                               'rotate_trans_flip_scale',
                               'rotate_trans_flip_scale_iter',
                               'rotate_trans_scale_iter',
                               'rotate_translate', 'rotate_translate_flip',
                               'rotate_translate_flip_iterative',
                               'rotate_translate_flip_resample',
                               'rotate_translate_iterative',
                               'rotate_translate_resample',
                               'rotate_translate_scale', 'rotate_translate_tree',
                               'rotational', 'rotational_iterative', 'rtf_exhaustive',
                               'rtf_slow_exhaustive', 'scale', 'symalign',
                               'symalignquat', 'translational', 'None'],
                      label='type', default=ALN_ROTATE_TRANSLATE_TREE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classalignParams', StringParam,
                      default='flip=0', label='params')
        line = group.addLine('classaligncmp: ')
        line.addParam('classaligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classaligncmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('Second stage aligner (clsavg)')
        line = group.addLine('classralign: ')
        line.addParam('classralignType', EnumParam,
                      choices=['None', 'refine'],
                      label='type', default=RALN_NONE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classralignParams', StringParam,
                      default='', label='params')
        line = group.addLine('classraligncmp: ')
        line.addParam('classraligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classraligncmpParams', StringParam,
                      default='', label='params')

        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._createIterTemplates()
        self._insertFunctionStep('convertImagesStep')
        args = self._prepareParams()
        self._insertFunctionStep('refineStep', args)
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------
    def convertImagesStep(self):
        partSet = self._getInputParticles()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)

        if partSet.hasCTF():
            program = getEmanProgram('e2ctf.py')
            acq = partSet.getAcquisition()

            args = " --voltage %d" % acq.getVoltage()
            args += " --cs %f" % acq.getSphericalAberration()
            args += " --ac %f" % (100 * acq.getAmplitudeContrast())
            args += " --threads=%d" % self.numberOfThreads.get()
            if not partSet.isPhaseFlipped():
                args += " --phaseflip"
            args += " --computesf --apix %f" % partSet.getSamplingRate()
            args += " --allparticles --autofit --curdefocusfix --storeparm -v 8"
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

        program = getEmanProgram('e2buildsets.py')
        args = " --setname=inputSet__ctf_flip --allparticles --minhisnr=-1"
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

        if self.inputClassAvg.hasValue():
            avgs = self.inputClassAvg.get()
            outputFn = self._getFileName('initialAvgSet')
            convertReferences(avgs, outputFn)

    def refineStep(self, args):
        """ Run the EMAN program to refine 2d. """
        program = getEmanProgram('e2refine2d.py')
        # mpi and threads are handled by EMAN itself
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        partSet = self._getInputParticles()
        classes2D = self._createSetOfClasses2D(partSet)
        self._fillClassesFromIter(classes2D, self._lastIter())

        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(self.inputParticles, classes2D)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        validateVersion(self, errors)

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Input Particles: %s" % self.getObjectTag('inputParticles'))
            summary.append("Classified into *%d* classes." % self.numberOfClassAvg)
            summary.append("Output set: %s" % self.getObjectTag('outputClasses'))

        summary.append('\n\n*Note:* final class averages produced by EMAN are '
                       'not aligned, while the particle inside each class are.')
        return summary

    def _methods(self):
        methods = "We classified input particles %s (%d items) " % (
            self.getObjectTag('inputParticles'),
            self._getInputParticles().getSize())
        methods += "into %d classes using e2refine2d.py " % self.numberOfClassAvg
        return [methods]

    #--------------------------- UTILS functions --------------------------------------------
    def _prepareParams(self):
        args1 = " --input=%s" % self._getParticlesStack()
        if self.inputClassAvg.hasValue():
            args1 += " --initial=%s" % self._getBaseName('initialAvgSet')
        args2 = self._commonParams()
        args = args1 + args2

        return args

    def _commonParams(self):
        args = " --ncls=%(ncls)d --iter=%(numberOfIterations)d --nbasisfp=%(nbasisfp)d"
        args += " --naliref=%(naliref)d"
        args += " --classkeep=%(classKeep)f --classiter=%(classiter)d "
        args += " --classaverager=%s" % self.getEnumText('classAveragerType')

        if self.doNormProj:
            args += " --normproj"
        if self.doFastSeed:
            args += " --fastseed"
        if self.shrink > 1:
            args += " --shrink %d" % self.shrink.get()
        if self.classKeepSig:
            args += " --classkeepsig"
        if self.classRefSf:
            args += " --classrefsf"

        if isNewVersion() and self.doAutomask:
            args += " --automask"

        for param in ['simcmp', 'simalign', 'simralign', 'classnormproc',
                      'classcmp', 'classalign', 'center',
                      'classaligncmp', 'classralign', 'classraligncmp']:
            args += self._getOptsString(param)

        if self.numberOfMpi > 1:
            args += " --parallel=mpi:%(mpis)d:%(scratch)s"
        else:
            args += " --parallel=thread:%(threads)d"

        params = {'ncls': self.numberOfClassAvg.get(),
                  'numberOfIterations': self.numberOfIterations.get(),
                  'nbasisfp': self.nbasisfp.get(),
                  'naliref': self.numberOfAlignRef.get(),
                  'classKeep': self.classKeep.get(),
                  'classiter': self.classIter.get(),
                  'threads': self.numberOfThreads.get(),
                  'mpis': self.numberOfMpi.get(),
                  'scratch': SCRATCHDIR
                  }
        args = args % params

        return args

    def _getBaseName(self, key, **args):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key, **args))

    def _getParticlesStack(self):
            return self._getFileName("partFlipSet")

    def _iterTextFile(self, iterN):
        f = open(self._getFileName('results', iter=iterN))

        for line in f:
            if '#' not in line:
                yield map(float, line.split())

        f.close()

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 2 digits iteration number

        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _getIterClasses(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration data.star file.
        """
        data_classes = self._getFileName('classes_scipion', iter=it)

        if clean:
            cleanPath(data_classes)

        if not exists(data_classes):
            clsSet = em.SetOfClasses2D(filename=data_classes)
            clsSet.setImages(self.inputParticles.get())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _fillClassesFromIter(self, clsSet, iterN):
        self._execEmanProcess(iterN)
        params = {'orderBy' : ['_micId', 'id'],
                  'direction' : 'ASC'
                  }
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=self._iterTextFile(iterN),
                             iterParams=params)

    def _execEmanProcess(self, iterN):
        clsFn = self._getFileName("cls", iter=iterN)
        classesFn = self._getFileName("classes", iter=iterN)

        proc = createEmanProcess(args='read %s %s %s %s 2d'
                                 % (self._getParticlesStack(), clsFn, classesFn,
                                    self._getBaseName('results', iter=iterN)),
                                 direc=self._getExtraPath())
        proc.wait()

        self._classesInfo = {}  # store classes info, indexed by class id
        for classId in range(self.numberOfClassAvg.get()):
            self._classesInfo[classId + 1] = (classId + 1,
                                              self._getExtraPath(classesFn))

    def _getOptsString(self, option):
        optionType = "optionType = self.getEnumText('" + option + "Type')"
        optionParams = 'optionParams = self.' + option + 'Params.get()'
        exec(optionType)
        exec(optionParams)

        if optionType == 'None':
            return ''
        if optionParams != '':
            argStr = ' --%s=%s:%s' % (option, optionType, optionParams)
        else:
            argStr = ' --%s=%s' % (option, optionType)

        return argStr

    def _updateParticle(self, item, row):
        if row[1] == 1:  # enabled
            item.setClassId(row[2] + 1)
            item.setTransform(rowToAlignment(row[3:], em.ALIGN_2D))
        else:
            setattr(item, "_appendItem", False)

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn = self._classesInfo[classId]
            item.setAlignment2D()
            item.getRepresentative().setLocation(classId, fn)
