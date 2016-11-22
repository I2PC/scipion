# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (jgomez@cnb.csic.es)
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

from os.path import basename, abspath

from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


# Some utility functions to import movies that are used in several tests.
class TestXmippBase(BaseTest):
    @classmethod
    def setData(cls):
        cls.dataset = DataSet.getDataSet('movies')
        cls.movie1 = cls.dataset.getFile('qbeta/qbeta.mrc')
        cls.movie2 = cls.dataset.getFile('cct/cct_1.em')
    
    @classmethod
    def runImportMovie(cls, pattern, samplingRate, voltage, scannedPixelSize,
                       magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """

        kwargs = {
            'filesPath': pattern,
            'magnification': magnification,
            'voltage': voltage,
            'sphericalAberration': sphericalAberration
        }

        # We have two options: pass the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        if samplingRate is not None:
            kwargs.update({'samplingRateMode': 0,
                           'samplingRate': samplingRate})
        else:
            kwargs.update({'samplingRateMode': 1,
                           'scannedPixelSize': scannedPixelSize})

        cls.protImport = cls.newProtocol(ProtImportMovies, **kwargs)
        cls.proj.launchProtocol(cls.protImport, wait=True)

        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ",
                            cls.protImport.getErrorMessage())

        # Check that input movies have been imported (a better way to do this?)
        if cls.protImport.outputMovies is None:
            raise Exception('Import of movies: %s, failed, '
                            'outputMovies is None.' % pattern)

        return cls.protImport
    
    @classmethod
    def runImportMovie1(cls, pattern):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(pattern, samplingRate=1.14, voltage=300,
                                  sphericalAberration=2.26,
                                  scannedPixelSize=None, magnification=50000)
    
    @classmethod
    def runImportMovie2(cls, pattern):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(pattern, samplingRate=1.4, voltage=300,
                                  sphericalAberration=2.7, scannedPixelSize=None,
                                  magnification=61000)


class TestOFAlignment(TestXmippBase):
    """This class check if the preprocessing micrographs protocol
    in Xmipp works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportMovie1(cls.movie1)
        cls.protImport2 = cls.runImportMovie2(cls.movie2)
    
    def runOFProtocol(self, movies, label="Default", saveMic=True,
                      saveMovie=False, useAlign=False):
        protOF = XmippProtOFAlignment(doSaveAveMic=saveMic,
                                      doSaveMovie=saveMovie,
                                      useAlignment=useAlign)
        protOF.setObjLabel(label)
        protOF.inputMovies.set(movies)
        self.launchProtocol(protOF)
        return protOF
    
    def testAlignOF1(self):
        protOF1 = self.runOFProtocol(self.protImport1.outputMovies,
                                     label="Movie MRC")
        self.assertIsNotNone(protOF1.outputMicrographs,
                             "SetOfMicrographs has not been created.")
    
    def testAlignOF2(self):
        protOF2 = self.runOFProtocol(self.protImport2.outputMovies,
                                     label="Movie EM")
        self.assertIsNotNone(protOF2.outputMicrographs,
                             "SetOfMicrographs has not been created.")
    
    def testAlignOFSaveMovieAndMic(self):
        protOF3 = self.runOFProtocol(self.protImport1.outputMovies,
                                     label="Save Movie", saveMovie=True)
        self.assertIsNotNone(protOF3.outputMovies,
                             "SetOfMovies has not been created.")
    
    def testAlignOFSaveMovieNoMic(self):
        protOF4 = self.runOFProtocol(self.protImport1.outputMovies,
                                     label="Save Movie", saveMic=False, 
                                     saveMovie=True)
        self.assertIsNotNone(protOF4.outputMovies,
                             "SetOfMovies has not been created.")
    
    def testAlignOFWAlignment(self):
        prot = XmippProtMovieCorr(doSaveAveMic=False)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)
        
        protOF5 = self.runOFProtocol(prot.outputMovies,
                                     label="Movie w Alignment",
                                     saveMic=False, 
                                     saveMovie=True)
        self.assertIsNotNone(protOF5.outputMovies,
                             "SetOfMovies has not been created.")



class TestCorrelationAlignment(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('movies')

    @classmethod
    def runImportMovies(cls, pattern, **kwargs):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        params = {'samplingRate': 1.14,
                  'voltage': 300,
                  'sphericalAberration': 2.7,
                  'magnification': 50000,
                  'scannedPixelSize': None,
                  'filesPattern': pattern
                  }
        if 'samplingRate' not in kwargs:
            del params['samplingRate']
            params['samplingRateMode'] = 0
        else:
            params['samplingRateMode'] = 1

        params.update(kwargs)

        protImport = cls.newProtocol(ProtImportMovies, **params)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        cls.protImport1 = cls.runImportMovies(cls.ds.getFile('qbeta/qbeta.mrc'),
                                              magnification=50000)
        cls.protImport2 = cls.runImportMovies(cls.ds.getFile('cct/cct_1.em'),
                                              magnification=61000)

    def _checkMicrographs(self, protocol):
        self.assertIsNotNone(getattr(protocol, 'outputMicrographs', None),
                             "Output SetOfMicrographs were not created.")

    def _checkAlignment(self, movie, goldRange, goldRoi):
        alignment = movie.getAlignment()
        range = alignment.getRange()
        msgRange = "Alignment range must be %s (%s) and it is %s (%s)"
        self.assertEqual(goldRange, range, msgRange
                         % (goldRange, range, type(goldRange), type(range)))
        roi = alignment.getRoi()
        msgRoi = "Alignment ROI must be %s (%s) and it is %s (%s)"
        self.assertEqual(goldRoi, roi,
                         msgRoi % (goldRoi, roi, type(goldRoi), type(roi)))

    def test_qbeta(self):
        prot = self.newProtocol(XmippProtMovieCorr)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1,7), [0, 0, 0, 0])

    def test_cct(self):
        prot = self.newProtocol(XmippProtMovieCorr,
                                doSaveMovie=True)
        prot.inputMovies.set(self.protImport2.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1,7), [0, 0, 0, 0])

    def test_qbeta_SkipCrop(self):
        prot = self.newProtocol(XmippProtMovieCorr,
                                alignFrame0=3, alignFrameN=5,
                                sumFrame0=3, sumFrameN=5,
                                cropOffsetX=10, cropOffsetY=10)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (3,5), [10, 10, 0, 0])


class TestAverageMovie(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('movies')

    @classmethod
    def runImportMovies(cls, pattern, **kwargs):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        params = {'samplingRate': 1.14,
                  'voltage': 300,
                  'sphericalAberration': 2.7,
                  'magnification': 50000,
                  'scannedPixelSize': None,
                  'filesPattern': pattern
                  }
        if 'samplingRate' not in kwargs:
            del params['samplingRate']
            params['samplingRateMode'] = 0
        else:
            params['samplingRateMode'] = 1

        params.update(kwargs)

        protImport = cls.newProtocol(ProtImportMovies, **params)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        cls.protImport1 = cls.runImportMovies(cls.ds.getFile('qbeta/qbeta.mrc'),
                                              magnification=50000)
        cls.protImport2 = cls.runImportMovies(cls.ds.getFile('cct/cct_1.em'),
                                              magnification=61000)
    
    def _checkMicrographs(self, protocol, goldDimensions):
        self.assertIsNotNone(getattr(protocol, 'outputMicrographs', None),
                             "Output SetOfMicrographs were not created.")
        mic = protocol.outputMicrographs[1]
        x, y, _ = mic.getDim()
        dims = (x, y)
        msgError = "The dimensions must be %s and it is %s"
        self.assertEqual(goldDimensions, dims, msgError % (goldDimensions, dims))

    def _checkAlignment(self, movie, goldRange, goldRoi):
        alignment = movie.getAlignment()
        range = alignment.getRange()
        msgRange = "Alignment range must be %s %s and it is %s (%s)"
        self.assertEqual(goldRange, range,
                         msgRange % (goldRange, range, type(goldRange), type(range)))
        roi = alignment.getRoi()
        msgRoi = "Alignment ROI must be %s (%s) and it is %s (%s)"
        self.assertEqual(goldRoi, roi,
                         msgRoi % (goldRoi, roi, type(goldRoi), type(roi)))
    
    def test_qbeta(self):
        prot = self.newProtocol(XmippProtMovieCorr,
                                alignFrame0=3, alignFrameN=5,
                                cropOffsetX=10, cropOffsetY=10,
                                doSaveAveMic=False)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)
        
        self._checkAlignment(prot.outputMovies[1],
                             (3,5), [10, 10, 0, 0])
        
        protAverage = self.newProtocol(XmippProtMovieAverage,
                                       cropRegion=1)
        protAverage.inputMovies.set(prot.outputMovies)
        protAverage.setObjLabel('average w alignment info')
        self.launchProtocol(protAverage)
        
        self._checkMicrographs(protAverage, (4086,4086))
        protAverage2 = self.newProtocol(XmippProtMovieAverage,
                                        sumFrame0=1, sumFrameN=1,
                                       cropRegion=2)
        protAverage2.inputMovies.set(prot.outputMovies)
        protAverage2.setObjLabel('average w alignment')
        self.launchProtocol(protAverage2)

        self._checkMicrographs(protAverage2, (4096,4096))

    def test_cct(self):
        protAverage = self.newProtocol(XmippProtMovieAverage,
                                       cropRegion=2,
                                       cropOffsetX=10, cropOffsetY=10,
                                       cropDimX=1500, cropDimY=1500)
        protAverage.inputMovies.set(self.protImport2.outputMovies)
        protAverage.setObjLabel('average imported movies')
        self.launchProtocol(protAverage)

        self._checkMicrographs(protAverage, (1500,1500))

    def test_cct2(self):
        protAverage = self.newProtocol(XmippProtMovieAverage,
                                       cropRegion=1)
        protAverage.inputMovies.set(self.protImport2.outputMovies)
        protAverage.setObjLabel('average imported movies')
        self.launchProtocol(protAverage)

        self._checkMicrographs(protAverage, (4096,4096))


class TestEstimateGain(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        ds = DataSet.getDataSet('movies')

        # Reduce input movie size to speed-up gain computation
        ih = ImageHandler()
        inputFn = ds.getFile('ribo/Falcon_2012_06_12-14_33_35_0_movie.mrcs')
        outputFn = cls.proj.getTmpPath(abspath(basename(inputFn)))

        frameImg = ih.createImage()
        xdim, ydim, zdim, ndim = ih.getDimensions(inputFn)
        n = max(zdim, ndim) / 2 # also half of the frames
        print "Scaling movie: %s -> %s" % (inputFn, outputFn)
        pwutils.cleanPath(outputFn)
        for i in range(1, n+1):
            frameImg.read((i, inputFn))
            frameImg.scale(xdim/2, ydim/2)
            frameImg.write((i, outputFn))

        args = cls.getArgs(outputFn)
        cls.protImport = cls.newProtocol(ProtImportMovies, **args)
        cls.launchProtocol(cls.protImport)

    @classmethod
    def getArgs(self, filesPath, pattern=''):
        return {'importFrom': ProtImportMovies.IMPORT_FROM_FILES,
                'filesPath': filesPath,
                'filesPattern': pattern,
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 300,
                'samplingRate': 3.54 * 2
                }

    def test_estimate(self):
        protGain = self.newProtocol(XmippProtMovieGain,
                                    objLabel='estimate gain')

        p = Pointer(self.protImport.outputMovies, extended=1)
        protGain.inputMovies.append(p)

        self.launchProtocol(protGain)