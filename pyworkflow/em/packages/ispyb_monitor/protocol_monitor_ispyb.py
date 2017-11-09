# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es) [1]
# *              Kevin Savage (kevin.savage@diamond.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] Diamond Light Source, Ltd
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

from os import environ
from os.path import realpath, join, dirname, exists, basename, abspath
from collections import OrderedDict
import math
from datetime import datetime

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.em import ImageHandler
from pyworkflow.em.protocol import ProtMonitor, Monitor, PrintNotifier
from pyworkflow.em.protocol import ProtImportMovies, ProtAlignMovies, ProtCTFMicrographs
from pyworkflow.gui import getPILImage
from pyworkflow.protocol.constants import STATUS_RUNNING

from sys import float_info

class ProtMonitorISPyB(ProtMonitor):
    """ Monitor to communicated with ISPyB system at Diamond.
    """
    _label = 'monitor to ISPyB'
    _lastUpdateVersion = VERSION_1_1

    def _defineParams(self, form):
        ProtMonitor._defineParams(self, form)

        group = form.addGroup('Experiment')
        group.addParam('visit', params.StringParam,
                      label="Visit",
                      help="Visit")

        form.addParam('db', params.EnumParam,
                      choices=["production", "devel", "test"],
                      label="Database",
                      help="Select which ISPyB database you want to use.")

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        inputProtocols = self.getInputProtocols()

        monitor = MonitorISPyB(self, workingDir=self._getPath(),
                               samplingInterval=self.samplingInterval.get(),
                               monitorTime=10000,
                               inputProtocols=inputProtocols,
                               visit=self.visit.get(),
                               dbconf=self.db.get(),
                               project=self.getProject())

        monitor.addNotifier(PrintNotifier())
        monitor.loop()

class MonitorISPyB(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocol, **kwargs):
        Monitor.__init__(self, **kwargs)
        self.protocol = protocol
        self.dataCollection = OrderedDict()
        self.movies = OrderedDict()
        self.motion_corrections = OrderedDict()
        self.motion_correction_drift = OrderedDict()
        self.ctfs = OrderedDict()
        self.numberOfFrames = None
        self.imageGenerator = None
        self.dcId = None
        self.imageNumber = 1
        self.previousParams = None
        self.visit = kwargs['visit']
        self.dbconf = kwargs['dbconf']
        self.project = kwargs['project']
        self.inputProtocols = self._sortInputProtocols(kwargs['inputProtocols'])
        self.ispybDb = ISPyBdb(experimentParams={'visit': self.visit})
    @staticmethod
    def _sortInputProtocols(protList):
        # we need sorted input protocols in order to process objects correctly
        movieProts = []
        alignProts = []
        ctfProts   = []
        for p in protList:
            if isinstance(p, ProtImportMovies):
                movieProts.append(p)
            elif isinstance(p, ProtAlignMovies):
                alignProts.append(p)
            elif isinstance(p, ProtCTFMicrographs):
                ctfProts.append(p)
        sortedProts = movieProts + alignProts + ctfProts
        return sortedProts

    def step(self):
        self.info("MonitorISPyB: only one step")

        prots = [self.getUpdatedProtocol(p) for p in self.inputProtocols]
        finished = [] # store True if protocol not running
        updateImageIds = []  # Store obj ids that have changes
        updateAlignIds = []  # Store obj ids that have changes
        updateCTFIds = []  # Store obj ids that have changes

        for prot in prots:
            self.info("protocol: %s" % prot.getRunName())
            if isinstance(prot, ProtImportMovies) and hasattr(prot, 'outputMovies'):
                self.create_movie_params(prot, updateImageIds)
            elif isinstance(prot, ProtAlignMovies) and hasattr(prot, 'outputMicrographs'):
                self.update_align_params(prot, updateAlignIds)
            elif isinstance(prot, ProtCTFMicrographs) and hasattr(prot, 'outputCTF'):
                self.update_ctf_params(prot, updateCTFIds)

            finished.append(prot.getStatus() != STATUS_RUNNING)

        dcParams = self.ispybDb.get_data_collection_params()
        self.safe_update(dcParams, self.dataCollection)
        if self.dcId:
            self.safe_update(dcParams, self.previousParams)
            dcParams['id'] = self.dcId
            dcParams['endtime'] = self.now()
            self.info("writing datacollection: %s" + str(dcParams))
            self.ispybDb.update_data_collection(dcParams)
        else:
            dcParams['starttime'] = self.now()
            dcParams['endtime'] = self.now()
            self.info("writing datacollection: %s" + str(dcParams))
            self.dcId = self.ispybDb.update_data_collection(dcParams)

        self.previousParams = dcParams

        for itemId in set(updateAlignIds):
            if 'autoProcProgramId' not in self.motion_corrections[itemId]:
                program = self.ispybDb.get_program_params()
                self.safe_update(program, self.motion_corrections[itemId])
                program['starttime'] = self.now()
                program_id = self.ispybDb.update_program(program)
                self.motion_corrections[itemId]['autoProcProgramId'] = program_id

            if 'imageNumber' not in self.motion_corrections[itemId]:
                self.motion_corrections[itemId]['imageNumber'] = self.imageNumber
                self.imageNumber += 1

            motionParams = self.ispybDb.get_motion_correction_params()
            self.safe_update(motionParams, self.motion_corrections[itemId])
            motionParams['dataCollectionId'] = self.dcId
            self.info("writing motion correction: %s" + str(motionParams))
            motionCorrectionId = self.ispybDb.update_motion_correction(motionParams)
            self.info("wrote motion correction: %s" + str(motionCorrectionId))
            self.motion_corrections[itemId]['motionCorrectionId'] = motionCorrectionId

            drift = self.motion_correction_drift[itemId]
            xs = drift['xs']
            ys = drift['ys']
            for (i, (x, y)) in enumerate(zip(xs, ys)):
                drift_params = self.ispybDb.get_motion_correction_drift_params()
                drift_params['motionCorrectionId'] = motionCorrectionId
                drift_params['frameNumber'] = i + 1
                drift_params['deltaX'] = x
                drift_params['deltaY'] = y
                self.ispybDb.update_motion_correction_drift(drift_params)

        for itemId in set(updateCTFIds):
            if 'autoProcProgramId' not in self.ctfs[itemId]:
                program = self.ispybDb.get_program_params()
                self.safe_update(program, self.ctfs[itemId])
                program['starttime'] = self.now()
                program_id = self.ispybDb.update_program(program)
                self.ctfs[itemId]['autoProcProgramId'] = program_id

            ctfParams = self.ispybDb.get_ctf_params()
            self.safe_update(ctfParams, self.ctfs[itemId])
            ctfParams['motionCorrectionId'] = self.motion_corrections[itemId]['motionCorrectionId']
            self.info("writing ctf: %s" + str(ctfParams))
            ctfId = self.ispybDb.update_ctf(ctfParams)
            self.info("wrote ctf: %s" + str(ctfId))
            self.ctfs[itemId]['ctfId'] = ctfId

        if all(finished):
            self.info("All finished, closing ISPyBDb connection")
            self.ispybDb.disconnect()

        return all(finished)



    def now(self):
        return datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')

    def iter_updated_set(self, objSet):
        objSet.load()
        objSet.loadAllProperties()
        for obj in objSet:
            yield obj
        objSet.close()

    def find_ispyb_path(self, input_file):
        """ Given a visit, find the path where png images should be stored. """
        if pwutils.envVarOn('SCIPIONBOX_ISPYB_ON'):
            p = realpath(join(self.project.path, input_file))
            while p and not p.endswith(self.visit):
                p = dirname(p)
            return join(p, '.ispyb')
        else:
            return self.protocol._getExtraPath()

    def create_movie_params(self, prot, updateIds):

        for movie in self.iter_updated_set(prot.outputMovies):
            movieId = movie.getObjId()
            if movieId in self.movies:  # this movie has been processed, skip
                continue
            movieFn = movie.getFileName()
            if self.numberOfFrames is None:
                self.numberOfFrames = movie.getNumberOfFrames()
                images_path = self.find_ispyb_path(movieFn)
                self.imageGenerator = ImageGenerator(self.project.path,
                                                     images_path,
                                                     bigThumb=True,
                                                     smallThumb=512)
            acquisition = movie.getAcquisition()

            self.movies[movieId] = {
                'experimenttype': 'single particle',
                'run_status': 'DataCollection Successful',
                'imgdir': abspath(dirname(movieFn)),
                'imgsuffix': pwutils.getExt(movieFn),
                'file_template': pwutils.removeBaseExt(movieFn) + '#####' + pwutils.getExt(movieFn),
                'file_location': abspath(dirname(movieFn)),
                'filename': movieFn,
                'n_passes': self.numberOfFrames,
                'magnification': acquisition.getMagnification(),
                'total_absorbed_dose': acquisition.getDoseInitial() + (acquisition.getDosePerFrame() * self.numberOfFrames),
                'wavelength': self.convert_volts_to_debroglie_wavelength(acquisition.getVoltage())
            }
            self.dataCollection.update(self.movies[movieId])
            updateIds.append(movieId)

    def safe_update(self, target, source):
        for key in source:
            if source[key] is not None:
                try:
                    target[key] = source[key]
                except KeyError:
                    pass

    @staticmethod
    def convert_volts_to_debroglie_wavelength(volts):
        rest_mass_of_electron = 9.10938356e-31
        planks_constant = 6.62607004e-34
        speed_of_light = 2.99792458e8
        charge_of_electron = 1.60217662e-19

        u = volts * 1e3
        return (planks_constant /
                (2 * rest_mass_of_electron * charge_of_electron * u + ((charge_of_electron * u / speed_of_light)** 2))**0.5
            ) * 1e10


    def update_align_params(self, prot, updateIds):
        for mic in self.iter_updated_set(prot.outputMicrographs):
            micId = mic.getObjId()
            if self.movies.get(micId, None) is not None:
                if micId in self.motion_corrections and 'comments' in self.motion_corrections[micId]:  # skip if we already have align info
                    continue

                micFn = mic.getFileName()
                renderable_image = self.imageGenerator.generate_image(micFn, micFn)

                xs, ys = prot._getMovieShifts(mic)
                totalMotion = sum(map(lambda p: (p[0]**2 + p[1]**2) ** 0.5, zip(xs, ys)))

                if micId not in self.motion_corrections:
                    self.motion_corrections[micId] = {}

                self.motion_corrections[micId].update({
                    'firstFrame': prot.alignFrame0.get(),
                    'lastFrame': prot.alignFrameN.get(),
                    'micrographSnapshotFullPath': renderable_image,
                    'micrographFullPath': micFn,
                    'driftPlotFullPath': getattr(prot, '_getPlotGlobal', lambda x: None)(mic),
                    'totalMotion': totalMotion,
                    'averageMotionPerFrame': abs(totalMotion / max(len(xs), 1)),
                    'comments': 'aligned',
                    'patchesUsed': (prot.patchX.get() + prot.patchY.get())/2,
                    'programs': getattr(prot, '_label', lambda x: None),
                    'status': 1
                })

                if micId not in self.motion_correction_drift:
                    self.motion_correction_drift[micId] = {}

                self.motion_correction_drift[micId].update({
                    'xs': xs,
                    'ys': ys
                })

                self.dataCollection['xtal_snapshot1'] = renderable_image
                print('%d has new align info' % micId)
                updateIds.append(micId)

    def update_ctf_params(self, prot, updateIds):
        for ctf in self.iter_updated_set(prot.outputCTF):
            micId = ctf.getObjId()
            if self.motion_corrections.get(micId, None) is not None:
                if micId in self.ctfs and 'status' in self.ctfs[micId]:  # skip if we already have align info
                    continue
                micFn = ctf.getMicrograph().getFileName()
                psdName = pwutils.replaceBaseExt(micFn, 'psd.jpeg')
                psdFn = ctf.getPsdFile()
                psdPng = self.imageGenerator.generate_image(psdFn, psdName)

                if micId not in self.ctfs:
                    self.ctfs[micId] = {}

                self.ctfs[micId].update({
                    'estimatedDefocus': ((ctf.getDefocusV() + ctf.getDefocusU())/2),
                    'astigmatism': ctf.getDefocusRatio(),
                    'estimatedResolution': ctf.getResolution(),
                    'astigmatismAngle': ctf.getDefocusAngle(),
                    'fftTheoreticalFullPath': psdPng,
                    'programs': getattr(prot, '_label', lambda x: None),
                    'status': 1
                })
                print('%d has new ctf info' % micId)
                updateIds.append(micId)


class ImageGenerator:
    def __init__(self, project_path, images_path,
                 bigThumb=None, smallThumb=None):
        self.project_path = project_path
        self.images_path = images_path
        self.ih = ImageHandler()
        self.img = self.ih.createImage()
        self.bigThumb = bigThumb
        self.smallThumb = smallThumb

    def generate_image(self, input_file, outputName=None):
        output_root = join(self.images_path, basename(outputName))
        output_file = output_root + '.jpg'

        print "Generating image: ", output_file

        if not exists(output_file):
            from PIL import Image
            self.img.read(join(self.project_path, input_file))
            pimg = getPILImage(self.img)

            pwutils.makeFilePath(output_file)
            if self.bigThumb:
                pimg.save(output_file, "JPEG")

            if self.smallThumb:
                pimg.thumbnail((self.smallThumb, self.smallThumb), Image.ANTIALIAS)
                pimg.save(output_root + 't.jpg', "JPEG")

        return output_file


def _loadMeanShifts(self, movie):
    alignMd = md.MetaData(self._getOutputShifts(movie))
    meanX = alignMd.getColumnValues(md.MDL_OPTICALFLOW_MEANX)
    meanY = alignMd.getColumnValues(md.MDL_OPTICALFLOW_MEANY)

    return meanX, meanY


def _saveAlignmentPlots(self, movie):
    """ Compute alignment shifts plot and save to file as a png image. """
    meanX, meanY = self._loadMeanShifts(movie)
    plotter = createAlignmentPlot(meanX, meanY)
    plotter.savefig(self._getPlotCart(movie))


def createAlignmentPlot(meanX, meanY):
    """ Create a plotter with the cumulative shift per frame. """
    sumMeanX = []
    sumMeanY = []
    figureSize = (8, 6)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()

    preX = 0.0
    preY = 0.0
    sumMeanX.append(0.0)
    sumMeanY.append(0.0)
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Cartesian representation')
    ax.set_xlabel('Drift x (pixels)')
    ax.set_ylabel('Drift y (pixels)')
    ax.plot(0, 0, 'yo-')
    i = 1
    for x, y in izip(meanX, meanY):
        preX += x
        preY += y
        sumMeanX.append(preX)
        sumMeanY.append(preY)
        #ax.plot(preX, preY, 'yo-')
        ax.text(preX-0.02, preY+0.02, str(i))
        i += 1

    ax.plot(sumMeanX, sumMeanY, color='b')
    ax.plot(sumMeanX, sumMeanY, 'yo')

    plotter.tightLayout()

    return plotter


class ISPyBdb:
    """ This is a Facade to provide access to the ispyb_api to store movies."""
    def __init__(self, experimentParams):
        self.experimentParams = experimentParams

        from ispyb import factory

        config = environ['ISPYB_CONFIG']
        self.dbconnection = factory.create_connection(config)
        self.core = factory.create_data_area(factory.DataAreaType.CORE, self.dbconnection)
        self.mxacquisition = factory.create_data_area(factory.DataAreaType.MXACQUISITION, self.dbconnection)
        self.emacquisition = factory.create_data_area(factory.DataAreaType.EMACQUISITION, self.dbconnection)
        self.mxdatareduction = factory.create_data_area(factory.DataAreaType.MXPROCESSING, self.dbconnection)

        self._create_group()

    def _create_group(self):
        self.visit_id = self.core.retrieve_visit_id(self.experimentParams['visit'])
        params = self.mxacquisition.get_data_collection_group_params()
        params['parentid'] = self.visit_id
        self.group_id = self.mxacquisition.insert_data_collection_group(self.convert_float_types(params).values())

    def get_data_collection_params(self):
        params = self.mxacquisition.get_data_collection_params()
        params['parentid'] = self.group_id
        params['visitid'] = self.visit_id
        return params

    def update_data_collection(self, params):
        return self.mxacquisition.insert_data_collection(self.convert_float_types(params).values())

    def get_image_params(self):
        return self.mxacquisition.get_image_params()

    def update_image(self, params):
        return self.mxacquisition.update_image(self.convert_float_types(params).values())

    def get_motion_correction_params(self):
        return self.emacquisition.get_motion_correction_params()

    def update_motion_correction(self, params):
        return self.emacquisition.insert_motion_correction(self.convert_float_types(params).values())

    def get_motion_correction_drift_params(self):
        return self.get_motion_correction_drift_params()

    def update_motion_correction_drift(self, params):
        raise self.update_motion_correction_drift(self.convert_float_types(params).values())

    def get_ctf_params(self):
        return self.emacquisition.get_ctf_params()

    def update_ctf(self, params):
        return self.emacquisition.insert_ctf(self.convert_float_types(params).values())

    def get_program_params(self):
        return self.mxdatareduction.get_program_params()

    def update_program(self, params):
        return self.mxdatareduction.insert_program(self.convert_float_types(params).values())

    def convert_float_types(self, params):
        for k in params:
            try:
                v = params[k]
                if math.isnan(v):
                    params[k] = None
                elif math.isinf(v):
                    params[k] = math.copysign(float_info.max, v)
            except TypeError:
                pass
        return params

    def disconnect(self):
        if self.dbconnection:
            self.dbconnection.disconnect()
