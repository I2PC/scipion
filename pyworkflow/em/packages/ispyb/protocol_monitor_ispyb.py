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

from os.path import realpath, join, dirname, exists, basename
from collections import OrderedDict

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.em import ImageHandler
from pyworkflow.em.protocol import ProtMonitor, Monitor, PrintNotifier
from pyworkflow.em.protocol import ProtImportMovies, ProtAlignMovies, ProtCTFMicrographs
from pyworkflow.gui import getPILImage

from ispyb_proxy import ISPyBProxy


class ProtMonitorISPyB(ProtMonitor):
    """ Monitor to communicated with ISPyB system at Diamond.
    """
    _label = 'monitor to ISPyB'

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
        monitor = MonitorISPyB(self, workingDir=self._getPath(),
                               samplingInterval=self.samplingInterval.get(),
                               monitorTime=100)

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
        self.allIds = OrderedDict()
        self.numberOfFrames = None
        self.imageGenerator = None
        self.visit = self.protocol.visit.get()
        self.project = self.protocol.getProject()

    def step(self):
        self.info("MonitorISPyB: only one step")

        prot = self.protocol

        proxy = ISPyBProxy(["prod", "dev", "test"][prot.db.get()],
                           experimentParams={'visit': prot.visit.get()})


        runs = [p.get() for p in self.protocol.inputProtocols]
        g = self.project.getGraphFromRuns(runs)

        nodes = g.getRoot().iterChildsBreadth()

        allParams = OrderedDict()

        for n in nodes:
            prot = n.run
            self.info("protocol: %s" % prot.getRunName())

            if isinstance(prot, ProtImportMovies):
                self.create_movie_params(prot, allParams)
            elif isinstance(prot, ProtAlignMovies) and hasattr(prot, 'outputMicrographs'):
                self.update_align_params(prot, allParams)
            elif isinstance(prot, ProtCTFMicrographs):
                self.update_ctf_params(prot, allParams)

        for itemId, params in allParams.iteritems():
            ispybId = proxy.sendMovieParams(params)
            # Use -1 as a trick when ispyb is not really used and id is None
            self.allIds[itemId] = ispybId or -1

        self.info("Closing proxy")
        proxy.close()

        return False

    def iter_updated_set(self, objSet):
        objSet.load()
        objSet.loadAllProperties()
        for obj in objSet:
            yield obj
        objSet.close()

    def find_ispyb_path(self, input_file):
        """ Given a visit, find the path where png images should be stored. """
        if pwutils.envVarOn('SCIPIONBOX_ISPYB_ON'):
            p = realpath(join(self.project_path, input_file))
            while p and not p.endswith(self.visit):
                p = dirname(p)
            return join(p, '.ispyb')
        else:
            return self.protocol._getExtraPath()

    def create_movie_params(self, prot, allParams):

        for movie in self.iter_updated_set(prot.outputMovies):
            movieFn = movie.getFileName()
            if self.numberOfFrames is None:
                self.numberOfFrames = movie.getNumberOfFrames()
                images_path = self.find_ispyb_path(movieFn)
                self.imageGenerator = ImageGenerator(self.project.path,
                                                     images_path,
                                                     smallThumb=512)

            movieId = movie.getObjId()

            allParams[movieId] = {
                'id': self.allIds.get(movieId, None),
                'imgdir': dirname(movieFn),
                'imgprefix': pwutils.removeBaseExt(movieFn),
                'imgsuffix': pwutils.getExt(movieFn),
                'file_template': movieFn,
                'n_images': self.numberOfFrames
             }


    def update_align_params(self, prot, allParams):
        for mic in self.iter_updated_set(prot.outputMicrographs):
            micFn = mic.getFileName()
            renderable_image = self.imageGenerator.generate_image(micFn, micFn)

            allParams[mic.getObjId()].update({
                'comments': 'aligned',
                'xtal_snapshot1':renderable_image
            })

    def update_ctf_params(self, prot, allParams):
        for ctf in self.iter_updated_set(prot.outputCTF):
            micFn = ctf.getMicrograph().getFileName()
            psdName = pwutils.replaceBaseExt(micFn, 'psd.png')
            psdFn = ctf.getPsdFile()
            psdPng = self.imageGenerator.generate_image(psdFn, psdName)
            allParams[ctf.getObjId()].update({
            'min_defocus': ctf.getDefocusU(),
            'max_defocus': ctf.getDefocusV(),
            'amount_astigmatism': ctf.getDefocusRatio()
            })


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
        output_file = output_root + '.png'

        print "Generating image: ", output_file

        if not exists(output_file):
            from PIL import Image
            self.img.read(join(self.project_path, input_file))
            pimg = getPILImage(self.img)

            pwutils.makeFilePath(output_file)
            if self.bigThumb:
                pimg.save(output_file, "PNG")

            if self.smallThumb:
                pimg.thumbnail((self.smallThumb, self.smallThumb), Image.ANTIALIAS)
                pimg.save(output_root + 't.png', "PNG")

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


class FileNotifier():
    def __init__(self, filename):
        self.f = open(filename, 'w')

    def notify(self, title, message):
        print >> self.f, title, message
        self.f.flush()