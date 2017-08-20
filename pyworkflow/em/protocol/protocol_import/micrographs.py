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

import os
from os.path import join, basename
import re
from datetime import timedelta
import socket
import select
import shlex

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message
from pyworkflow.em import ImageHandler
from pyworkflow.em.constants import SAMPLING_FROM_IMAGE, SAMPLING_FROM_SCANNER

from images import ProtImportImages



class ProtImportMicBase(ProtImportImages):
    """ Just to have a base class to both 
    ProtImportMicrographs and ProtImportMovies
    """
    _checkStacks = False

    def _defineAcquisitionParams(self, form):
        group = ProtImportImages._defineAcquisitionParams(self, form)

        group.addParam('samplingRateMode', params.EnumParam,
                       choices=[Message.LABEL_SAMP_MODE_1,
                                Message.LABEL_SAMP_MODE_2],
                       default=SAMPLING_FROM_IMAGE,
                       label=Message.LABEL_SAMP_MODE,
                       help=Message.TEXT_SAMP_MODE)

        group.addParam('samplingRate', params.FloatParam,  default=1.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_IMAGE, 
                       label=Message.LABEL_SAMP_RATE,
                       help=Message.TEXT_SAMP_RATE)

        group.addParam('scannedPixelSize', params.FloatParam, default=7.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER,
                       label=Message.LABEL_SCANNED,
                       help='')
        return group
        
    def setSamplingRate(self, micSet):
        """ Set the sampling rate to the given set. """
        if self.samplingRateMode == SAMPLING_FROM_IMAGE:
            micSet.setSamplingRate(self.samplingRate.get())
        else:
            micSet.setScannedPixelSize(self.scannedPixelSize.get())

    def _acquisitionWizardCondition(self):
        """ By default this wizard will appears only when we import from
        a format that is not from files.
        But movie-import also can have a wizard to read from FEI xml files. """
        return 'True'

    def loadAcquisitionInfo(self):
        """ Return a proper acquistionInfo (dict)
        or an error message (str).
        """
        if self.importFrom != self.IMPORT_FROM_FILES:
            return ProtImportImages.loadAcquisitionInfo(self)

        result = "Could not find acquistion information"

        for fileName, fileId in self.iterFiles():
            baseName = pwutils.removeExt(fileName)
            xml1 = baseName.replace('_frames', '.xml')
            if os.path.exists(xml1):
                result = self._parseXML(xml1)
            else:
                xml2 = baseName + '.xml'
                result = self._parseXML(xml2)

        return result

    def _parseXML(self, fileName):
        """ Parse micrograph XML files from FEI. """
        import xml.etree.ElementTree as ET

        # get context
        context = iter(ET.iterparse(fileName,
                                    events=('start', 'end')))

        labels = {'AccelerationVoltage': 'voltage',
                  'InstrumentModel': 'InstrumentModel',
                  'NominalMagnification': 'magnification'}

        # acq['amplitudeContrast'] = None
        # acq['sphericalAberration'] = None
        acq = {}

        def get(key, elem):
            acq[labels[key]] = elem.text

        pixelSize = False

        for event, elem in context:

            if event == 'start':
                if 'pixelSize' in elem.tag:
                    print "started: pixelSize"
                    pixelSize = True

            elif event == 'end':
                for l in labels:
                    if '}%s' % l in elem.tag:
                        get(l, elem)

                if '}numericValue' in elem.tag and pixelSize:
                    acq['samplingRate'] = float(elem.text) * 10e+09  # Convert to A
                    pixelSize = False

            else:
                raise Exception("Unknown event type %s" % event)

        # Correct for units conversion
        acq['voltage'] = float(acq['voltage']) / 1000.

        return acq


    
class ProtImportMicrographs(ProtImportMicBase):
    """Protocol to import a set of micrographs to the project"""
    _label = 'import micrographs'
    _outputClassName = 'SetOfMicrographs' 
    
    IMPORT_FROM_EMX = 1
    IMPORT_FROM_XMIPP3 = 2
    IMPORT_FROM_SCIPION = 3

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        choices = ProtImportImages._getImportChoices(self)
        return choices + ['emx', 'xmipp3', 'scipion']
    
    def _defineImportParams(self, form):
        """ Just redefine to put some import parameters
        before the acquisition related parameters.
        """
        form.addParam('emxFile', params.FileParam,
              condition = '(importFrom == %d)' % self.IMPORT_FROM_EMX,
              label='Input EMX file',
              help="Select the EMX file containing micrographs information.\n"
                   "See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]")
        
        form.addParam('mdFile', params.FileParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_XMIPP3,
                      label='Micrographs metadata file',
                      help="Select the micrographs Xmipp metadata file.\n"
                           "It is usually a _micrograph.xmd_ file result\n"
                           "from import, preprocess or downsample protocols.")
        
        form.addParam('sqliteFile', params.FileParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_SCIPION,
                      label='Micrographs sqlite file',
                      help="Select the micrographs sqlite file.\n")
    
    #--------------------------- INSERT functions ---------------------------------------------------
    def _insertAllSteps(self):
        importFrom = self.importFrom.get()
        ci = self.getImportClass()
        
        if ci is None:
            ProtImportMicBase._insertAllSteps(self)
        else:
            self._insertFunctionStep('importMicrographsStep', importFrom,
                                     self.importFilePath)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importMicrographsStep(self, importFrom, *args):
        ci = self.getImportClass()
        ci.importMicrographs()
        
        summary = "Import from *%s* file:\n" % self.getEnumText('importFrom')
        summary += self.importFilePath + '\n'
        
        if self.hasAttribute('outputParticles'):
            particles = self.outputParticles
            summary += '   Particles: *%d* (ctf=%s, alignment=%s)\n' % (particles.getSize(),
                                                                        particles.hasCTF(),
                                                                        particles.getAlignment())
                                                                      
        if self.hasAttribute('outputCoordinates'): # EMX files can contain only Coordinates information
            summary += '   Coordinates: *%d* \n' % (self.outputCoordinates.getSize())
            
        if self.hasAttribute('outputMicrographs'): # EMX files can contain only Coordinates information
            summary += '   Micrographs: *%d* \n' % (self.outputMicrographs.getSize())
        
        if self.copyFiles:
            summary += '\n_WARNING_: Binary files copied into project (extra disk space)'
            
        self.summaryVar.set(summary)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        from pyworkflow.em.convert import ImageHandler
        ci = self.getImportClass()
        if ci is None:
            errors = ProtImportMicBase._validate(self)
            for micFn, _ in self.iterFiles():
                imgh = ImageHandler()
                if imgh.isImageFile(micFn):
                    _, _, z, n = imgh.getDimensions(micFn)
                    if n > 1 or z > 1:
                        errors.append("The protocol not support micrographs stored in stacks. "
                                      "If you want to obtain your micrographs individually, "
                                      "you can run the following command:\n"
                                      "scipion run scipion_directory/scripts/split_stacks.py --files *your files* --ext *extension*")
                # JMRT: only check the first image, for large dataset
                # even reading the header can take a while
                break 
            return errors
            
        else:
            return ci.validateMicrographs()
    
    def _summary(self):
        if self.importFrom == self.IMPORT_FROM_FILES:
            return ProtImportMicBase._summary(self)
        else:
            return [self.summaryVar.get('No summary information.')]
    
    #--------------------------- UTILS functions -------------------------------
    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        if self.importFrom == self.IMPORT_FROM_EMX:
            from pyworkflow.em.packages.emxlib import EmxImport
            self.importFilePath = self.emxFile.get('').strip()
            return EmxImport(self, self.importFilePath)
        elif self.importFrom == self.IMPORT_FROM_XMIPP3:
            from pyworkflow.em.packages.xmipp3.dataimport import XmippImport
            self.importFilePath = self.mdFile.get('').strip()
            return XmippImport(self, self.mdFile.get())
        elif self.importFrom == self.IMPORT_FROM_SCIPION:
            from dataimport import ScipionImport
            self.importFilePath = self.sqliteFile.get('').strip()
            return ScipionImport(self, self.importFilePath) 
        else:
            self.importFilePath = ''
            return None       


class ProtImportMovies(ProtImportMicBase):
    """ Protocol to import a set of movies (from direct detector cameras)
    to the project.
    """
    _label = 'import movies'
    _outputClassName = 'SetOfMovies'

    def __init__(self, **kwargs):
        ProtImportMicBase.__init__(self, **kwargs)
        self.serverSocket = None
        self.connectionList = None

    def _defineAcquisitionParams(self, form):
        group = ProtImportMicBase._defineAcquisitionParams(self, form)

        line = group.addLine('Dose (e/A^2)',
                             help="Initial accumulated dose (usually 0) and "
                                  "dose per frame. ")

        line.addParam('doseInitial', params.FloatParam, default=0,
                      label='Initial')

        line.addParam('dosePerFrame', params.FloatParam, default=None,
                      allowsNull=True,
                      label='Per frame')

        form.addParam('gainFile', params.FileParam,
                      label='Gain image', 
                      help='A gain reference related to a set of movies'
                           ' for gain correction')

        form.addParam('darkFile', params.FileParam,
                      label='Dark image', 
                      help='A dark image related to a set of movies')

    def _defineParams(self, form):
        ProtImportMicBase._defineParams(self, form)

        form.addSection('Frames')

        streamingConditioned = "dataStreaming"
        framesCondition = "inputIndividualFrames"

        form.addParam('inputIndividualFrames', params.BooleanParam,
                      default=False,
                      label="Input individual frames?",
                      help="Select Yes if movies are acquired in individual "
                           "frame files. ")

        form.addParam('numberOfIndividualFrames', params.IntParam,
                      condition=framesCondition,
                      label='Number of frames',
                      help='Provide how many frames are per movie. ')

        form.addParam('stackFrames', params.BooleanParam,
                      default=False, condition=framesCondition,
                      label="Create movie stacks?",
                      help="Select Yes if you want to create a new stack for "
                           "each movies with its frames. ")

        # This is not working so for now its hidden
        form.addParam('writeMoviesInProject', params.BooleanParam,
                      default=False, condition=framesCondition + " and stackFrames",
                      label="Write stacks in the project folder?",
                      help="If Yes, the created stack files will be written "
                           "in the project folder. By default the movies will "
                           "be written in the same place where input frames are.")
        # form.addParam('writeMoviesInProject', params.HiddenBooleanParam,
        #               default=False, condition=framesCondition + " and stackFrames")


        form.addParam('movieSuffix', params.StringParam,
                      default='_frames.mrcs',
                      condition=framesCondition + " and stackFrames",
                      label="Movie suffix",
                      help="Suffix added to the output movie filename."
                           "Use the extension to select the format ("
                           "e.g., .mrcs, .stk)")

        form.addParam('deleteFrames', params.BooleanParam,
                      default=False,
                      condition=framesCondition + " and stackFrames",
                      label="Delete frame files?",
                      help="Select Yes if you want to remove the individual "
                           "frame files after creating the movie stack. ")

        streamingSection = form.getSection('Streaming')
        streamingSection.addParam('streamingSocket', params.BooleanParam, default=False,
                                  condition=streamingConditioned,
                                  expertLevel=params.LEVEL_ADVANCED,
                                  label="Use streaming socket",
                                  help="Use a socket to discover new files instead of polling\n"
                                       "your directory.\n")

        streamingSection.addParam('socketPort', params.IntParam, default=5000,
                                  condition=streamingConditioned + ' and streamingSocket',
                                  expertLevel=params.LEVEL_ADVANCED,
                                  label="Socket port",
                                  help="Port to use for the streaming socket.\n")

    # --------------------------- INSERT functions ---------------------------------------------------
    def _insertAllSteps(self):
        # Only the import movies has property 'inputIndividualFrames'
        # so let's query in a non-intrusive manner
        inputIndividualFrames = getattr(self, 'inputIndividualFrames', False)

        if self.dataStreaming or inputIndividualFrames:
            if self.streamingSocket:
                self.launchSocket()
            funcName = 'importImagesStreamStep'
        else:
            funcName = 'importImagesStep'

        self._insertFunctionStep(funcName, self.getPattern(),
                                 self.voltage.get(),
                                 self.sphericalAberration.get(),
                                 self.amplitudeContrast.get(),
                                 self.magnification.get())

    # --------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        """Overwriting to skip file validation if streaming with socket"""
        if self.streamingSocket:
            errors = []
        else:
            errors = ProtImportMicBase._validate(self)
        return errors

    # --------------------------- UTILS functions -------------------------------
    def setSamplingRate(self, movieSet):
        ProtImportMicBase.setSamplingRate(self, movieSet)
        movieSet.setGain(self.gainFile.get())
        movieSet.setDark(self.darkFile.get())
        acq = movieSet.getAcquisition()
        acq.setDoseInitial(self.doseInitial.get())
        acq.setDosePerFrame(self.dosePerFrame.get())

    def _setupFirstImage(self, movie, imgSet):
        # Create a movie object to read dimensions
        dimMovie = movie.clone()
        movieFn = movie.getFileName()

        def decompress(program, args, ext, nExt):
            movieFolder = self._getTmpPath()
            movieName = basename(movie.getFileName())
            movieTmpLink = join(movieFolder, movieName)
            pwutils.cleanPath(movieTmpLink)
            pwutils.createAbsLink(os.path.abspath(movieFn), movieTmpLink)
            self.runJob(program, args % movieName, cwd=movieFolder)
            dimMovie.setFileName(movieTmpLink.replace(ext, nExt))

        if movieFn.endswith('bz2'):
            decompress('bzip2', '-d -f %s', '.bz2', '')

        elif movieFn.endswith('tbz'):
            decompress('tar', 'jxf %s', '.tbz', '.mrc')

        dim = dimMovie.getDim()
        range = [1, dim[2], 1]

        movie.setFramesRange(range)
        imgSet.setDim(dim)
        imgSet.setFramesRange(range)

    def launchSocket(self):
        host = ''  # Where do we get this?!!
        serverSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        serverSocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        serverSocket.bind((host, self.socketPort))
        serverSocket.listen(10)
        serverSocket.setblocking(0)
        self.serverSocket = serverSocket
        self.connectionList = [serverSocket]
        self.info("Socket started on port " + str(self.socketPort))
        return serverSocket

    def iterFilenamesFromSocket(self):
        recv_buffer = 4096  # Advisable to keep it as an exponent of 2
        read_sockets, wr_sockets, err_sockets = select.select(self.connectionList, [], [], 0)
        for sock in read_sockets:
            if sock is self.serverSocket:
                # New connection received through self.serverSocket
                sockfd, addr = self.serverSocket.accept()
                sockfd.setblocking(0)
                self.connectionList.append(sockfd)
                self.debug("Client (%s, %s) connected" % addr)
                continue
            else:
                # Data received from a client, process it
                try:
                    # In Windows, sometimes when a TCP program closes abruptly,
                    # a "Connection reset by peer" exception will be thrown
                    data = sock.recv(recv_buffer)
                    if data:
                        files = shlex.split(data)
                        self.debug("Data received in socket:")
                        self.debug(files)
                        for fileName in files:
                            if os.path.exists(fileName):
                                if fileName in self.importedFiles:
                                    self._spreadMessage('WARNING: Not importing, already imported file %s \n' % fileName,
                                                        sock)
                                    continue
                                else:
                                    self._spreadMessage('OK: Importing file %s \n' % fileName, sock)
                                    fileId = None
                                    yield fileName, fileId
                            else:
                                self._spreadMessage('WARNING: Not importing, path does not exist %s \n' % fileName,
                                                    sock)
                                continue
                    else:
                        continue
                # client disconnected, remove from socket list
                except Exception as e:
                    self.debug("Exception reading socket!!")
                    self.debug(str(e))
                    sock.close()
                    self.connectionList.remove(sock)
                    continue
        return

    def _spreadMessage(self, message, sock):
        try:
            if sock is not None:
                sock.send(message)
            self.debug(message)
        except:
            pass

    def iterNewInputFiles(self):
        """ In the case of importing movies, we want to override this method
        for the case when input are individual frames and we want to create
        movie stacks before importing.
        The frames pattern should contains a part delimited by $.
        The id expression with # is not supported for simplicity.
        """

        if not (self.inputIndividualFrames and self.stackFrames):
            # In this case behave just as
            if self.streamingSocket:
                for fileName, fileId in self.iterFilenamesFromSocket():
                    yield fileName, fileId
            else:
                for fileName, fileId in ProtImportMicBase.iterNewInputFiles(self):
                    yield fileName, fileId
            return

        if self.dataStreaming:
            if self.streamingSocket:
                filePaths = [f[0] for f in self.iterFilenamesFromSocket()]
            else:
                # Consider only the files that are not changed in the fileTime delta
                # if processing data in streaming
                fileTimeout = timedelta(seconds=self.fileTimeout.get())
                filePaths = [f for f in self.getMatchFiles()
                             if not self.fileModified(f, fileTimeout)]
        else:
            filePaths = self.getMatchFiles()

        frameRegex = re.compile("(?P<prefix>.+[^\d]+)(?P<frameid>\d+)")
        # Group all frames for each movie
        # Key of the dictionary will be the common prefix and the value
        # will be a list with all frames in that movie
        frameDict = {}

        for fileName in filePaths:
            fnNoExt = pwutils.removeExt(fileName)

            match = frameRegex.match(fnNoExt)

            if match is None:
                raise Exception("Incorrect match of frame files pattern!")

            d = match.groupdict()
            prefix = d['prefix']
            frameid = int(d['frameid'])

            if prefix not in frameDict:
                frameDict[prefix] = []

            frameDict[prefix].append((frameid, fileName))

        suffix = self.movieSuffix.get()
        ih = ImageHandler()

        for movieFn in self.createdStacks:
            if movieFn not in self.importedFiles:
                yield movieFn, None

        def checkMovie():
            for k, v in frameDict.iteritems():
                movieFn = k + suffix

                if self.writeMoviesInProject:
                    movieFn = self._getExtraPath(os.path.basename(movieFn))

                if (movieFn not in self.importedFiles and
                    movieFn not in self.createdStacks and
                    len(v) == self.numberOfIndividualFrames):
                    movieOut = movieFn

                    if movieOut.endswith("mrc"):
                        movieOut += ":mrcs"

                    # By default we will write the movie stacks
                    # unless we are in continue mode and the file exists
                    writeMovie = True
                    if (self.isContinued() and os.path.exists(movieFn)):
                        self.info("Skipping movie stack: %s, seems to be done"
                                  % movieFn)
                        writeMovie = False

                    if writeMovie:
                        self.info("Writing movie stack: %s" % movieFn)
                        # Remove the output file if exists
                        pwutils.cleanPath(movieFn)

                        for i, frame in enumerate(sorted(v, key=lambda x: x[0])):
                            frameFn = frame[1] # Frame name stored previously
                            ih.convert(frameFn, (i+1, movieOut))

                            if self.deleteFrames:
                                pwutils.cleanPath(frameFn)

                    # Now return the newly created movie file as imported file
                    self.createdStacks.add(movieFn)
                    return
        checkMovie()

    def ignoreCopy(self, source, dest):
        pass

    def getCopyOrLink(self):
        if (self.inputIndividualFrames and self.stackFrames and
            self.writeMoviesInProject):
            return self.ignoreCopy
        else:
            return ProtImportMicBase.getCopyOrLink(self)

    def _cleanUp(self):
        if self.streamingSocket:
            self.debug('Closing socket...')
            self.serverSocket.shutdown(socket.SHUT_RDWR)
            self.serverSocket.close()
