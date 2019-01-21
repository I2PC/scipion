# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os
from os.path import join
from glob import glob
import re
from datetime import timedelta, datetime

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.utils.path import expandPattern, copyFile, createAbsLink
from pyworkflow.em.protocol import EMProtocol

class ProtImport(EMProtocol):
    """ Base class for other all Import protocols. """


class ProtImportFiles(ProtImport):
    """ Base class for other Import protocols. 
    All imports protocols will have:
    1) Several options to import from (_getImportOptions function)
    2) First option will always be "from files". (for this option 
      files with a given pattern will be retrieved  and the ### will 
      be used to mark an ID part from the filename.
      - For each file a function to process it will be called
        (_importFile(fileName, fileId))
    """
    IMPORT_FROM_FILES = 0
    BLACKLIST_REGEXPS = 0

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        importChoices = self._getImportChoices()
        filesCondition = self._getFilesCondition()

        form.addSection(label='Import')

        if len(importChoices) > 1:  # not only from files
            form.addParam('importFrom', params.EnumParam,
                          choices=importChoices, default=self._getDefaultChoice(),
                          label='Import from',
                          help='Select the type of import.')
        else:
            form.addHidden('importFrom', params.EnumParam,
                           choices=importChoices, default=self.IMPORT_FROM_FILES,
                           label='Import from',
                           help='Select the type of import.')
        form.addParam('filesPath', params.PathParam,
                      condition=filesCondition,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/Particles/data/day??_micrographs/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/Particles/data/day*_micrographs/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/Particles/data/day#_micrographs/\n"
                           "'#' represents one digit that will be used as "
                           "micrograph ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern',
                      condition=filesCondition,
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.")
        form.addParam('copyFiles', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Copy files?",
                      help="By default the files are not copied into the "
                           "project to avoid data duplication and to save "
                           "disk space. Instead of copying, symbolic links are "
                           "created pointing to original files. This approach "
                           "has the drawback that if the project is moved to "
                           "another computer, the links need to be restored.")

        self._defineImportParams(form)

        self._defineAcquisitionParams(form)

        form.addSection('Streaming')

        form.addParam('dataStreaming', params.BooleanParam, default=False,
                      label="Process data in streaming?",
                      help="Select this option if you want import data as it is "
                           "generated and process on the fly by next protocols. "
                           "In this case the protocol will keep running to check "
                           "new files and will update the output Set, which can "
                           "be used right away by next steps.")

        form.addParam('timeout', params.IntParam, default=43200,
                      condition='dataStreaming',
                      label="Timeout (secs)",
                      help="Interval of time (in seconds) after which, if no new file "
                           "is detected, the protocol will end. When finished, "
                           "the output Set will be closed and no more data will be "
                           "added to it. \n"
                           "Note 1:  The default value is  high (12 hours) to avoid "
                           "the protocol finishes during the aqcuisition of the "
                           "microscpe. You can also stop it from right click and press "
                           "STOP_STREAMING.\n"
                           "Note 2: If you're using individual frames when importing "
                           "movies, the timeout won't be refreshed until a whole "
                           "movie is stacked.")

        form.addParam('fileTimeout', params.IntParam, default=30,
                      condition='dataStreaming',
                      label="File timeout (secs)",
                      help="Interval of time (in seconds) after which, if a file has "
                           "not changed, we consider it as a new file. \n")

        self._defineBlacklistParams(form)

    def _defineImportParams(self, form):
        """ Override to add options related to the different types
        of import that are allowed by each protocol.
        """
        pass

    def _getBlacklistSetClass(self):
        """ Returns the class to be blacklisted by this protocol.
        """
        return "SetOfImages"

    def _defineBlacklistParams(self, form):
        """ Options to blacklist certain items when launching the
        import protocol.
        """
        form.addSection(label="Blacklist")
        form.addParam("blacklistSet", params.PointerParam,
                      pointerClass=self._getBlacklistSetClass(),
                      allowsNull=True,
                      label="Blacklist Set",
                      help="Files on this set will not be imported")
        form.addParam('blacklistDateFrom', params.StringParam,
                      label="Blacklist from date",
                      allowsNull=True,
                      help="Files acquired after this date will not be imported. "
                           "Must follow format: YYYY-mm-dd HH:MM:SS \n"
                           "e.g: 2019-01-14 14:18:05")
        form.addParam('blacklistDateTo', params.StringParam,
                      label="Blacklist to date",
                      allowsNull=True,
                      help="Files acquired before this date will not be imported. "
                           "Must follow format: YYYY-mm-dd HH:MM:SS \n"
                           "e.g: 2019-01-14 14:18:05")
        form.addParam('useRegexps', params.EnumParam,
                      default=self.BLACKLIST_REGEXPS,
                      choices=['RegExps', 'File names'],
                      label='Blacklist file type',
                      help="Choose RegExp if the black list file contains regular expressions. Set to File Names if "
                           "the black list file contains file names.")
        form.addParam('blacklistFile', params.FileParam,
                      label="Blacklist File",
                      allowsNull=True,
                      help="Blacklist everything included in this file. If Use RegExps is True,"
                           "lines will be interpreted as regular expressions. E.g: \n"
                           "(.*)GRID_0[1-5](.*)\n"
                           "(.*)/GRID_10/Falcon_2019_01_14-16_(.*)\n"
                           "If Use RegExps is False, lines will be interpreted as file names. E.g.\n"
                           "/path/to/GRID_10/Falcon_2019_01_14-16_51_20_0_movie.mrcs\n"
                           "/path/to/GRID_10/Falcon_2019_01_14-16_55_40_0_movie.mrcs"
                      )

    def _defineAcquisitionParams(self, form):
        """ Override to add options related to acquisition info.
        """
        pass

    def _getDefaultChoice(self):
        return self.IMPORT_FROM_FILES

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        if self.importFrom == self.IMPORT_FROM_FILES:
            if not self.getPattern():
                errors.append("The path and pattern can not be both empty!!!")
            else:
                # Just check the number of files matching the pattern
                self.getMatchFiles()
                if self.numberOfFiles == 0:
                    errors.append("There are no files matching the pattern %s"
                                  % self.getPattern())

        dates = [self.blacklistDateFrom.get(), self.blacklistDateTo.get()]
        parsedDates = []
        for d in dates:
            if d:
                try:
                    parsedDates.append(datetime.strptime(self.blacklistDateTo.get(), "%Y-%m-%d %H:%M:%S"))
                except ValueError as e:
                    errors.append("Bad date formatting in blacklist date %s: %s" % (d, e))

        if len(parsedDates) == 2 and parsedDates[0] > parsedDates[1]:
            errors.append("Wrong blacklist dates: date from must be earlier than date to")
        return errors

    # --------------------------- BASE methods to be overriden ------------------
    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        return ['files']

    def _getFilesCondition(self):
        """ Return an string representing the condition
        when to display the files path and pattern to grab
        files.
        """
        return '(importFrom == %d)' % self.IMPORT_FROM_FILES

    # --------------------------- UTILS functions -------------------------------
    def getPattern(self):
        """ Expand the pattern using environ vars or username
        and also replacing special character # by digit matching.
        """
        self._idRegex = None
        filesPath = self.filesPath.get('').strip()
        filesPattern = self.filesPattern.get('').strip()

        if filesPattern:
            fullPattern = join(filesPath, filesPattern)
        else:
            fullPattern = filesPath

        pattern = expandPattern(fullPattern.replace("$", ""))
        match = re.match('[^#]*(#+)[^#]*', pattern)

        if match is not None:
            g = match.group(1)
            n = len(g)
            # prepare regex pattern - place ids, handle *, handle ?
            idregex = pattern.replace(g, '(%s)' % ('[0-9]' * n))
            idregex = idregex.replace('*', '.*')
            idregex = idregex.replace('?', '.')
            self._idRegex = re.compile(idregex)
            pattern = pattern.replace(g, '[0-9]' * n)

        return pattern

    def getMatchFiles(self, pattern=None):
        """ Return a sorted list with the paths of files that
        matched the pattern.
        """
        if pattern is None:
            pattern = self.getPattern()

        filePaths = glob(pattern)
        filePaths.sort()
        self.numberOfFiles = len(filePaths)

        return filePaths

    def getCopyOrLink(self):
        # Set a function to copyFile or createLink
        # depending in the user selected option 
        if self.copyFiles:
            return copyFile
        else:
            return createAbsLink

    def fileModified(self, fileName, fileTimeout):
        """ Check if the fileName modification time is less
        than a given timeout.
        Params:
            fileName: input filename that will be checked.
            fileTimeout: timeout """
        self.debug('Checking file: %s' % fileName)
        mTime = datetime.fromtimestamp(os.path.getmtime(fileName))
        delta = datetime.now() - mTime
        self.debug('   Modification time: %s' % pwutils.prettyTime(mTime))
        self.debug('   Delta: %s' % pwutils.prettyDelta(delta))

        return delta < fileTimeout

    def _getUniqueFileName(self, fileName):
        """To be overwritten by subclasses"""
        return os.path.basename(fileName)

    def getItemsToBlacklistFromFile(self):
        if not hasattr(self, '_fileItemsToBlacklist'):
            blacklistfile = self.blacklistFile.get()
            blacklistItems = set()
            if blacklistfile:
                with open(blacklistfile, 'r') as f:
                    for blacklistedItem in f:
                        blacklistedItem = blacklistedItem.strip()
                        blacklistItems.add(blacklistedItem)
            self._fileItemsToBlacklist = blacklistItems

        return self._fileItemsToBlacklist

    def getBlacklistedItems(self):
        if not hasattr(self, '_blacklistedItems'):
            self._blacklistedItems = set()
        return self._blacklistedItems

    def isBlacklisted(self, fileName):
        # check if already blacklisted
        blacklistedItems = self.getBlacklistedItems()
        if fileName in blacklistedItems:
            return True

        # Blacklisted by set
        blacklistSet = self.blacklistSet.get()
        if blacklistSet is not None:
            for img in blacklistSet:
                blacklistFileName = img.getFileName()
                if ((os.path.islink(blacklistFileName)
                     and fileName == os.readlink(blacklistFileName))
                        or (self._getUniqueFileName(fileName) == os.path.basename(blacklistFileName))):
                    print("Blacklist warning: %s is blacklisted by the input set" % fileName)
                    blacklistedItems.add(fileName)
                    return True

        # Blacklisted by date
        blacklistDateFrom = self.blacklistDateFrom.get()
        blacklistDateTo = self.blacklistDateTo.get()
        doDateBlacklist = blacklistDateFrom is not None and blacklistDateTo is not None
        if doDateBlacklist:
            fileDate = datetime.fromtimestamp(os.path.getmtime(fileName))
            if blacklistDateFrom:
                parsedDateFrom = datetime.strptime(blacklistDateFrom, "%Y-%m-%d %H:%M:%S")
                if blacklistDateTo:
                    parsedDateTo = datetime.strptime(blacklistDateTo, "%Y-%m-%d %H:%M:%S")
                    if parsedDateFrom <= fileDate <= parsedDateTo:
                        print("Blacklist warning: %s is blacklisted by date" % fileName)
                        blacklistedItems.add(fileName)
                        return True
                else:
                    if parsedDateFrom <= fileDate:
                        print("Blacklist warning: %s is blacklisted by date" % fileName)
                        blacklistedItems.add(fileName)
                        return True

            elif blacklistDateTo:
                parsedDateTo = datetime.strptime(blacklistDateTo, "%Y-%m-%d %H:%M:%S")
                if fileDate <= parsedDateTo:
                    print("Blacklist warning: %s is blacklisted by date" % fileName)
                    blacklistedItems.add(fileName)
                    return True

        # Blacklisted by file
        items2blacklist = self.getItemsToBlacklistFromFile()
        for item2blacklist in items2blacklist:
            if self.useRegexps.get() == self.BLACKLIST_REGEXPS:
                if re.match(item2blacklist, fileName):
                    print("Blacklist warning: %s matched blacklist regexp %s"
                          % (fileName, item2blacklist))
                    blacklistedItems.add(fileName)
                    return True
            elif fileName in item2blacklist:
                print("Blacklist warning: %s is blacklisted " % fileName)
                blacklistedItems.add(fileName)
                return True

    def iterFiles(self):
        """ Iterate through the files matched with the pattern.
        Provide the fileName and fileId.
        """
        filePaths = self.getMatchFiles()

        for fileName in filePaths:
            if self.isBlacklisted(fileName):
                continue
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise Exception("File '%s' doesn't match the pattern '%s'"
                                    % (fileName, self.getPattern()))

                fileId = int(match.group(1))

            else:
                fileId = None

            yield fileName, fileId
