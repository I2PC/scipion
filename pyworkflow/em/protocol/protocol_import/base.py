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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
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

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        importChoices = self._getImportChoices()
        filesCondition = self._getFilesCondition()
        
        form.addSection(label='Import')

        if len(importChoices) > 1: # not only from files
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
                           "The path can also contain wildcards to select\n"
                           "from several folders.\n\n"
                           "For example:\n"
                           "  ~/Particles/\n"
                           "  data/day??_micrographs/")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern', 
                      condition=filesCondition,
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.")

        form.addParam('copyFiles', params.BooleanParam, default=False, 
                      expertLevel=params.LEVEL_ADVANCED, 
                      label="Copy files?",
                      help="By default the files are not copied into the\n"
                           "project to avoid data duplication and to save\n"
                           "disk space. Instead of copying, symbolic links are\n"
                           "created pointing to original files. This approach\n"
                           "has the drawback that if the project is moved to\n"
                           "another computer, the links need to be restored.\n")

        self._defineImportParams(form)

        self._defineAcquisitionParams(form)

        form.addSection('Streaming')
        
        form.addParam('dataStreaming', params.BooleanParam, default=False,
              label="Process data in streaming?",
              help="Select this option if you want import data as it is\n"
                   "generated and process on the fly by next protocols."
                   "In this case the protocol will keep running to check \n"
                   "new files and will update the output Set, which can \n"
                   "be used right away by next steps.\n")
        
        form.addParam('timeout', params.IntParam, default=7200,
              condition='dataStreaming',
              label="Timeout (secs)",
              help="Interval of time (in seconds) after which, if no new \n"
                   "file is detected, the protocol will ends.\n"
                   "When finished, the output Set will be closed and\n"
                   "not more data will be added to it. \n")

        form.addParam('fileTimeout', params.IntParam, default=30,
              condition='dataStreaming',
              label="File timeout (secs)",
              help="Interval of time (in seconds) after which, if a file \n"
                   "have not changed, we consider as a new file. \n")
        
        form.addParam('endTokenFile', params.StringParam, default=None,
              condition='dataStreaming',
              label="End token file",
              help="Specify an ending file if you want to have more control\n"
                   "about when to stop the import of files.")

    def _defineImportParams(self, form):
        """ Override to add options related to the different types
        of import that are allowed by each protocol.
        """
        pass

    def _defineAcquisitionParams(self, form):
        """ Override to add options related to acquisition info.
        """
        pass

    def _getDefaultChoice(self):
        return  self.IMPORT_FROM_FILES
    
    #--------------------------- INFO functions --------------------------------
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
            
        return errors
    
    #--------------------------- BASE methods to be overriden ------------------
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
    
    #--------------------------- UTILS functions -------------------------------
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
            self._idRegex = re.compile(pattern.replace(g, '(%s)' % ('\d'*n)))
            pattern = pattern.replace(g, '[0-9]'*n)
        
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

    def iterFiles(self):
        """ Iterate through the files matched with the pattern.
        Provide the fileName and fileId.
        """
        filePaths = self.getMatchFiles()
        
        for fileName in filePaths:
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


