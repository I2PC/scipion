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
"""
Some Spider protocol base classes.
"""

from pyworkflow.em import EMProtocol, ProtClassify2D
from pyworkflow.utils import removeExt, replaceBaseExt, runJob

from ..spider import loadEnvironment, getScript, END_HEADER, REGEX_KEYFRL, REGEX_KEYVALUE, SPIDER
from ..convert import writeSetOfImages



class SpiderProtocol(EMProtocol):
    """ Sub-class of EMProtocol to group some common Spider utils. """
            
    def convertInput(self, attrName, stackFn, selFn):
        """ Convert from an input pointer of SetOfImages to Spider.
        Params:
            attrName: the attribute name of the input pointer
            stackFn: the name of the stack for converted images
            selFn: the name of the selection file.
        """
        imgSetPointer = getattr(self, attrName)
        writeSetOfImages(imgSetPointer.get(), stackFn, selFn)
        
    def _getFileName(self, key):
        """ Give a key, append the extension
        and prefix the protocol working dir. 
        """
        template = '%(' + key + ')s.%(ext)s'
        
        return self._getPath(template % self._params)
    
    def getExt(self):
        """ Return the extension used in the script,
        stored in a dictionary called self._params. 
        """
        return self._params['ext']
    
    def __substituteVar(self, match, paramsDict, lineTemplate):
        if match and match.groupdict()['var'] in paramsDict:
            d = match.groupdict()
            d['value'] = paramsDict[d['var']]
            return lineTemplate % d
        return None
    
    def runScript(self, inputScript, ext, paramsDict):
        """ This function will create a valid Spider script
        by copying the template and replacing the values in dictionary.
        After the new file is read, the Spider interpreter is invoked.
        """
        loadEnvironment()
        self._enterWorkingDir()
        
        outputScript = replaceBaseExt(inputScript, ext)
    
        fIn = open(getScript(inputScript), 'r')
        fOut = open(outputScript, 'w')
        inHeader = True # After the end of header, not more value replacement
        inFrL = False
        
        for i, line in enumerate(fIn):
            if END_HEADER in line:
                inHeader = False
            if inHeader:
                try:
                    newLine = self.__substituteVar(REGEX_KEYVALUE.match(line), paramsDict, 
                                              "%(var)s%(s1)s=%(s2)s%(value)s%(rest)s\n")
                    if newLine is None and inFrL:
                        newLine = self.__substituteVar(REGEX_KEYFRL.match(line), paramsDict, 
                                                  "%(var)s%(value)s%(rest)s\n")
                    if newLine:
                        line = newLine
                except Exception, ex:
                    print ex, "on line (%d): %s" % (i+1, line)
                    raise ex
                inFrL = line.lower().startswith("fr l")
            fOut.write(line)
        fIn.close()
        fOut.close()    
    
        scriptName = removeExt(outputScript)
        log = getattr(self, '_log', None)        
        runJob(log, SPIDER, "%(ext)s @%(scriptName)s" % locals())
        
        self._leaveWorkingDir()
    
    
class SpiderProtClassify(ProtClassify2D, SpiderProtocol):
    """ Diday's method, using 'CL CLA' 
    """
    def __init__(self, **kwargs):
        ProtClassify2D.__init__(self, **kwargs)
        SpiderProtocol.__init__(self, **kwargs)
        
    def getClassDir(self):
        return self._params['classDir']