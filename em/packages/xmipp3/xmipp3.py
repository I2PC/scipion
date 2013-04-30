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
This sub-package will contains Xmipp3.0 specific protocols
"""

import os

def getXmippPath(*paths):
    '''Return the path the the Xmipp installation folder
    if a subfolder is provided, will be concatenated to the path'''
    if os.environ.has_key('XMIPP_HOME'):
        return os.path.join(os.environ['XMIPP_HOME'], *paths)  
    else:
        raise Exception('XMIPP_HOME environment variable not set')
    
    
class XmippProtocol():
    """ This class groups some common functionalities that
    share some Xmipp protocols, like converting steps.
    """
    
    def convertInputToXmipp(self, inputName, convertFunc, resultFn):
        """ This step can be used whenever a convertion is needed.
        It will receive the inputName and get this attribute from self,
        invoke the convert function and check the result files if
        convertion was done (otherwise the input was already in Xmipp format).
        """
        inputAttr = getattr(self, inputName)
        inputXmipp = convertFunc(inputAttr, resultFn)
        
        if inputXmipp != inputAttr:
            self._insertChild(inputName + 'Xmipp', inputXmipp)
            return [resultFn] # validate resultFn was produced if converted
        
    def getConvertedInput(self, inputName):
        """ Retrieve the converted input, it can be the case that
        it is the same as input, when not convertion was done. 
        """
        return getattr(self, inputName + 'Xmipp', getattr(self, inputName))
        

