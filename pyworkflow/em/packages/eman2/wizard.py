# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
"""
This module implement some wizards
"""

import os

import pyworkflow as pw
from pyworkflow.em.wizard import EmWizard
from protocol_autopick import SparxGaussianProtPicking
from pyworkflow.em import CoordinatesObjectView
from pyworkflow.em.showj import CLASSIFIER
from pyworkflow.em.packages.xmipp3 import writeSetOfMicrographs
from pyworkflow.utils import makePath, cleanPath
from pyworkflow.utils.utils import readProperties

#===============================================================================
# PICKER
#===============================================================================

class SparxGaussianPickerWizard(EmWizard):
    _targets = [(SparxGaussianProtPicking, ['boxSize', 'lowerThreshold', 'higherThreshold', 'gaussWidth'])]
    
    
    def show(self, form):
        autopickProt = form.protocol
        micSet = autopickProt.getInputMicrographs()
        if not micSet:
            print 'must specify input micrographs'
            return
        project = autopickProt.getProject()
        micfn = micSet.getFileName()
        coordsDir = project.getTmpPath(micSet.getName())
        cleanPath(coordsDir)
        makePath(coordsDir)
        
        pickerProps = os.path.join(coordsDir, 'picker.conf')
        f = open(pickerProps, "w")
        params = ['boxSize', 'lowerThreshold', 'higherThreshold', 'gaussWidth']
        args = {
          "params": ','.join(params),
          "preprocess" : "%s sxprocess.py" % pw.getScipionScript(),
          "picker" : "%s e2boxer.py" % pw.getScipionScript(),
          "convert" : pw.join('apps', 'pw_convert.py'),
          'coordsDir': coordsDir,
          'micsSqlite': micSet.getFileName(),
          "boxSize": autopickProt.boxSize,
          "lowerThreshold": autopickProt.lowerThreshold,
          "higherThreshold": autopickProt.higherThreshold,
          "gaussWidth": autopickProt.gaussWidth,
          "extraParams":autopickProt.extraParams
          }


        f.write("""
        parameters = %(params)s
        boxSize.value = %(boxSize)s
        boxSize.label = Box Size
        boxSize.help = some help
        lowerThreshold.value =  %(lowerThreshold)s
        lowerThreshold.label = Lower Threshold
        lowerThreshold.help = some help
        higherThreshold.help = some help
        higherThreshold.value =  %(higherThreshold)s
        higherThreshold.label = Higher Threshold
        gaussWidth.help = some help
        gaussWidth.value =  %(gaussWidth)s
        gaussWidth.label = Gauss Width
        runDir = %(coordsDir)s
        preprocessCommand = %(preprocess)s demoparms --makedb=thr_low=%%(lowerThreshold):thr_hi=%%(higherThreshold):boxsize=%%(boxSize):gauss_width=%%(gaussWidth):%(extraParams)s
        autopickCommand = %(picker)s --gauss_autoboxer=demoparms --write_dbbox --boxsize=%%(boxSize) --norm=normalize.ramp.normvar %%(micrograph) 
        convertCommand = %(convert)s --coordinates --from eman2 --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s
        """ % args)
        f.close()
        process = CoordinatesObjectView(project, micfn, coordsDir, autopickProt,
                                        pickerProps=pickerProps).show()
        process.wait()
        myprops = readProperties(pickerProps)

        if myprops['applyChanges'] == 'true':
            for param in params:
                form.setVar(param, myprops[param + '.value'])

