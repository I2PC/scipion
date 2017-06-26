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
from pyworkflow.em import CoordinatesObjectView
from pyworkflow.utils import makePath, cleanPath, readProperties

from protocol_dogpicker import DogPickerProtPicking


#===============================================================================
# PICKER
#===============================================================================

class DogPickerWizard(EmWizard):
    _targets = [(DogPickerProtPicking, ['diameter', 'threshold'])]


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
        # Get current values of the properties
#         micfn = os.path.join(coordsDir, 'micrographs.xmd')
#         writeSetOfMicrographs(micSet, micfn)
        dogpickerProps = os.path.join(coordsDir, 'picker.conf')
        f = open(dogpickerProps, "w")

        args = {
          "dogpicker" : os.path.join(os.environ['DOGPICKER_HOME'], "ApDogPicker.py"),
          "convert" : pw.join('apps', 'pw_convert.py'),
          'coordsDir': coordsDir,
          'micsSqlite': micSet.getFileName(),
          "diameter": autopickProt.diameter,
          "threshold": autopickProt.threshold,
          "apix": micSet.getSamplingRate()
          }


        f.write("""
        parameters = diameter,threshold
        diameter.value = %(diameter)s
        diameter.label = Diameter
        diameter.help = some help
        threshold.value =  %(threshold)s
        threshold.label = Threshold
        threshold.help = some help
        autopickCommand = %(dogpicker)s  --thresh=%%(threshold) --diam=%%(diameter) --apix=%(apix)s  --image=%%(micrograph) --outfile=%(coordsDir)s/%%(micrographName).txt 
        convertCommand = %(convert)s --coordinates --from dogpicker --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s
        """ % args)
        f.close()
        process = CoordinatesObjectView(project, micfn, coordsDir, autopickProt, pickerProps=dogpickerProps).show()
        process.wait()
        # Check if the wizard changes were accepted or just canceled

        myprops = readProperties(dogpickerProps)

        if myprops['applyChanges'] == 'true':
            form.setVar('diameter', myprops['diameter.value'])
            form.setVar('threshold', myprops['threshold.value'])

