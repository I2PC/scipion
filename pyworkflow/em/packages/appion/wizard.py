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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module implement some wizards
"""

import os
from pyworkflow.em.wizard import EmWizard
from protocol_picking import DogPickerProtPicking
from pyworkflow.em import CoordinatesObjectView
from pyworkflow.em.showj import CLASSIFIER
from pyworkflow.em.packages.xmipp3 import writeSetOfMicrographs
#===============================================================================
# PICKER
#===============================================================================

class DogPickerWizard(EmWizard):
    _targets = [(DogPickerProtPicking, ['diameter', 'threshold'])]
    
    
    def show(self, form):
        autopickProt = form.protocol
        micSet = autopickProt.getInputMicrographs()
        project = autopickProt.getProject()
        micfn = project.getTmpPath(micSet.getName() + '_micrographs.xmd')
        writeSetOfMicrographs(micSet, micfn)
        coordsDir = autopickProt.getProject().getTmpPath()
        # Get current values of the properties
        dogpickerConf = os.path.join(coordsDir, 'dogpicker.conf')
        f = open(dogpickerConf, "w")

        args = {
          "dogpicker" : os.path.join(os.environ['DOGPICKER_HOME'], "ApDogPicker.py"),
          "diameter": autopickProt.diameter,
          "threshold": autopickProt.threshold,
          "apix": micSet.getSamplingRate(),
          "dogpickerOut": os.path.join(coordsDir, 'dogpicker.out')
          }


        f.write("""
        parameters = diameter,threshold
        diameter.value = %(diameter)s
        diameter.label = Diameter
        diameter.help = some help
        threshold.value =  %(threshold)s
        threshold.label = Threshold
        threshold.help = some help
        command = %(dogpicker)s  --thresh=%%(threshold) --diam=%%(diameter) --apix=%(apix)s  --image=%%(micrograph) --outfile=%(dogpickerOut)s 
        """ % args)
        print "Launching picking GUI..."
        CoordinatesObjectView(project, micfn, coordsDir, autopickProt, pickerConf=dogpickerConf).show()
        

