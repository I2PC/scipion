# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
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
Since the Projection Matching protocol of Xmipp 3 has a very large
form definition, we have separated in this sub-module.
"""

# Functions outside th loop loop for xmipp_projection_matching
def insertExecuteCtfGroups(self):
    #...
    self._insertRunJobStep('xmipp_ctf_group') #...

def insertInitAngularReferenceFile(self):
    #...
    self._insertRunJobStep('') #...

# Functions in loop for xmipp_projection_matching

def insertAngularProjectLibrary(self):
    #...
    self._insertRunJobStep('xmipp_angular_project_library') #...

def insertProjectionMatching(self):
    #...
    self._insertRunJobStep('xmipp_angular_projection_matching') #...

def insertAssignImagesToReferences(self):
    #...
    self._insertRunJobStep('') #...

def insertAngularClassAverage(self):
    #...
    self._insertRunJobStep('xmipp_angular_class_average') #...

def insertReconstruction(self):
    #...
    self._insertRunJobStep('xmipp_reconstruct_fourier') #...

def insertComputeResolution(self):
    #...
    self._insertRunJobStep('xmipp_resolution_fsc') #...

def insertFilterVolume(self):
    #...
    self._insertRunJobStep('xmipp_transform_filter') #...

