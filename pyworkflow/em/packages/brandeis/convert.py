# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
# *
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
This module contains converter functions that will serve to:
1. Write from base classes to Grigorieff packages specific files
2. Read from Grigo packs files to base classes
"""

from brandeis import *



def readSetOfClasses3D(classes3DSet, fileparList, volumeList):
    """read from frealign .par.
    """
    imgSet = classes3DSet.getImages()
    
    for ref, volFn in enumerate(volumeList):
        class3D = Class3D()
        class3D.setObjId(ref+1)
        vol = Volume()
        vol.copyObjId(class3D)
        vol.setLocation(volFn)
        
        class3D.setRepresentative(vol)
        classes3DSet.append(class3D)
        
        file1 = fileparList[ref]
        f1 = open(file1)
        for l in f1:
            if not l.startswith('C'):
                values = l.split()
                prob = float(values[11])
                if prob > 0:
                    objId = int(values[7])
                    img = imgSet[objId]
                    class3D.append(img)
        f1.close()
        
        # Check if write function is necessary
        class3D.write()
        

def parseCtffindOutput(filename):
    """ Retrieve defocus U, V and angle from the 
    output file of the ctffind3 execution.
    """
    f = open(filename)
    result = None
    for line in f:
        if 'Final Values' in line:
            # Take DefocusU, DefocusV and Angle as a tuple
            # that are the first three values in the line
            result = tuple(map(float, line.split()[:3]))
            break
    f.close()
    return result
