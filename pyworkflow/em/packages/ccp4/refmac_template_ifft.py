
# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
#                Marta Martinez (mmmtnez@cnb.csic.es)
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
template_ifft="""#!/bin/sh
# This  script will run REFMAC for a model against an EM map. 
# # 
# Prerequistes:
#	A model in PDB or mmcif format. The CRYST card must agree with the box size of the map.
#   The box size can be obtained using MAPDUMP from CCP4.
#	A map, or maps, in .mrc or .map format.
#	External restraint files.
#	Electron scattering file: atomsf_electron.lib.
#
# The script starts by calculating complex structure factors around a given radius taken from the input model.
# The radius is given in Angstrom and is can be user-defined (SFCALC MRAD 3). 
#

# set generateMaskedVolume to True if you want to see the mask using chimera
generateMaskedVolume=%(MASKED_VOLUME)s

#Directory of output files (extra folder)
OUTPUTDIR=%(OUTPUTDIR)s

#CCP4 PATH
PATHCCP4=%(CCP4_HOME)s

PATHMRCENV=$PATHCCP4/setup-scripts/ccp4.setup-sh
PATHMRCBIN=$PATHCCP4/bin
. $PATHMRCENV

#transform the mask (originally in fourier space) so we can see it (in real space)
#this step is not needed but allow users to check that the PDB file and the 3Dmap are
#in the same coordinate system


fft HKLIN ${OUTPUTDIR}_masked_fs.mtz MAPOUT ${OUTPUTDIR}masked_fs.map << END-fft > ${OUTPUTDIR}ifft.log
    LABIN F1=Fout0 PHI=Pout0 
    SCALE F1 1.0 300.0 
    RESOLUTION %(RESOMAX)d
    GRID %(XDIM)d %(YDIM)d %(ZDIM)d  
    END 
END-fft

#done
"""