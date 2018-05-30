
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
template_map_mtz="""#!/bin/sh
# This  script will run REFMAC for a model against an EM map. 
# # 
# Prerequistes:
#	A model in PDB or mmcif format. The CRYST card must agree with the box size 
#   of the map.
#   The box size can be obtained using MAPDUMP from CCP4.
#	A map, or maps, in .mrc or .map format.
#	External restraint files.
#	Electron scattering file: atomsf_electron.lib.
#
# The script starts by calculating complex structure factors around a given
# radius taken from the input model.
# The radius is given in Angstrom and is can be user-defined (SFCALC MRAD 3). 
#

#refmac binary
refmac=%(REFMAC_BIN)s

#pdbset binary
pdbset=%(PDBSET_BIN)s

#PATH to PDB file
PDBDIR=%(PDBDIR)s

#PDB molecule without extension or path just the basename
MOL=%(PDBFILE)s

PDBFILE=${MOL}.pdb

#3D MAP FILENAME
MAPFILE=%(MAPFILE)s

#CCP4 PATH
PATHCCP4=%(CCP4_HOME)s

#Directory of output files (extra folder)
OUTPUTDIR=%(OUTPUTDIR)s

# Delete some temporary files. Otherwise if the script is executed
# two times there will be conflicts
RM='rm -f'
${RM} ${OUTPUTDIR}map2mtz.mtz ${OUTPUTDIR}map2mtz.log 
${RM} ${OUTPUTDIR}_orig_data_start.txt
${RM} ${OUTPUTDIR}pdbset.pdb
${RM} ${OUTPUTDIR}refmac-refined.mtz ${OUTPUTDIR}refmac-refined.pdb

#################################################################
# Nothing to be changed below (unless you know what you are doing)
##################################################################
#load ccp4 enviroment
###TODO: MOVE THIS INITIALIZATION TO convert.getEnv
PATHMRCBIN=$PATHCCP4/bin
PATHMRCENV=$PATHMRCBIN/ccp4.setup-sh
. $PATHMRCENV

# create a mask by calculating complex structure factors around a given radius 
# taken from the input model 
pdb_in=${PDBDIR}/${PDBFILE}

echo we will work on $pdb_in
$refmac MAPIN ${MAPFILE} \\
    HKLOUT ${OUTPUTDIR}map2mtz.mtz \\
    > ${OUTPUTDIR}map2mtz.log \\
    << eof
    MODE SFCALC 
    RESO %(RESOMAX)f
    SOURCE EM MB
    END
eof

$pdbset XYZIN  $pdb_in \\
    XYZOUT ${OUTPUTDIR}pdbset.pdb  > ${OUTPUTDIR}pdbset.log \\
    << eof
    CELL %(Xlength)f %(Ylength)f %(Zlength)f  90.0 90.0 90.0
    END
eof

#done
"""

