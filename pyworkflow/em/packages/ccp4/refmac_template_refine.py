
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
template_refine="""#!/bin/sh
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

# 3D map
MAPFILE=%(MAPFILE)s

# PDB molecule without extension or path just the basename
MOL=%(PDBFILE)s

PDBFILE=${MOL}.pdb

# Directory of output files (extra folder)
OUTPUTDIR=%(OUTPUTDIR)s

# CCP4 PATH
PATHCCP4=%(CCP4_HOME)s

PATHMRCBIN=$PATHCCP4/bin
PATHMRCENV=$PATHMRCBIN/ccp4.setup-sh
. $PATHMRCENV

#refmac binary
refmac=%(REFMAC_BIN)s

# Delete some temporary files. If not and the script is executed
# two times there will be conflicts

# refine the structure

$refmac HKLIN ${OUTPUTDIR}map2mtz.mtz \\
        XYZIN  ${OUTPUTDIR}pdbset.pdb \\
        XYZOUT ${OUTPUTDIR}refmac-refined.pdb\\
        atomsf ${PATHCCP4}/lib/data/atomsf_electron.lib \\
        > ${OUTPUTDIR}refine.log <<EOF  

LABIN FP=Fout0 PHIB=Pout0

# set resolution limits. The FSC=0.143 cutoff is recommended. 
RESO = %(RESOMIN)f  %(RESOMAX)f

# set B-factors at prior to refinement:
BFACTOR_SET=%(BFACTOR_SET)s

# specify number of refinement cycles:
NCYCLE = %(NCYCLE)d

# set refinement weight:
WEIGHT MATRIX = %(WEIGHT MATRIX)s

#specify any other keyword:
##MAKE CISP NO
##MAKE SS NO
MAKE HYDR NO #remove to use hydrogens
SOLVENT NO

source EM MB

###specify external restraints below###

ridge dist sigma 0.01
ridge dist dmax 4.2

END
EOF
          echo Output ${OUTPUTDIR}refmac-refined.pdb 

      #fi
#done
"""