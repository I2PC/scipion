
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

# set generateMaskedVolume to True if you want to see the mask using chimera
generateMaskedVolume=%(MASKED_VOLUME)s

# PDB molecule without extension or path just the basename
MOL=%(PDBFILE)s

PDBFILE=${MOL}.pdb

# Directory of output files (extra folder)
OUTPUTDIR=%(OUTPUTDIR)s

# CCP4 PATH
PATHCCP4=%(CCP4_HOME)s

PATHMRCENV=$PATHCCP4/setup-scripts/ccp4.setup-sh
PATHMRCBIN=$PATHCCP4/bin
. $PATHMRCENV

#refmac binary
refmac=%(REFMAC_BIN)s

# suffix for output files 
molecule_id=${MOL}
# Delete some temporary files. If not and the script is executed
# two times there will be conflicts

# refine the structure

$refmac HKLIN ${OUTPUTDIR}_masked_fs.mtz \\
        HKLOUT ${OUTPUTDIR}$molecule_id-refined.mtz \\
        XYZIN  ${OUTPUTDIR}$molecule_id-initial.pdb \\
        XYZOUT ${OUTPUTDIR}$molecule_id-refined.pdb\\
        atomsf ${PATHCCP4}/bin/atomsf_electron.lib \\
        > ${OUTPUTDIR}refine.log <<EOF  

LABIN FP=Fout0 PHIB=Pout0

# set resolution limits. The FSC=0.143 cutoff is recommended. 
RESO = (%(RESOMIN)f  %(RESOMAX)f)

# set B-factors at prior to refinement:
BFACTOR_SET=%(BFACTOR_SET)s

# specify map sharpening to be used during refinement:
REFI_SHARPEN=%(REFI_SHARPEN)s

# specify number of refinement cycles:
NCYCLE = %(NCYCLE)d

# set refinement weight:
WEIGHT MATRIX = %(WEIGHT MATRIX)f

#specify any other keyword:
##MAKE CISP NO
##MAKE SS NO
#MAKE HYDR N #remove to use hydrogens

source em 
#pdbout format mmcif #remove to output mmcif
@shifts.txt


###specify external restraints below###

EXTERNAL USE MAIN
EXTERNAL DMAX 4.2
EXTERNAL WEIGHT SCALE 5
EXTERNAL WEIGHT GMWT 0.1
@external.restraints

MONI DIST 1000000
END
EOF
          echo Output ${OUTPUTDIR}$molecule_id-refined.pdb  
          #fi
          #cd - 
      #fi
#done
"""