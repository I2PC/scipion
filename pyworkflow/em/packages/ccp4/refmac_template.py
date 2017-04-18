
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
Template=
"""

#
# This script will run REFMAC for a model against an EM map. 
# In this example the model is A/molecule-A.pdb
# 
# The latest version of REFMAC is required from CCP4 (www.ccp4.ac.uk)
# 
# Prerequistes:
#	A model in PDB or mmcif format. The CRYST card must agree with the box size of the map. The box size can be obtained using MAPDUMP from CCP4.
#	A map, or maps, in .mrc or .map format.
#	External restraint files.
#	Electron scattering file: atomsf_electron.lib.
#
# The script starts by calculating complex structure factors around a given radius taken from the input model.
# The radius is given in Angstrom and is can be user-defined (SFCALC MRAD 3). 
#

REFMACPATH=/usr/local/CCP4/
# set generateMaskedVolume to True if you want to see the mask using chimera
generateMaskedVolume=true

#refmac binary
refmac=${REFMACPATH}/refmacgfortran

#PATH to PDB file
PDBDIR=./

#PDB molecule without extension
MOL= %(mol)s
PDBFILE=${MOL}.pdb

#3D MAP FILENAME ###(extension .mrc? always)
MAPFILE= %(vol)s

#CCP4 PATH
PATHCCP4=/usr/local/CCP4/ccp4-6.5

#CHIMERA command
CHIMERA=/usr/local/bin/chimera

#suffix for output files 
molecule_id=${MOL}
# Delete some temporaly files. If not and the script is executed
# two times there will be conflicts
if [ "$generateMaskedVolume" = true ] ; then
    rm  average_for_refmac.mtz refmac-makes-sfs.log orig_data_rescut.txt
    rm $molecule_id-initial.pdb masked_fs.mtz shifts.txt
    echo "GENERATE MASKED VOLUME"
else
    echo "DO NOT GENERATE MASKED VOLUME"
fi
rm $molecule_id-refined.mtz $molecule_id-refined.pdb
rm $molecule_id-refmac-refines.log ifft.log
rm  masked_fs.map

#just show what is in the directory
ls

#################################################################
#Nothing to be changed below (unless you know what you are doing)
#load ccp4 enviroment
PATHMRCENV=$PATHCCP4/setup-scripts/ccp4.setup-sh
PATHMRCBIN=$PATHCCP4/bin

echo $PATHMRCENV
#####echo "source must be done outside the script"
source $PATHMRCENV

#for molecule_id in [A]######
#do
     pdb_in=${PDBDIR}/${PDBFILE}
     if [ "$generateMaskedVolume" = true ] ; then
      #if [ -e ${pdb_in} ] ; then
          echo we will work on $pdb_in
          #cd $molecule_id
          #if [ ! -e  refmac-makes-sfs.log ] ; then 
   	     $refmac MAPIN ${MAPFILE} \
                     HKLOUT average_for_refmac.mtz \
		     XYZIN $pdb_in \
		     XYZOUT $molecule_id-initial.pdb \
                     > refmac-makes-sfs.log \
            << eof
               MODE SFCALC 
               sfcalc mapradius 3
	       SFCALC mradius 3
               sfcalc shift
               END
eof
     fi    
#check volume

fft HKLIN masked_fs.mtz MAPOUT masked_fs.map << END-fft > ifft.log 
LABIN F1=Fout0 PHI=Pout0 
SCALE F1 1.0 300.0 
RESOLUTION 3.0 
GRID 200 200 200 
END 
END-fft
${CHIMERA} masked_fs.map &

#end check
                $refmac HKLIN masked_fs.mtz \
                        HKLOUT $molecule_id-refined.mtz \
                        XYZIN  $molecule_id-initial.pdb \
                        XYZOUT $molecule_id-refined.pdb\
			atomsf ${PATHCCP4}/bin/atomsf_electron.lib \
			> $molecule_id-refmac-refines.log <<EOF  

LABIN FP=Fout0 PHIB=Pout0

# set resolution limit. The FSC=0.143 cutoff is recommended. 
RESO 3.5

# set B-factors at 40 prior to refinement:
##BFACtor SET 40

# specify map sharpening to be used during refinement:
##REFI sharpen -50

# specify number of refinement cycles:
NCYCLE 30

# set refinement weight:
WEIGht MATRix 0.01

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
          echo Output $molecule_id-refined.pdb  
          #fi
          #cd - 
      #fi
#done


"""