/***************************************************************************
 *
 * Authors:    Sjors Scheres                     scheres@cnb.csic.es
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "angular_neighbourhood.h"

#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_image.h>

// Read arguments ==========================================================
void ProgAngularNeighbourhood::readParams()
{
    fn_sel = getParam("-i1");
    fn_ref = getParam("-i2");
    fn_out = getParam("-o");
    maxdist = getDoubleParam("--dist");
    fn_sym = getParam("--sym");
    if (fn_sym != "")
        SL.readSymmetryFile(fn_sym);
    check_mirrors=checkParam("--check_mirrors");
    DF2.read(fn_ref);
    DF1.read(fn_sel);
}

// Show ====================================================================
void ProgAngularNeighbourhood::show()
{
    std::cerr
    << "Selfile                      : " << fn_sel  << std::endl
    << "Neighbourhoods docfile       : " << fn_ref  << std::endl
    << "Output root                  : " << fn_out  << std::endl
    << "Max. neighbour distance      : " << maxdist << std::endl
    << "Symmetry file                : " << fn_sym  << std::endl
    ;
}

// usage ===================================================================
void ProgAngularNeighbourhood::defineParams()
{
    addUsageLine("Computes the images in Metadata1 that are within the neighbourhood of");
    addUsageLine("images in Metadata2. ");
    addParamsLine("   -i1 <metadata>          : Metadata1 with images and their angles");
    addParamsLine("   -i2 <metadata>          : Metadata2 with images and their angles");
    addParamsLine("   -o <metadata>           : Output metadata ");
    addParamsLine("                           :+For each image in Metadata2, there is a block ");
    addParamsLine("                           :+in the output metadata corresponding to all ");
    addParamsLine("                           :+images in Metadata1 within its neighbourhood.");
    addParamsLine("  [--dist  <d=10>        ] : Maximum neighbourhood distance ");
    addParamsLine("  [--sym <symmetry=\"\"> ] : Symmetry description if any");
    addParamsLine("                             :+The definition of the symmetry is described at [[transform_symmetrize_v3][transform_symmetrize]]");
    addParamsLine("  [--check_mirrors]        : Check mirrors");
    addSeeAlsoLine("angular_distance");
}

// Compute Projection Neighbourhood -----------------------------------------
void ProgAngularNeighbourhood::run()
{
    show();
    fn_out.deleteFile();
    MetaData SF_out;
    FileName fn2, fn1;

    std::cerr << "Calculating ...\n";
    FOR_ALL_OBJECTS_IN_METADATA(DF2)
    {
        // Read reference projection direction
        double rot2, tilt2;
        DF2.getValue(MDL_IMAGE,fn2,__iter.objId);
        DF2.getValue(MDL_ANGLE_ROT,rot2,__iter.objId);
        DF2.getValue(MDL_ANGLE_TILT,tilt2,__iter.objId);
        FOR_ALL_OBJECTS_IN_METADATA(DF1)
        {
            // Read assigned angles from document file
            double rot1, tilt1;
            DF1.getValue(MDL_ANGLE_ROT,rot1,__iter.objId);
            DF1.getValue(MDL_ANGLE_TILT,tilt1,__iter.objId);
            double psi1=0.;
            double dist = SL.computeDistance(rot2, tilt2, 0., rot1, tilt1, psi1,
                                             true, check_mirrors, false);
            if (dist <= maxdist)
            {
                DF1.getValue(MDL_IMAGE,fn1,__iter.objId);
                SF_out.setValue(MDL_IMAGE, fn1, SF_out.addObject());
            }
        }

        // finished reading all particles for this neighbourhood
        SF_out.write((String)"neighbourhoodOf_"+fn2+"@"+fn_out,MD_APPEND);
        SF_out.clear();
    }
}
