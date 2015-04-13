/***************************************************************************
 *
 * Authors:    Roberto Marabini      roberto@cnb.csic.es (2015)
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

#include "angular_break_symmetry.h"

#include <data/args.h>
#include <data/histogram.h>

// Read arguments ==========================================================
void ProgAngularBreakSymmetry::readParams()
{
    fn_ang = getParam("-i");
    fn_out  = getParam("-o");
    fn_sym  = getParam("--sym");
}

// Show ====================================================================
void ProgAngularBreakSymmetry::show()
{
    std::cout
    << "Angular metadata File (IN): " << fn_ang       << std::endl
    << "Angular metadata (OUT)    : " << fn_out    << std::endl
    << "Symmetry file    : " << fn_sym        << std::endl
    ;
}

// usage ===================================================================
void ProgAngularBreakSymmetry::defineParams()
{
    addUsageLine("randomize angle input file ");
    addParamsLine("   -i <Metadata>        : Angular document file");
    addParamsLine("   -o <Metadata>        : Output file");
    addParamsLine("   --sym <symmetry>     : Symmetry file");
}

// Produce side information ================================================
void ProgAngularBreakSymmetry::produce_side_info()
{
    SL.readSymmetryFile(fn_sym);
    DF.read(fn_ang);
}

//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgAngularBreakSymmetry::run()
{
    produce_side_info();

    MDRow row;
    FOR_ALL_OBJECTS_IN_METADATA(DF)
    {
        // Read input data
        double rot,  tilt,  psi;

        DF.getValue(MDL_ANGLE_ROT,  rot,__iter.objId);
        DF.getValue(MDL_ANGLE_TILT, tilt,__iter.objId);
        DF.getValue(MDL_ANGLE_PSI,  psi,__iter.objId);


        // Apply random rotation within symmetry matrices
        SL.breakSymmetry(rot, tilt, psi,
                         rot, tilt, psi);

        // Bring both angles to a normalized set
        //rot2  = realWRAP(rot2, -180, 180);
        //tilt2 = realWRAP(tilt2, -180, 180);
        //psi2  = realWRAP(psi2, -180, 180);

        DF.setValue(MDL_ANGLE_ROT, rot,__iter.objId);
        DF.setValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        DF.setValue(MDL_ANGLE_PSI, psi,__iter.objId);
    }

    DF.write(fn_out);
}

