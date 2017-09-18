/***************************************************************************
 *
 * Authors:    Tomas Majtner            tmajtner@cnb.csic.es (2017)
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

#include "evaluate_coordinates.h"

#include <data/xmipp_funcs.h>
#include <data/mask.h>
#include <data/filters.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>


// Read arguments ==========================================================
void ProgEvaluateCoordinates::readParams()
{
    fnGt = getParam("-g");
    fnEval = getParam("-e");
    errMargin = getIntParam("-t");
}

// Show ====================================================================
void ProgEvaluateCoordinates::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Ground truth coordinates:          " << fnGt         << std::endl
    << "Coordinates to evaluate:           " << fnEval       << std::endl
    << "Tolerance of center misplacement:  " << errMargin    << std::endl
    ;
}

// Usage ===================================================================
void ProgEvaluateCoordinates::defineParams()
{
    addUsageLine("Clusters a set of images");
    addParamsLine("  -g <selfile>     : Selfile containing ground truth coordinates");
    addParamsLine("  -e <selfile>     : Selfile containing coordinates to evaluate");
    addParamsLine("  [-t <int=10>]    : Tolerance of center misplacement");
}

void ProgEvaluateCoordinates::run()
{
    MDRow row;
    int evalXCoor, evalYCoor, gtXCoor, gtYCoor;
    int truePosivites = 0, totalEval = 0, totalGT = 0;

    for (int m = 1; m < 51; m++)
    {
        FileName micGT=formatString("mic_KLH_Dataset_I_Blind_Test_%04d@%s", m, fnGt.c_str());
        FileName micEval=formatString("mic_KLH_Dataset_I_Blind_Test_%04d@%s", m, fnEval.c_str());
        GT.read(micGT);
        Eval.read(micEval);

        int iterE = 1;
        FOR_ALL_OBJECTS_IN_METADATA(Eval)
        {
            Eval.getRow(row, iterE);
            row.getValue(MDL_XCOOR, evalXCoor);
            row.getValue(MDL_YCOOR, evalYCoor);

            int iterG = 1;
            FOR_ALL_OBJECTS_IN_METADATA(GT)
            {
                GT.getRow(row, iterG);
                row.getValue(MDL_XCOOR, gtXCoor);
                row.getValue(MDL_YCOOR, gtYCoor);

                if ((evalXCoor > (gtXCoor-errMargin-1)) && (evalXCoor < (gtXCoor+errMargin+1)) &&
                    (evalYCoor > (gtYCoor-errMargin-1)) && (evalYCoor < (gtYCoor+errMargin+1)))
                    truePosivites++;

                iterG++;
            }
            iterE++;
        }
        totalEval += Eval.size();
        totalGT += GT.size();
    }

    std::cout << std::endl;
    std::cout << "True positives (correctly picked particles): " << truePosivites << std::endl;
    std::cout << "False positives (incorrectly picked particles): " << totalEval - truePosivites << std::endl;
    std::cout << "False negatives (missed particles): " << totalGT - truePosivites << std::endl << std::endl;

    std::cout << "How many relevant particles are picked:" << std::endl;
    std::cout << "True positive rate: " << (double) truePosivites / totalGT << std::endl << std::endl;

    std::cout << "Ratio between wrongly picked particles and all picked particles:" << std::endl;
    std::cout << "False positive rate: " << (double) (totalEval - truePosivites) / totalEval << std::endl;
    std::cout << std::endl;
}