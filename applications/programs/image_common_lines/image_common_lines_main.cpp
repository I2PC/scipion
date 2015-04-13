/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <reconstruction/common_lines.h>

RUN_XMIPP_PROGRAM(ProgCommonLine)
#ifdef NEVERDEFINED
int main(int argc, char *argv[])
{
	try
    {
        int k = atoi(argv[1]);
        Matrix2D<double> quaternions;
        std::cerr << "DEBUG_JM: calling randomQuaternions" << std::endl;
        randomQuaternions(k, quaternions);

        //std::cerr << "DEBUG_JM: quaternions: " << quaternions << std::endl;
        Matrix2D<double> clMatrix, clCorr;
        size_t L = 10e+15;
        L = 3600;

        commonlineMatrixCheat(quaternions, L, clMatrix, clCorr);
        clMatrix.resize(100,100);
        clMatrix.read("Yoel/commonlines.txt");

        std::cerr << "DEBUG_JM: clMatrix: " << clMatrix << std::endl;

        int k1(0), k2(1), k3(2);
        DMatrix sMatrix, R;

        //  tripletRotationMatrix(clMatrix, L, 0, 1, 2, R);
        //  std::cerr << "DEBUG_JM: R: " << R << std::endl;
        std::cerr << "DEBUG_JM: calling syncMatrix" << std::endl;
        computeSyncMatrix(clMatrix, L, sMatrix, &quaternions);

        //std::cerr << "DEBUG_JM: sMatrix: " << sMatrix << std::endl;

        rotationsFromSyncMatrix(sMatrix, &quaternions);
    }
    catch (XmippError &xe)
    {
        std::cout << xe;
        exit(1);
    }
}
#endif
