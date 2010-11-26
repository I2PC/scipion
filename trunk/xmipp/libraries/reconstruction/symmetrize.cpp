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

#include "symmetrize.h"

#include <data/args.h>

/* Read parameters --------------------------------------------------------- */
void ProgSymmetrize::readParams()
{
    fn_in  = getParam("-i");
    fn_out = getParam("-o");
    fn_sym = getParam("--sym");
    do_not_generate_subgroup = checkParam("--no_group");
    wrap = !checkParam("--dont_wrap");
    useBsplines = checkParam("--splines");
}

/* Usage ------------------------------------------------------------------- */
void ProgSymmetrize::defineParams()
{
    addUsageLine("Symmetrize a 3D volume ");
    addUsageLine("Example of use: Sample at i3 symmetry and the volume is not wrapped");
    addUsageLine("   xmipp_symmetrize -i input.vol --sym i3 --dont_wrap");
    addParamsLine("    -i <file_in>        : input 3D Xmipp file");
    addParamsLine("   [-o <file_out=\"\">]      : if no name is given then the input file is rewritten");
    addParamsLine("    --sym <sym_file>     : symmetry file (see the manual)");
    addParamsLine("   [--no_group]          : do not generate symmetry subgroup");
    addParamsLine("   [--dont_wrap]         : by default, the volume is wrapped");
    addParamsLine("   [--splines]           : by default, trilinear interpolation");
}

/* Symmetrize ------------------------------------------------------- */
//#define DEBUG
void symmetrize(const SymList &SL, MultidimArray<double> &V_in, MultidimArray<double> &V_out,
                                int Splinedegree,bool wrap, bool show_progress, bool do_outside_avg)
{
    Matrix2D<double> L(4, 4), R(4, 4); // A matrix from the list
    MultidimArray<double> V_aux, V_aux2;
    Matrix1D<double> sh(3);
    double dum, avg = 0.;

    if (do_outside_avg)
    {
        MultidimArray<int> mask;
        int rad;
        mask.resize(V_in);
        mask.setXmippOrigin();
        rad = XMIPP_MIN(V_in.ydim, V_in.xdim);
        rad = XMIPP_MIN(rad, V_in.zdim);
        BinaryCircularMask(mask, rad / 2, OUTSIDE_MASK);
        computeStats_within_binary_mask(mask, V_in, dum, dum, avg, dum);
    }
    V_out = V_in;

    if (show_progress)
    {
        std::cerr << "Symmetrizing ...\n";
        init_progress_bar(SL.SymsNo());
    }

    for (int i = 0; i < SL.SymsNo(); i++)
    {
        SL.get_matrices(i, L, R);

        SL.get_shift(i, sh);
        R(3, 0) = sh(0) * V_aux.xdim;//colNumber();
        R(3, 1) = sh(1) * V_aux.ydim;//rowNumber();
        R(3, 2) = sh(2) * V_aux.zdim;//sliceNumber();

        applyGeometry(Splinedegree,V_aux, V_in, R.transpose(), IS_NOT_INV, wrap, avg);
        //#define DEBUG
#ifdef DEBUG

        V_aux.write((std::string)"PPPsym_" + integerToString(i) + ".vol");
#endif

        /* *** CO: I am not very sure about the reason for this, but it
           seems to work */
        // applyGeometry(V_aux2(),L,V_aux(),IS_NOT_INV,prm.wrap);
        arrayByArray(V_out, V_aux, V_out, '+');
        if (show_progress)
            progress_bar(i);
    }
    if (show_progress)
        progress_bar(SL.SymsNo());
    arrayByScalar(V_out, SL.SymsNo() + 1.0f, V_out, '/');
}
#undef DEBUG

/* Main program ------------------------------------------------------------ */
void ProgSymmetrize::run()
{
    SymList           SL;
    Image<double>     V_in;
    Image<double>     V_out;

    double accuracy = (do_not_generate_subgroup) ? -1 : 1e-6;
    SL.read_sym_file(fn_sym, accuracy);
    std:: cout << "Number of symmetries: " << SL.SymsNo() << std::endl;
    V_in.read(fn_in);

    //std::cerr << prm;
    if (!useBsplines)
        symmetrize(SL, V_in(), V_out(), LINEAR, wrap, true,false);
    else
        symmetrize(SL, V_in(), V_out(), 3, wrap, true,true);
    if (fn_out == "")
        fn_out = V_in.name();
    V_out.write(fn_out);
}


