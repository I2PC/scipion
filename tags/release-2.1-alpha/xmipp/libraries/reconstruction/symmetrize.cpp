/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "symmetrize.h"

#include <data/args.h>

/* Read parameters --------------------------------------------------------- */
void Symmetrize_Parameters::read(int argc, char **argv)
{
    fn_in  = getParameter(argc, argv, "-i");
    fn_out = getParameter(argc, argv, "-o", "");
    fn_sym = getParameter(argc, argv, "-sym");
    do_not_generate_subgroup = checkParameter(argc, argv, "-no_group");
    wrap = !checkParameter(argc, argv, "-dont_wrap");
    useBsplines = checkParameter(argc, argv, "-splines");
}

/* Usage ------------------------------------------------------------------- */
void Symmetrize_Parameters::usage()
{
    std::cout
        << "Usage: symmetrize [Purpose and Parameters]\n"
        << "Purpose: Symmetrize a 3D volume\n"
        << "Parameter Values: (notice space before value)\n"
        << "    -i <file_in>        : input 3D Xmipp file\n"
        << "   [-o <file_out>]      : if no name is given then the input file is\n"
        << "                          rewritten\n"
        << "    -sym <sym_file>     : symmetry file (see the manual)\n"
        << "   [-no_group]          : do not generate symmetry subgroup\n"
        << "   [-dont_wrap]         : by default, the volume is wrapped\n"
        << "   [-splines]           : by default, trilinear interpolation\n";
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &out, const Symmetrize_Parameters &prm)
{
    out << "File in:       " << prm.fn_in  << std::endl
        << "File out:      " << prm.fn_out << std::endl
        << "Symmetry file: " << prm.fn_sym << std::endl
        << "Generate group:" << !prm.do_not_generate_subgroup << std::endl
        << "Wrapping:      ";
    print(out, prm.wrap); out << std::endl;
    out << "Splines:       ";
    print(out, prm.useBsplines); out << std::endl;
    return out;
}

/* Really symmetrize ------------------------------------------------------- */
//#define DEBUG
void symmetrize(const SymList &SL, VolumeXmipp &V_in, VolumeXmipp &V_out,
                bool wrap, bool show_progress)
{
    Matrix2D<double> L(4, 4), R(4, 4); // A matrix from the list
    VolumeXmipp V_aux, V_aux2;
    Matrix1D<double> sh(3);
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
        R(3, 0) = sh(0) * V_aux().colNumber();
        R(3, 1) = sh(1) * V_aux().rowNumber();
        R(3, 2) = sh(2) * V_aux().sliceNumber();

        /* *** CO: I don't know why the compiler doesn't allow me
           to reuse V_in !!!, this is very memory wasting */
        applyGeometry(V_aux(), R.transpose(), V_in(), IS_NOT_INV, wrap);
#ifdef DEBUG
        V_aux.write((std::string)"PPPsym_" + integerToString(i) + ".vol");
#endif

        /* *** CO: I am not very sure about the reason for this, but it
           seems to work */
        // applyGeometry(V_aux2(),L,V_aux(),IS_NOT_INV,prm.wrap);
        arrayByArray(V_out(), V_aux(), V_out(), '+');
        if (show_progress) progress_bar(i);
    }
    if (show_progress) progress_bar(SL.SymsNo());
    arrayByScalar(V_out(), SL.SymsNo() + 1.0f, V_out(), '/');
}
#undef DEBUG

/* Really symmetrize using Bsplines ------------------------------------------------ */
//#define DEBUG
void symmetrize_Bspline(const SymList &SL, VolumeXmipp &V_in, VolumeXmipp &V_out,
                        int Splinedegree, bool wrap, bool do_outside_avg)
{

    Matrix2D<double> L(4, 4), R(4, 4); // A matrix from the list
    VolumeXmipp V_aux, V_aux2;
    Matrix1D<double> sh(3);
    double dum, avg = 0.;

    if (do_outside_avg)
    {
        Matrix3D<int> mask;
        int rad;
        mask.resize(V_in());
        mask.setXmippOrigin();
        rad = XMIPP_MIN(V_in().rowNumber(), V_in().colNumber());
        rad = XMIPP_MIN(rad, V_in().sliceNumber());
        BinarySphericalMask(mask, rad / 2, OUTSIDE_MASK);
        computeStats_within_binary_mask(mask, V_in(), dum, dum, avg, dum);
    }

    V_out = V_in;

    for (int i = 0; i < SL.SymsNo(); i++)
    {
        SL.get_matrices(i, L, R);

        SL.get_shift(i, sh);
        R(3, 0) = sh(0) * V_aux().colNumber();
        R(3, 1) = sh(1) * V_aux().rowNumber();
        R(3, 2) = sh(2) * V_aux().sliceNumber();

        applyGeometryBSpline(V_aux(), R.transpose(), V_in(), Splinedegree, IS_NOT_INV, wrap, avg);
        arrayByArray(V_out(), V_aux(), V_out(), '+');

#ifdef DEBUG
        V_aux.write((std::string)"PPPsym_" + integerToString(i) + ".vol");
#endif

    }
    arrayByScalar(V_out(), SL.SymsNo() + 1.0f, V_out(), '/');

}
#undef DEBUG

/* Main program ------------------------------------------------------------ */
void ROUT_symmetrize(const Symmetrize_Parameters &prm)
{
    SymList         SL;
    VolumeXmipp     V_in;
    VolumeXmipp     V_out;
    FileName        fn_out;

    double accuracy = (prm.do_not_generate_subgroup) ? -1 : 1e-6;
    SL.read_sym_file(prm.fn_sym, accuracy);
    std:: cout << "Number of symmetries: " << SL.SymsNo() << std::endl;
    V_in.read(prm.fn_in);

    std::cerr << prm;
    if (!prm.useBsplines)
        symmetrize(SL, V_in, V_out, prm.wrap, true);
    else
        symmetrize_Bspline(SL, V_in, V_out, 3, prm.wrap, true);
    if (prm.fn_out == "") fn_out = V_in.name();
    else fn_out = prm.fn_out;
    V_out.write(fn_out);
}
