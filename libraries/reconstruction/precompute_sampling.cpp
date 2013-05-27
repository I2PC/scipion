/***************************************************************************
 *
 * Authors:
 *
 * Roberto Marabini
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

#include "precompute_sampling.h"

/* Empty constructor ------------------------------------------------------- */
Prog_Sampling_Parameters::Prog_Sampling_Parameters()
{
    /** sampling object 1 by default*/
    mysampling.setSampling(1);
}


/* Read parameters --------------------------------------------------------- */
void Prog_Sampling_Parameters::read(int argc, const char **argv)
{
    sampling_file_root = getParameter(argc, argv, "-o");
    fn_sym = getParameter(argc, argv, "-sym");
    //symmetry = getParameter(argc, argv, "-symmetry", "cn");
    //sym_order = textToInteger(getParameter(argc, argv, "-sym_order", "1"));
    sampling = textToFloat(getParameter(argc, argv, "-sampling_rate", "5"));
    neighborhood = textToFloat(getParameter(argc, argv, "-neighborhood", "1"));
    max_tilt_angle = textToFloat(getParameter(argc, argv, "-max_tilt_angle","91"));
    min_tilt_angle = textToFloat(getParameter(argc, argv, "-min_tilt_angle","-91"));
}

/* Usage ------------------------------------------------------------------- */
void Prog_Sampling_Parameters::usage()
{
    std::cerr << "precompute_sampling\n"
    << "   -o root_file_name           : Root for output files\n"
    << "  [-sym cn]   :One of the 17 possible symmetries in\n"
    << "                                single particle electronmicroscopy\n"
    << "                                i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv, dnh, t, td, th, o, oh, i, ih\n"
    << "                               : where n may change from 1 to 99\n"
    << "  [-sampling_rate 5]           : Distance in degrees between sampling points\n"
    << "  [-neighborhood 1]            : A sampling point is neighbor if closer than this value in degrees\n"
    << "  [-max_tilt_angle  91]        : maximum tilt angle in degrees\n"
    << "  [-min_tilt_angle -91]        : minimum tilt angle in degrees\n"
    << "\n"
    << "Example of use: Sample at 2degres and compute neighboor at "
    << " 5 degrees for c6 symmetry\n"
    << "   xmipp_precompute_sampling -o out -sym c6 "
    << " -sampling_rate 2 -neighborhood 5 -max_tilt_angle 70 "
    << " -min_tilt_angle 50 \n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Sampling_Parameters::show()
{
    std::cout
    << "Sampling rate:      " << sampling    << std::endl
    << "output files root:  " << sampling_file_root << std::endl
    << "symmetry group:            " << fn_sym << std::endl
    //<< "symmetry order:     " << sym_order << std::endl
    << "neighborhood:       " << neighborhood << std::endl
    << "max_tilt_angle:     " << max_tilt_angle << std::endl
    << "min_tilt_angle:     " << min_tilt_angle << std::endl
    ;
}



/* Run --------------------------------------------------------------------- */
void Prog_Sampling_Parameters::run()
{
    show();
    mysampling.setSampling(sampling);
    mysampling.setNeighborhoodRadius(neighborhood);
    mysampling.computeSamplingPoints(false,max_tilt_angle,min_tilt_angle);
    //mysampling.computeSamplingPoints(false);
    mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order);
    mysampling.removeRedundantPoints(symmetry, sym_order);
    mysampling.createAsymUnitFile(sampling_file_root);
    //mysampling.computeNeighbors();
    //#define DEBUG6
#ifdef DEBUG6
    for (int i = 0; i < mysampling.no_redundant_sampling_points_vector.size(); i++)
    {
        std::cout << mysampling.no_redundant_sampling_points_vector[i].transpose() << " 1.1 2.2222 " << std::endl;
        //std::cout << mysampling.no_redundant_sampling_points_angles[i].transpose() << " 1.21 1 " << std::endl;
    }
#endif
#undef DEBUG6
    //#define DEBUG6
#ifdef DEBUG6
    {   /* for sphere coverage only R is important.
                  for psi L and R are important
                */
        double rot,  tilt,  psi;
        double rotp, tiltp, psip;
        Matrix1D<double>  row(3);
        Matrix2D<double>  L(4, 4), R(4, 4);
        //std::cerr << "mysampling.SL.SymsNo():" << mysampling.SL.SymsNo() << std::endl;
        for (int i = 0; i < mysampling.no_redundant_sampling_points_vector.size(); i++)
        {
            if (i == 50)
            {
                rot = XX(mysampling.no_redundant_sampling_points_angles[i]);
                tilt = YY(mysampling.no_redundant_sampling_points_angles[i]);
                psi = ZZ(mysampling.no_redundant_sampling_points_angles[i]);
                std::cerr << 1 << " " << 3 << " " << rot << " " << tilt << " " << psi << std::endl;
                //std::cerr << mysampling.no_redundant_sampling_points_vector[i].transpose() << " 1 1" << std::endl;
            }
            for (int isym = 0; isym < mysampling.SL.SymsNo(); isym++)
            {
                //std::cerr << L << R << std::endl;
                mysampling.SL.get_matrices(isym, L, R);
                R.resize(3, 3); // as only the relative orientation
                row = R * mysampling.no_redundant_sampling_points_vector[i];
                std::cout << row.transpose() << " 1 " << isym + 2 << std::endl;
                if (i == 50)
                {
                    //std::cerr << row.transpose() << " 1 " << isym +2 << std::endl;
                    L.resize(3, 3); // Erase last row and column
                    //std::cerr << L << R << std::endl;
                    Euler_apply_transf(L, R, rot, tilt, psi, rotp, tiltp, psip);
                    Euler_direction(rotp, tiltp, psip, row);
                    std::cerr << isym + 2 << " " << 3 << " " << rotp << " " << tiltp << " " << psip << std::endl;
                }
            }
        }
    }
#endif
#undef DEBUG6
    //#define DEBUG6
#ifdef DEBUG6
    for (int i = 0;
         i < mysampling.sampling_points_vector.size();
         i++)
        std::cout  <<  mysampling.sampling_points_vector[i].transpose()  << " 1 1 " << std::endl;
#endif
#undef DEBUG6
    //#define DEBUG6
#ifdef DEBUG6
    for (int i = 0; i < mysampling.no_redundant_sampling_points_vector.size();i++)
    {
        std::cout  << mysampling.no_redundant_sampling_points_vector[i].transpose() << " 1.2 3 " << std::endl;
        for (int j = 0; j < mysampling.my_neighbors[i].size();j++)
            std::cout  <<
            mysampling.no_redundant_sampling_points_vector[mysampling.my_neighbors[i][j]].transpose()
            << " 1.1 2 " << std::endl;
    }
#endif
#undef DEBUG6
}

