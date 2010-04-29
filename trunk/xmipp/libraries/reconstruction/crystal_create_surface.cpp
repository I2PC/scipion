/***************************************************************************
 *
 * Authors:
 *
 * Javier Rodríguez Falces (jrodriguez@cnb.csic.es)
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

#include "crystal_create_surface.h"
#include "fourier_filter.h"
#include <data/args.h>
#include <data/fft.h>
#include <data/funcs.h>
#include <data/docfile.h>
#include <data/integration.h>
#include <data/error.h>
#include <data/blobs.h>

#include <fstream>

/* Default constructor ------------------------------------------------------- */
Prog_create_surface::Prog_create_surface()
{
    fn_in = "";
    init_random_generator(); // random initializate
    a.resize(2);
    b.resize(2);
}


/* Usage ------------------------------------------------------------------- */
void Prog_create_surface::usage()
{
    std::cout << "   -i <pdb file>                    : File to read data from\n"
              << "   -o <pdb file>                    : File to write result data in\n"
              << "   -f <string>                      : Name of the surface to map\n";
}

/* Read parameters --------------------------------------------------------- */
void Prog_create_surface::read(int argc, char **argv)
{
    fn_in = getParameter(argc, argv, "-i");
    fn_out = getParameter(argc, argv, "-o", "");
    option = getParameter(argc, argv, "-f", ""); // 1:parabole, 2:cosine
    if (fn_out == "") fn_out = fn_in.without_extension() + ".out";
}

void Prog_create_surface::read_input_file()
{
    FILE    *fh_param;
    char    line[201];
    int     lineNo = 0;
    char    *auxstr;
    if ((fh_param = fopen(fn_in.c_str(), "r")) == NULL)
        REPORT_ERROR(3005,
                     (std::string)"Prog_Project_Parameters::read: There is a problem "
                     "opening the file " + fn_in);

    while (fgets(line, 200, fh_param) != NULL)
    {
        if (line[0] == 0)    continue;
        if (line[0] == '#')  continue;
        if (line[0] == '\n') continue;
        switch (lineNo)
        {
        case 0:
            a.X() = textToFloat(firstToken(line), 3007,
                                "Prog_Project_Parameters::read: Error in lattice vector a");
            a.Y() = textToFloat(nextToken(), 3007,
                                "Prog_Project_Parameters::read: Error in lattice vector a");
            lineNo = 1;
            break;
        case 1:
            b.X() = textToFloat(firstToken(line), 3007,
                                "Prog_Project_Parameters::read: Error in lattice vector b");
            b.Y() = textToFloat(nextToken(), 3007,
                                "Prog_Project_Parameters::read: Error in lattice vector b");
            lineNo = 2;
            break;
        case 2:
            hmin = textToInteger(firstToken(line), 3007,
                                 "Prog_Project_Parameters::read: Error in the index hmax of the lattice");
            hmax = textToInteger(nextToken(), 3007,
                                 "Prog_Project_Parameters::read: Error in the index hmin of the lattice");
            lineNo = 3;
            break;
        case 3:
            kmin = textToInteger(firstToken(line), 3007,
                                 "Prog_Project_Parameters::read: Error in the index kmax of the lattice");
            kmax = textToInteger(nextToken(), 3007,
                                 "Prog_Project_Parameters::read: Error in the index kmin of the lattice");
            lineNo = 4;
            break;
        } /* switch end */
    } /* while end */
    if (lineNo != 4)
        REPORT_ERROR(3007, (std::string)"Prog_Project_Parameters::read: I "
                     "couldn't read all parameters from file " + fn_in);

    fclose(fh_param);
}



// Function to generate the cosine surface
class Func1: public doubleFunction
{
    //actual function to be integrated
public:           //This should be in testinteg
    double x;
    double cte1, cte2, contribution_x;
    virtual double operator()()    //overloads pure virtual
    {
        return sqrt(contribution_x + cte1*sin(cte2*x)*sin(cte2*x));
    }
};

// Function to generate the parabolic surface
class Func2: public doubleFunction
{
    //actual function to be integrated
public:          //This should be in testinteg
    double x;
    double Aperture;
    virtual double operator()()    //overloads pure virtual
    {
        return sqrt(1 + 4.*Aperture*Aperture*x*x);
    }
};




// The maping_function maps each (h,k) point in a plane surface into a parabole and/or a cosine surface
void Prog_create_surface::maping_function(Matrix1D<double > a, Matrix1D<double > b, int h, int k, surface_coordinates &result)
{

    int cont;
    double uplimit, lowlimit, middle, vx;
    double integralt, inte_low = 0, inte_high, cteA, cteB, coordx;
    Matrix1D<double> vec(2), normalized_vec(2), normal_vec(3);
    Func1 cosine;   //cosine surface
    Func2 parabole; //parabole surface

//*
//For each (h,k) point we calculate its "lattice" vector
//*
    vec = h * a + k * b;

//*
//setting of constants and variables depending on the surface implemented
//*
    if (option.compare("cosine") == 0)
    {
        coordx = fabs(vec.X());
        cteA = 100;
        cteB = 0.01;
        vx = vec.X();
        cosine.cte1 = fabs(cteA * cteA * cteB * cteB * vx * vx);
        cosine.cte2 = fabs(cteB * vx);
        cosine.contribution_x = fabs(vx * vx);
    }
    else if (option.compare("parabole") == 0)
    {
        coordx = vec.module();
        normalized_vec = vec;
        normalized_vec.selfNormalize();
        parabole.Aperture = 0.001;
    }
    else
    {
        std::cerr << "function " << option << " not implemented" << std::endl;
        exit(1);
    }
//*
// initialization values
//*
    uplimit = fabs(10 * coordx);
    lowlimit = 0;
    middle = (uplimit + lowlimit) / 2;
    inte_high = middle;
    cont = 0;
    while (1)
    {
        //*
        // Calculate the integral over a parabolic or cosine surface
        //*
        if (option.compare("cosine") == 0)
        {
            Trapeze Trap(cosine, cosine.x, inte_low, inte_high);
            integralt = Trap();
        }
        else if (option.compare("parabole") == 0)
        {
            Trapeze Trap(parabole, parabole.x, inte_low, inte_high);
            integralt = Trap();
        }

        //*
        // Dicotomic search
        //*
        if (integralt < coordx)
        {
            lowlimit = lowlimit + (uplimit - lowlimit) / 2;
            uplimit = uplimit;
            middle  = lowlimit + (uplimit - lowlimit) / 2;
        }
        else
        {
            lowlimit = lowlimit;
            uplimit = lowlimit + (uplimit - lowlimit) / 2;
            middle  = lowlimit + (uplimit - lowlimit) / 2;
        }
        inte_high = middle;
        cont++;
        if ((fabs(coordx - integralt) < 0.05) || cont > 20)
            break;
    }

    if (option.compare("cosine") == 0)
    {
        // coordinates in the plane surface
        result.cox_ideal = vec.X();
        result.coy_ideal = vec.Y();
        result.coz_ideal = 0;
        // coordinates in the cosine surface
        result.cox_real = vx * inte_high;
        result.coy_real = vec.Y();
        result.coz_real = cteA * cos(cteB * result.cox_real);
        // normal coordinates
        normal_vec.X() = + cteA * cteB * sin(cteB * result.cox_real) ;
        normal_vec.Y() = 0 ;
        normal_vec.Z() = + 1 ;
        normal_vec.selfNormalize();
        result.Nx = normal_vec.X();
        result.Ny = normal_vec.Y();
        result.Nz = normal_vec.Z();
    }
    else if (option.compare("parabole") == 0)
    {
        // coordinates in the plane surface
        result.cox_ideal = vec.X();
        result.coy_ideal = vec.Y();
        result.coz_ideal = 0;
        // coordinates in the parabole surface
        result.cox_real = normalized_vec.X() * inte_high;
        result.coy_real = normalized_vec.Y() * inte_high;
        result.coz_real = parabole.Aperture * (result.cox_real * result.cox_real + result.coy_real * result.coy_real);
        // normal coordinates
        normal_vec.X() = -2. * parabole.Aperture * result.cox_real ;
        normal_vec.Y() = -2. * parabole.Aperture * result.coy_real ;
        normal_vec.Z() = + 1 ;
        normal_vec.selfNormalize();
        result.Nx = normal_vec.X();
        result.Ny = normal_vec.Y();
        result.Nz = normal_vec.Z();
    }
}


/* Run --------------------------------------------------------------------- */
void Prog_create_surface::run()
{

    double x_des = 0, y_des = 0, z_des = 0.0; //Noise parameters (u and std)
    double x_nor_des = 0, y_nor_des = 0, z_nor_des = 0.0; //Noise parameters (u and std)
    Matrix1D<double> data_line(10);

    DocFile  DF_report_standard;
    surface_coordinates result;
    DF_report_standard.append_comment("Headerinfo columns: rot tilt psi x y corr");

    for (int h = hmin; h <= hmax; h++)
    {
        for (int k = kmin; k <= kmax; k++)
        {

            maping_function(a, b, h, k, result);

            data_line(0) = h; // h
            data_line(1) = k; // k
            data_line(2) = result.cox_ideal; // x
            data_line(3) = result.coy_ideal; // y
            data_line(4) = result.cox_real - result.cox_ideal  + rnd_gaus(0, x_des); // Delta_x
            data_line(5) = result.coy_real - result.coy_ideal  + rnd_gaus(0, y_des); // Delta_y
            data_line(6) = result.coz_real + rnd_gaus(0, z_des); // Delta_z
            data_line(7) = result.Nx + rnd_gaus(0, x_nor_des); // N_x
            data_line(8) = result.Ny + rnd_gaus(0, y_nor_des); // N_y
            data_line(9) = result.Nz + rnd_gaus(0, z_nor_des); // N_z
            DF_report_standard.append_data_line(data_line);
        }
    }
    DF_report_standard.write(fn_out);
}
