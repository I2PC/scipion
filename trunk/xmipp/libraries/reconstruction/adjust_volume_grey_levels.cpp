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

#include "adjust_volume_grey_levels.h"
#include <data/projection.h>

#include <data/args.h>

/* Read parameters --------------------------------------------------------- */
void Prog_Adjust_Volume_Parameters::read(int argc, char **argv)
{
    fn_vol = getParameter(argc, argv, "-i");
    fn_sel = getParameter(argc, argv, "-sel");
    fn_out = getParameter(argc, argv, "-o", "");
    optimize = checkParameter(argc, argv, "-optimize");
    probb_eval = textToFloat(getParameter(argc, argv, "-probb_eval", "0.2"));
    produce_side_info();
}

/* Usage ------------------------------------------------------------------- */
void Prog_Adjust_Volume_Parameters::usage()
{
    std::cerr << "Usage: adjust_volume\n";
    std::cerr << "   -i <Volume>         : Input volume\n"
    << "   -sel MetaDataFile      : Set of projections\n"
    << "  [-o <Output Volume>] : By default, the input one\n"
    << "  [-optimize]          : Optimize\n"
    << "  [-probb_eval <p=0.2>]: Probability of being evaluated\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Adjust_Volume_Parameters::show()
{
    std::cout << "Input Volume:  " << fn_vol   << std::endl
    << "Input MetaDAtaFile: " << fn_sel   << std::endl
    << "Output Volume: " << fn_out   << std::endl
    << "Optimize:      " << optimize << std::endl
    ;
}


/* Produce side information ------------------------------------------------ */
void Prog_Adjust_Volume_Parameters::produce_side_info()
{
    // Read input volume
    VolumeXmipp IV;
    IV.read(fn_vol);
    V = IV();
    V.setXmippOrigin();

    // Read input metadataFile
    SF.read(fn_sel,NULL);
}

/* Goal function -----------------------------------------------------------  */
Prog_Adjust_Volume_Parameters *global_adjust_volume_prm;

double projection_mismatching(double *p, void *prm)
{
    return global_adjust_volume_prm->mismatching(p[1], p[2]);
}

//#define DEBUG
double Prog_Adjust_Volume_Parameters::mismatching(double a, double b)
{
    if (a <= 0) return 1e38;

    // Transform the volume
    Matrix3D<double> aux = V;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(aux) aux(k, i, j) = a * aux(k, i, j) + b;

    // Compute the mismatching
    double retval = 0;
    //SF.go_first_ACTIVE();
    int N = 0;
    //while (!SF.eof())
    long int ret=SF.firstObject();
    if(ret==MetaData::NO_OBJECTS_STORED)
    {
            std::cerr << "Empty inputFile File\n";
            exit(1);
    }
    do
    {
        // Read next image
        FileName fn = SF.image();
        if (fn=="") break;

        // Skip randomly some images
        double x = rnd_unif(0, 1);
        if (x > probb_eval) continue;
        N++;

        ImageXmipp I;
        I.read(fn);
        I().setXmippOrigin();

        // Project the auxiliary volume in the same direction
        Projection P;
        project_Volume(aux, P, YSIZE(I()), XSIZE(I()),
                       I.rot(), I.tilt(), I.psi());

        // Compute the difference
        Matrix2D<double> diff;
        diff = I() - P();
        retval += diff.sum2();

#ifdef DEBUG
        I.write("PPPexp.xmp");
        I() = P();
        I.write("PPPtheo.xmp");
        std::cout << "Difference=" << diff.sum2() << std::endl;
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
#endif
    }
    while (SF.nextObject()!= MetaData::NO_MORE_OBJECTS);

    return retval / N;
}
#undef DEBUG

/* Apply ------------------------------------------------------------------- */
void Prog_Adjust_Volume_Parameters::apply(Matrix3D<double> &out)
{
    // Compute the average power and average value of all the projections
    double sum = 0, sum2 = 0, N = 0;
    std::cerr << "Computing first estimate of the linear transformation ...\n";
    int imgno = SF.size();
    init_progress_bar(imgno);
    int i = 0;

    int projXdim, projYdim;
    long int ret=SF.firstObject();
    if(ret==MetaData::NO_OBJECTS_STORED)
    {
            std::cerr << "Empty inputFile File\n";
            exit(1);
    }
    do
    {
        // Read image
        FileName fn = SF.image();
        if (fn=="") break;
        ImageXmipp I;
        I.read(fn);
        projXdim = XSIZE(I());
        projYdim = YSIZE(I());

        // Compute the image statistics
        double avg, stddev, min, max;
        I().computeStats(avg, stddev, min, max);
        double Ni = projXdim * projYdim;
        sum += avg;
        sum2 += stddev * stddev;
        N += Ni;

        // End of loop
        i++;
        if (i % 10 == 0) progress_bar(i);
    }
    while (SF.nextObject()!= MetaData::NO_MORE_OBJECTS);
    progress_bar(imgno);
    std::cout << std::endl;

    // Statistics of the volume
    double avg0, stddev0, min0, max0;
    V.computeStats(avg0, stddev0, min0, max0);

    // First guess of the transformation parameters a*(x-vm)+b
    // r is the average length of a ray
    // The equations to solve are
    // v_i follows currently N(avg0,stddev0)
    // we want it to follow N(avgF,stddevF)
    // such that sum_{i=1}^r{v_i} follows N(avg_pict,stddev_pict)
    // On the other hand we know that
    //     sum_{i=1}^r{v_i} follows N(r avgF, sqrt(r)*stddevF)
    // Thus,
    //     avgF=avg_pict/r
    //     stddevF=stddev_pict/sqrt(r)
    // the equations of the transformation are
    //    y=ax+b
    //    avg_y     =a*avg_x+b=avgF
    //    stddev_y  =a*stddev_y=stddevF
    // Therefore
    //    a=stddevF/stddev0
    //    b=avgF-a*avg0
    double r = pow(MULTIDIM_SIZE(V), 1.0 / 3.0);
    double avg_pict = sum / imgno;
    double stddev_pict = sqrt(sum2 / imgno);
    double avgF = avg_pict / r;
    double stddevF = stddev_pict / sqrt(r);
    double a = stddevF / stddev0;
    double b = avgF - a * avg0;
    std::cout << "First Linear transformation: y=" << a << "*x+"  << b << std::endl;

    // Optimize
    if (optimize)
    {
        Matrix1D<double> p(2), steps(2);
        p(0) = a;
        p(1) = b;
        steps.initConstant(1);
        double ftol = 0.01, fret;
        int iter;
        global_adjust_volume_prm = this;
        powellOptimizer(p, 1, 2, &projection_mismatching, NULL,
            ftol, fret, iter, steps, true);
        a = p(0);
        b = p(1);
    }

    // Apply the transformation
    out = V;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(V) out(k, i, j) = a * V(k, i, j) + b;
}

/* Run --------------------------------------------------------------------- */
void Prog_Adjust_Volume_Parameters::run()
{
    VolumeXmipp out;
    apply(out());
    if (fn_out == "") out.write(fn_vol);
    else            out.write(fn_out);
}
