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
#include <data/numerical_tools.h>
#include <data/projection.h>
#include <data/xmipp_image.h>
#include <data/args.h>

void ProgAdjustVolume::defineParams()
{
    // Usage
    addUsageLine("Search for a linear transformation of the gray values of a volume so that the");
    addUsageLine("error between the theoretically projected images and the experimental images");
    addUsageLine("is minimized. This program must be used before computing the Volumetric SSNR");
    addUsageLine("if the reconstruction algorithm scales the output volume differently.");
    // See also
    addSeeAlsoLine("resolution_ssnr");

    // Examples
    addExampleLine("Adjust a volume to a set of images:", false);
    addExampleLine("xmipp_adjust_volume_grey_levels -i input_volume.vol -m experimental.sel -o output_volume.vol ");

    // Parameters
    addParamsLine(" -i <volume_file>      : Volume to adjust its range.");
    addParamsLine(" alias --input;");
    addParamsLine(" -m  <metadata_file>   : Set of projections of the volume.");
    addParamsLine(" alias --metadata;");
    addParamsLine(" [-o <volume_file=\"\">]    : Output adjusted volume. By default, the input one.");
    addParamsLine(" alias --output;");
    addParamsLine(" [--optimize]          : Optimize the linear transformation. By default, keep the initially computed.");
    addParamsLine(" [--probb_eval <p=0.2>] : The goal function is evaluated each time from a random set of projections.");
    addParamsLine("                       : Each image has a probability of p of being evaluated.");
}

void ProgAdjustVolume::readParams()
{
    fn_vol = getParam("-i");
    fn_sel = getParam("-m");
    fn_out = (checkParam("-o"))? getParam("-o") : fn_vol;
    optimize = checkParam("--optimize");
    probb_eval = getDoubleParam("--probb_eval");
    verbose = (checkParam("-v"))? getIntParam("-v"): false;
    tempFile = (STR_EQUAL(fn_vol.c_str(), fn_out.c_str()));
}

void ProgAdjustVolume::run()
{
    if (verbose)
        show();

    // Read input volume
    ImIn.read(fn_vol);
    V.alias(ImIn());
    V.setXmippOrigin();

    int Xdim, Ydim, Zdim;
    Xdim = XSIZE(V);
    Ydim = YSIZE(V);
    Zdim = ZSIZE(V);

    Image<float> ImOut;
    ImOut.mapFile2Write(Xdim, Ydim, Zdim, fn_out, tempFile);
    ImOut().setXmippOrigin();

    // Read input metadataFile
    SF.read(fn_sel,NULL);

    apply(ImOut());

    ImIn.clear();
    ImOut.write();
}

/* Goal function -----------------------------------------------------------  */
ProgAdjustVolume *globalAdjustVolumeProg;
double projectionMismatching(double *p, void *prm)
{
    return globalAdjustVolumeProg->mismatching(p[1], p[2]);
}

//#define DEBUG
double ProgAdjustVolume::mismatching(double a, double b)
{
    if (a <= 0)
        return 1e38;

    // Transform the volume
    MultidimArray<double> aux = V;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(aux) aux(k, i, j) = a * aux(k, i, j) + b;

    // Compute the mismatching
    double retval = 0;
    //SF.go_first_ACTIVE();
    int N = 0;
    //while (!SF.eof())
    if(SF.isEmpty())
    {
        std::cerr << "Empty inputFile File\n";
        exit(1);
    }
    Image<double> I;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        // Skip randomly some images
        double x = rnd_unif(0, 1);
        if (x > probb_eval)
            continue;
        N++;

        I.readApplyGeo(SF,__iter.objId);
        I().setXmippOrigin();

        // Project the auxiliary volume in the same direction
        Projection P;
        projectVolume(aux, P, YSIZE(I()), XSIZE(I()),
                      I.rot(), I.tilt(), I.psi());

        // Compute the difference
        MultidimArray<double> diff;
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
    //TODO Review this
    //    while (SF.nextObject()!= NO_MORE_OBJECTS)
    //        ;

    return retval / N;
}
#undef DEBUG

/* Apply ------------------------------------------------------------------- */
void ProgAdjustVolume::apply(MultidimArray<float> &out)
{
    // Compute the average power and average value of all the projections
    double sum = 0, sum2 = 0, N = 0;
    std::cerr << "Computing first estimate of the linear transformation ...\n";
    int imgno = SF.size();
    init_progress_bar(imgno);
    int i = 0;

    int projXdim, projYdim;
    if(SF.isEmpty())
    {
        std::cerr << "Empty inputFile File\n";
        exit(1);
    }
    Image<double> I;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        // Read image
        I.readApplyGeo(SF,__iter.objId);
        projXdim = XSIZE(I());
        projYdim = YSIZE(I());

        // Compute the image statistics
        double avg=0., stddev=0., min, max;
        I().computeStats(avg, stddev, min, max);
        double Ni = projXdim * projYdim;
        sum += avg;
        sum2 += stddev * stddev;
        N += Ni;

        // End of loop
        i++;
        if (i % 10 == 0)
            progress_bar(i);
    }
    progress_bar(imgno);
    std::cout << std::endl;

    // Statistics of the volume
    double avg0=0., stddev0=0., min0, max0;
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
        globalAdjustVolumeProg = this;
        powellOptimizer(p, 1, 2, &projectionMismatching, NULL,
                        ftol, fret, iter, steps, true);
        a = p(0);
        b = p(1);
    }

    // Apply the transformation
//    out = V;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(V) out(k, i, j) = a * V(k, i, j) + b;
}

/* Show -------------------------------------------------------------------- */
void ProgAdjustVolume::show()
{
    std::cout << "Input Volume:  " << fn_vol   << std::endl
    << "Input MetaDAtaFile: " << fn_sel   << std::endl
    << "Output Volume: " << fn_out   << std::endl
    << "Optimize:      " << optimize << std::endl;
}

