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

#include <data/args.h>
#include <data/image.h>
#include <data/filters.h>
#include <data/geometry.h>
#include <data/mask.h>

void Usage(const Mask_Params &m);

// Alignment parameters needed by fitness ----------------------------------
class AlignParams
{
public:
#define COVARIANCE       1
#define LEAST_SQUARES    2

    int alignment_method;

    Image<double> V1;
    Image<double> V2;
    Image<double> Vaux;
    const MultidimArray<int> *mask_ptr;
};

// Global parameters needed by fitness ------------------------------------
AlignParams prm;

// Apply transformation ---------------------------------------------------
void applyTransformation(const MultidimArray<double> &V2, 
                         MultidimArray<double> &Vaux,
                         double *p)
{
    Matrix1D<double> r(3);
    Matrix2D<double> A, Aaux;

    double greyScale = p[0];
    double greyShift = p[1];
    double rot       = p[2];
    double tilt      = p[3];
    double psi       = p[4];
    double scale     = p[5];
    ZZ(r)            = p[6];
    YY(r)            = p[7];
    XX(r)            = p[8];

    Euler_angles2matrix(rot, tilt, psi, A, true);
    translation3DMatrix(r,Aaux);
    A = A * Aaux;
    scale3DMatrix(vectorR3(scale, scale, scale),Aaux);
    A = A * Aaux;

    applyGeometry(LINEAR, Vaux, V2, A, IS_NOT_INV, WRAP);
    Vaux*=greyScale;
    Vaux+=greyShift;
}

// Fitness between two volumes --------------------------------------------
double fitness(double *p)
{
    applyTransformation(prm.V2(),prm.Vaux(),p);

    // Correlate
    double fit;
    switch (prm.alignment_method)
    {
    case (COVARIANCE):
                    fit = -correlation_index(prm.V1(), prm.Vaux(), prm.mask_ptr);
        break;
    case (LEAST_SQUARES):
                    fit = rms(prm.V1(), prm.Vaux(), prm.mask_ptr);
        break;
    }
    return fit;
}

double wrapperFitness(double *p, void *prm)
{
    return fitness(p+1);
}

int main(int argc, char **argv)
{
    FileName fn1, fn2;
    double   rot0, rotF, tilt0, tiltF, psi0, psiF;
    double   step_rot, step_tilt, step_psi;
    double   scale0, scaleF, step_scale;
    double   z0, zF, y0, yF, x0, xF, step_z, step_y, step_x;
    double   grey_scale0, grey_scaleF, step_grey;
    double   grey_shift0, grey_shiftF, step_grey_shift;
    int      tell;
    bool     apply;
    bool     mask_enabled, mean_in_mask;
    bool     usePowell, onlyShift;
    Mask_Params mask(INT_MASK);

    // Get parameters =======================================================
    try
    {
        fn1 = getParameter(argc, argv, "-i1");
        fn2 = getParameter(argc, argv, "-i2");
        getThreeDoubleParams(argc, argv, "-rot", rot0, rotF, step_rot, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-tilt", tilt0, tiltF, step_tilt, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-psi", psi0, psiF, step_psi, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-scale", scale0, scaleF, step_scale, 1, 1, 1);
        getThreeDoubleParams(argc, argv, "-grey_scale", grey_scale0, grey_scaleF,
                             step_grey, 1, 1, 1);
        getThreeDoubleParams(argc, argv, "-grey_shift", grey_shift0, grey_shiftF,
                             step_grey_shift, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-z", z0, zF, step_z, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-y", y0, yF, step_y, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-x", x0, xF, step_x, 0, 0, 1);

        mask_enabled = checkParameter(argc, argv, "-mask");
        mean_in_mask = checkParameter(argc, argv, "-mean_in_mask");
        if (mask_enabled)
            mask.read(argc, argv);

        usePowell = checkParameter(argc, argv, "-local");
        onlyShift = checkParameter(argc, argv, "-onlyShift");

        if (step_rot   == 0)
            step_rot = 1;
        if (step_tilt  == 0)
            step_tilt  = 1;
        if (step_psi   == 0)
            step_psi = 1;
        if (step_scale == 0)
            step_scale = 1;
        if (step_grey  == 0)
            step_grey  = 1;
        if (step_grey_shift  == 0)
            step_grey_shift  = 1;
        if (step_z     == 0)
            step_z = 1;
        if (step_y     == 0)
            step_y = 1;
        if (step_x     == 0)
            step_x = 1;
        tell = checkParameter(argc, argv, "-show_fit");
        apply = checkParameter(argc, argv, "-apply");

        if (checkParameter(argc, argv, "-covariance"))
        {
            prm.alignment_method = COVARIANCE;
        }
        else if (checkParameter(argc, argv, "-least_squares"))
        {
            prm.alignment_method = LEAST_SQUARES;
        }
        else
        {
            prm.alignment_method = COVARIANCE;
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE << std::endl;
        Usage(mask);
        exit(1);
    }

    // Main program =========================================================
    //#define DEBUG
    try
    {
        prm.V1.read(fn1);
        prm.V1().setXmippOrigin();
        prm.V2.read(fn2);
        prm.V2().setXmippOrigin();

        // Initialize best_fit
        double best_fit;
        Matrix1D<double> best_align(8);
        bool first = true;

        // Generate mask
        if (mask_enabled)
        {
            mask.generate_mask(prm.V1());
            prm.mask_ptr = &(mask.get_binary_mask());
        }
        else
            prm.mask_ptr = NULL;

        // Exhaustive search
        if (!usePowell)
        {
            // Count number of iterations
            int times = 1;
            if (!tell)
            {
                if (grey_scale0 != grey_scaleF)
                    times *= FLOOR(1 + (grey_scaleF - grey_scale0) / step_grey);
                if (grey_shift0 != grey_shiftF)
                    times *= FLOOR(1 + (grey_shiftF - grey_shift0) / step_grey_shift);
                if (rot0 != rotF)
                    times *= FLOOR(1 + (rotF - rot0) / step_rot);
                if (tilt0 != tiltF)
                    times *= FLOOR(1 + (tiltF - tilt0) / step_tilt);
                if (psi0 != psiF)
                    times *= FLOOR(1 + (psiF - psi0) / step_psi);
                if (scale0 != scaleF)
                    times *= FLOOR(1 + (scaleF - scale0) / step_scale);
                if (z0 != zF)
                    times *= FLOOR(1 + (zF - z0) / step_z);
                if (y0 != yF)
                    times *= FLOOR(1 + (yF - y0) / step_y);
                if (x0 != xF)
                    times *= FLOOR(1 + (xF - x0) / step_x);
                init_progress_bar(times);
            }
            else
                std::cout << "#grey_factor rot tilt psi scale z y x fitness\n";

            // Iterate
            int itime = 0;
            int step_time = CEIL((double)times / 60.0);
            Matrix1D<double> r(3);
            for (double grey_scale = grey_scale0; grey_scale <= grey_scaleF ; grey_scale += step_grey)
                for (double grey_shift = grey_shift0; grey_shift <= grey_shiftF ; grey_shift += step_grey_shift)
                    for (double rot = rot0; rot <= rotF ; rot += step_rot)
                        for (double tilt = tilt0; tilt <= tiltF ; tilt += step_tilt)
                            for (double psi = psi0; psi <= psiF ; psi += step_psi)
                                for (double scale = scale0; scale <= scaleF ; scale += step_scale)
                                    for (ZZ(r) = z0; ZZ(r) <= zF ; ZZ(r) += step_z)
                                        for (YY(r) = y0; YY(r) <= yF ; YY(r) += step_y)
                                            for (XX(r) = x0; XX(r) <= xF ; XX(r) += step_x)
                                            {
                                                // Form trial vector
                                                Matrix1D<double> trial(9);
                                                trial(0) = grey_scale;
                                                trial(1) = grey_shift;
                                                trial(2) = rot;
                                                trial(3) = tilt;
                                                trial(4) = psi;
                                                trial(5) = scale;
                                                trial(6) = ZZ(r);
                                                trial(7) = YY(r);
                                                trial(8) = XX(r);

                                                // Evaluate
                                                double fit =
                                                    fitness(MATRIX1D_ARRAY(trial));

                                                // The best?
                                                if (fit < best_fit || first)
                                                {
                                                    best_fit = fit;
                                                    best_align = trial;
                                                    first = false;
                                                }

                                                // Show fit
                                                if (tell)
                                                    std::cout << trial.transpose()
                                                    << fit << std::endl;
                                                else
                                                    if (++itime % step_time == 0)
                                                        progress_bar(itime);
                                            }
            if (!tell)
                progress_bar(times);
        }
        else
        {
            // Use Powell optimization
            Matrix1D<double> x(9), steps(9);
            double fitness;
            int iter;
            steps.initConstant(1);
            if (onlyShift)
                steps(0)=steps(1)=steps(2)=steps(3)=steps(4)=steps(5)=0;
            x(0)=(grey_scale0+grey_scaleF)/2;
            x(1)=(grey_shift0+grey_shiftF)/2;
            x(2)=(rot0+rotF)/2;
            x(3)=(tilt0+tiltF)/2;
            x(4)=(psi0+psiF)/2;
            x(5)=(scale0+scaleF)/2;
            x(6)=(z0+zF)/2;
            x(7)=(y0+yF)/2;
            x(8)=(x0+xF)/2;

            powellOptimizer(x,1,9,&wrapperFitness,NULL,
                            0.01,fitness,iter,steps,true);
            best_align=x;
            best_fit=fitness;
            first=false;
        }

        if (!first)
            std::cout << "The best correlation is for\n"
            << "Scale                  : " << best_align(5) << std::endl
            << "Translation (X,Y,Z)    : " << best_align(8)
            << " " << best_align(7) << " " << best_align(6)
            << std::endl
            << "Rotation (rot,tilt,psi): "
            << best_align(2) << " " << best_align(3) << " "
            << best_align(4) << std::endl
            << "Best grey scale       : " << best_align(0) << std::endl
            << "Best grey shift       : " << best_align(1) << std::endl
            << "Fitness value         : " << best_fit << std::endl;
        if (apply)
        {
            applyTransformation(prm.V2(),prm.Vaux(),MATRIX1D_ARRAY(best_align));
            prm.V2()=prm.Vaux();
            prm.V2.write();
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE << std::endl;
        return 1;
    }
    return 0;
}

void Usage(const Mask_Params &m)
{
    std::cerr << "Purpose: Align two volumes varying orientation, position and scale\n"
    << "Usage: align3D [options]\n"
    << "   -i1 <volume1>        : the first volume to align\n"
    << "   -i2 <volume2>        : the second one\n"
    << "  [-rot  <rot0>  <rotF>  <step_rot>  : in degrees\n"
    << "  [-tilt  <tilt0> <tiltF> <step_tilt> : in degrees\n"
    << "  [-psi  <psi0>  <psiF>  <step_psi>  : in degrees\n"
    << "  [-scale  <sc0>   <scF>   <step_sc>   : size scale margin\n"
    << "  [-grey_scale <sc0>   <scF>   <step_sc>   : grey scale margin\n"
    << "  [-grey_shift <sh0>   <shF>   <step_sh>   : grey shift margin\n"
    << "  [-z   <z0>  <zF>  <step_z>    : Z position in pixels\n"
    << "  [-y   <y0>  <yF>  <step_y>    : Y position in pixels\n"
    << "  [-x   <x0>  <xF>  <step_x>    : X position in pixels\n"
    << "  [-show_fit]         : Show fitness values\n"
    << "  [-apply]         : Apply best movement to -i2\n"
    << "  [-mean_in_mask]        : Use the means within the mask\n"
    << "  [-covariance]        : Covariance fitness criterion\n"
    << "  [-least_squares]        : LS fitness criterion\n"
    << "  [-local]                                 : Use local optimizer instead of\n"
    << "                                             exhaustive search\n"
    << "  [-onlyShift]                             : Only shift\n"
    ;
    m.usage();
    std::cout << std::endl;
}
