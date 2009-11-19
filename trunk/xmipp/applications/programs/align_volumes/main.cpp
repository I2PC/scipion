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

#include <data/args.h>
#include <data/volume.h>
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
       
       VolumeXmipp V1;
       VolumeXmipp V2;
       VolumeXmipp Vaux;
       const Matrix3D<int> *mask_ptr;
};

// Global parameters needed by fitness ------------------------------------
AlignParams prm;

// Apply transformation ---------------------------------------------------
void applyTransformation(const Matrix3D<double> &V2, Matrix3D<double> &Vaux,
   double *p) 
{
    Matrix1D<double> r(3);
    Matrix2D<double> A;

    double grey  = p[0];
    double rot   = p[1];
    double tilt  = p[2];
    double psi   = p[3];
    double scale = p[4];
    ZZ(r)        = p[5];
    YY(r)        = p[6];
    XX(r)        = p[7];

    Euler_angles2matrix(rot, tilt, psi, A);
    A.resize(4, 4);
    A(3, 3) = 1;
    A = A * translation3DMatrix(r);
    A = A * scale3DMatrix(vectorR3(scale, scale, scale));
    
    applyGeometryBSpline(Vaux, A, V2, 3, IS_NOT_INV, WRAP);
    Vaux*=grey;
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
    int      tell;
    bool     apply;
    bool     mask_enabled, mean_in_mask;
    bool     usePowell;
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
        getThreeDoubleParams(argc, argv, "-z", z0, zF, step_z, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-y", y0, yF, step_y, 0, 0, 1);
        getThreeDoubleParams(argc, argv, "-x", x0, xF, step_x, 0, 0, 1);

        mask_enabled = checkParameter(argc, argv, "-mask");
        mean_in_mask = checkParameter(argc, argv, "-mean_in_mask");
        if (mask_enabled)
            mask.read(argc, argv);

	usePowell = checkParameter(argc, argv, "-local");

        if (step_rot   == 0) step_rot	= 1;
        if (step_tilt  == 0) step_tilt  = 1;
        if (step_psi   == 0) step_psi	= 1;
        if (step_scale == 0) step_scale = 1;
        if (step_grey  == 0) step_grey  = 1;
        if (step_z     == 0) step_z	= 1;
        if (step_y     == 0) step_y	= 1;
        if (step_x     == 0) step_x	= 1;
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
            mask.generate_3Dmask(prm.V1());
            prm.mask_ptr = &(mask.get_binary_mask3D());
        }
        else prm.mask_ptr = NULL;

        // Exhaustive search
    	if (!usePowell)
	{
            // Count number of iterations
            int times = 1;
            if (!tell)
            {
        	if (grey_scale0 != grey_scaleF)
                    times *= FLOOR(1 + (grey_scaleF - grey_scale0) / step_grey);
        	if (rot0 != rotF) times *= FLOOR(1 + (rotF - rot0) / step_rot);
        	if (tilt0 != tiltF) times *= FLOOR(1 + (tiltF - tilt0) / step_tilt);
        	if (psi0 != psiF) times *= FLOOR(1 + (psiF - psi0) / step_psi);
        	if (scale0 != scaleF) times *= FLOOR(1 + (scaleF - scale0) / step_scale);
        	if (z0 != zF) times *= FLOOR(1 + (zF - z0) / step_z);
        	if (y0 != yF) times *= FLOOR(1 + (yF - y0) / step_y);
        	if (x0 != xF) times *= FLOOR(1 + (xF - x0) / step_x);
        	init_progress_bar(times);
            }
            else
        	std::cout << "#grey_factor rot tilt psi scale z y x fitness\n";

            // Iterate
            int itime = 0;
            int step_time = CEIL((double)times / 60.0);
	    Matrix1D<double> r(3);
            for (double grey = grey_scale0; grey <= grey_scaleF ; grey += step_grey)
        	for (double rot = rot0; rot <= rotF ; rot += step_rot)
                    for (double tilt = tilt0; tilt <= tiltF ; tilt += step_tilt)
                	for (double psi = psi0; psi <= psiF ; psi += step_psi)
                            for (double scale = scale0; scale <= scaleF ; scale += step_scale)
                        	for (ZZ(r) = z0; ZZ(r) <= zF ; ZZ(r) += step_z)
                                    for (YY(r) = y0; YY(r) <= yF ; YY(r) += step_y)
                                	for (XX(r) = x0; XX(r) <= xF ; XX(r) += step_x)
                                	{
					    // Form trial vector
					    Matrix1D<double> trial(8);
					    trial(0) = grey;
					    trial(1) = rot;
					    trial(2) = tilt;
					    trial(3) = psi;
					    trial(4) = scale;
					    trial(5) = ZZ(r);
					    trial(6) = YY(r);
					    trial(7) = XX(r);
					    
					    // Evaluate
					    double fit =
					       fitness(MULTIDIM_ARRAY(trial));

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
            if (!tell) progress_bar(times);
	}
	else
	{
	    // Use Powell optimization
	    Matrix1D<double> x(8), steps(8);
	    double fitness;
	    int iter;
	    steps.initConstant(1);
	    x(0)=(grey_scale0+grey_scaleF)/2;
	    x(1)=(rot0+rotF)/2;
	    x(2)=(tilt0+tiltF)/2;
	    x(3)=(psi0+psiF)/2;
	    x(4)=(scale0+scaleF)/2;
	    x(5)=(z0+zF)/2;
	    x(6)=(y0+yF)/2;
	    x(7)=(x0+xF)/2;
	    
	    powellOptimizer(x,1,8,&wrapperFitness,NULL,
                0.01,fitness,iter,steps,true);
	    best_align=x;
            first=false;
	}

        if (!first)
             std::cout << "The best correlation is for\n"
                       << "Scale                  : " << best_align(4) << std::endl
                       << "Translation (X,Y,Z)    : " << best_align(7)
		       << " " << best_align(6) << " " << best_align(5)
		       << std::endl
            	       << "Rotation (rot,tilt,psi): "
            	       << best_align(1) << " " << best_align(2) << " "
		       << best_align(3) << std::endl
            	       << "Best grey factor       : " << best_align(0) << std::endl
            	       << "Fitness value          : " << best_fit << std::endl;
        if (apply)
        {
	    applyTransformation(prm.V2(),prm.Vaux(),MULTIDIM_ARRAY(best_align));
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
	      << "   -i1 <volume1>			     : the first volume to align\n"
	      << "   -i2 <volume2>			     : the second one\n"
	      << "  [-rot	 <rot0>  <rotF>  <step_rot>  : in degrees\n"
	      << "  [-tilt	 <tilt0> <tiltF> <step_tilt> : in degrees\n"
	      << "  [-psi	 <psi0>  <psiF>  <step_psi>  : in degrees\n"
	      << "  [-scale	 <sc0>   <scF>   <step_sc>   : size scale margin\n"
	      << "  [-grey_scale <sc0>   <scF>   <step_sc>   : grey scale margin\n"
	      << "  [-z 	 <z0>	 <zF>	 <step_z>    : Z position in pixels\n"
	      << "  [-y 	 <y0>	 <yF>	 <step_y>    : Y position in pixels\n"
	      << "  [-x 	 <x0>	 <xF>	 <step_x>    : X position in pixels\n"
	      << "  [-show_fit] 			     : Show fitness values\n"
	      << "  [-apply]				     : Apply best movement to -i2\n"
	      << "  [-mean_in_mask]			     : Use the means within the mask\n"
	      << "  [-covariance]			     : Covariance fitness criterion\n"
	      << "  [-least_squares]			     : LS fitness criterion\n"
	      << "  [-local]                                 : Use local optimizer instead of\n"
	      << "                                             exhaustive search\n"
    ;
    m.usage();
    std::cout << std::endl;
}
