/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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

#include <data/volume.h>
#include <data/args.h>
#include <data/mask.h>
#include <data/filters.h>
#include <data/geometry.h>

#include <cstdio>

void Usage();
double evaluateSymmetry(double *p, void *prm);
VolumeXmipp volume, volume_sym, volume_aux;
Mask_Params mask_prm(INT_MASK);
int rot_sym;
bool useSplines;

int main(int argc, char **argv)
{
    FileName fn_input, fn_output;
    double   rot0,  rotF,  step_rot;
    double   tilt0, tiltF, step_tilt;
    bool     centerVolume;
	bool     local;
    bool     show;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = getParameter(argc, argv, "-i");
        fn_output = getParameter(argc, argv, "-o","");
        rot_sym = textToInteger(getParameter(argc, argv, "-rot_sym", "0"));
        centerVolume = checkParameter(argc, argv, "-center_volume");
        useSplines = checkParameter(argc, argv, "-useSplines");
        local = checkParameter(argc, argv, "-local");
        show = checkParameter(argc, argv, "-show");
        int i;
        if ((i = paremeterPosition(argc, argv, "-rot")) != -1)
        {
			if (!local)
			{
            	if (i + 3 >= argc)
                	REPORT_ERROR(1, "findcenter3D: Not enough parameters behind -rot");
            	rot0    = textToFloat(argv[i+1]);
            	rotF    = textToFloat(argv[i+2]);
            	step_rot = textToFloat(argv[i+3]);
			}
			else
				rot0=textToFloat(getParameter(argc, argv, "-rot"));
        }
        else
        {
            rot0 = 0;
            rotF = 355;
            step_rot = 5;
        }
        if ((i = paremeterPosition(argc, argv, "-tilt")) != -1)
        {
			if (!local)
			{
            	if (i + 3 >= argc)
                	REPORT_ERROR(1, "findcenter3D: Not enough parameters behind -tilt");
            	tilt0    = textToFloat(argv[i+1]);
            	tiltF    = textToFloat(argv[i+2]);
            	step_tilt = textToFloat(argv[i+3]);
			}
			else
				tilt0=textToFloat(getParameter(argc, argv, "-tilt"));
        }
        else
        {
            tilt0 = 0;
            tiltF = 90;
            step_tilt = 5;
        }
        mask_prm.read(argc, argv);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        mask_prm.usage();
        exit(1);
    }

    // Find Center and symmetry elements ------------------------------------
    try
    {
        // Read input volumes
        volume.read(fn_input);
        volume().setXmippOrigin();
        mask_prm.generate_3Dmask(volume());

        // Compute center of mass
        Matrix1D<double> centerOfMass;
        volume().centerOfMass(centerOfMass, &mask_prm.get_binary_mask3D());
        std::cout << "Center of mass (X,Y,Z)= " << centerOfMass.transpose() << std::endl;

        // Move origin to that center of mass
	if (useSplines)
	    volume().selfTranslateBSpline(3,-centerOfMass, DONT_WRAP);
	else
            volume().selfTranslate(-centerOfMass, DONT_WRAP);
        if (centerVolume) volume.write();

        // Look for the rotational symmetry axis
        if (rot_sym > 1)
        {
            double best_corr = 0, best_rot = rot0 - step_rot, best_tilt = tilt0 - step_tilt;
			if (!local)
			{
            	int maxsteps = FLOOR((rotF - rot0) / step_rot+1) *
	                	   FLOOR((tiltF - tilt0) / step_tilt +1);
            	std::cerr << "Searching symmetry axis ...\n";
            	if (!show) init_progress_bar(maxsteps);
            	int i = 0;
            	for (double rot = rot0; rot <= rotF; rot += step_rot)
                	for (double tilt = tilt0; tilt <= tiltF; tilt += step_tilt)
                	{
				    	Matrix1D<double> p(2);
						p(0)=rot;
						p(1)=tilt;
                    	double corr=-evaluateSymmetry(MULTIDIM_ARRAY(p)-1,NULL);
                    	if (corr > best_corr)
                    	{
                        	best_corr = corr;
                        	best_rot = rot;
                        	best_tilt = tilt;
                    	}

                    	// progress bar
                    	if ((i++) % XMIPP_MAX(maxsteps / 60, 1) == 0 && !show)
                            progress_bar(i);
                        if (show)
                            std::cout << "rot=" << rot << " tilt=" << tilt
                                      << " corr=" << corr << std::endl;
                	}
            	if (!show) progress_bar(maxsteps);
			}
            else
            {
                Matrix1D<double> p(2), steps(2);
				p(0)=rot0;
				p(1)=tilt0;
                double fitness;
                int iter;
                steps.initConstant(1);
                powellOptimizer(p,1,2,&evaluateSymmetry,NULL,0.01,
                   fitness,iter,steps,true);
                best_rot=p(0);
                best_tilt=p(1);
            }
            std::ofstream fh_out;
            if (fn_output!="")
            {
                fh_out.open(fn_output.c_str());
                if (!fh_out)
                    REPORT_ERROR(1,(std::string)"Cannot open "+fn_output+
                        " for output");
            }
            std::cout << "Symmetry axis (rot,tilt)= " << best_rot << " "
                      << best_tilt << " --> ";
            if (fn_output!="")
                fh_out << "Symmetry axis (rot,tilt)= " << best_rot << " "
                       << best_tilt << " --> ";
            Matrix2D<double> Euler;
            Matrix1D<double> sym_axis(3);
            Euler_angles2matrix(best_rot, best_tilt, 0, Euler);
            Euler.getRow(2, sym_axis);
            std::cout << sym_axis << std::endl;
            if (fn_output!="")
            {
                fh_out << sym_axis << std::endl;
                fh_out.close();
            }
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cerr << "Purpose:\n";
    std::cerr << "    Finds the 3D center of mass within a mask\n"
              << "    and a symmetry rotational axis passing through that center\n";

    std::cerr << "Usage: findcenter3D [options]" << std::endl
              << "    -i <volume>                         : volume to process\n"
              << "   [-center_volume]                     : save the centered volume\n"
              << "   [-useSplines]                        : use BSplines(3) for the interpolations\n"
              << "   [-rot_sym <n=0>]                     : order of the rotational axis\n"
              << "   [-rot  <rot0=0>  <rotF=355> <step=5>]: limits for rotational axis\n"
              << "   [-tilt <tilt0=0> <tiltF=90> <step=5>]: limits for rotational axis\n"
              << "   [-local]                             : perform a local search\n"
              << "                                          in this case use -rot rot0 -tilt tilt0\n"
              << "   [-show]                              : show correlation for all trials\n"
    ;
}

/* Evaluate symmetry ------------------------------------------------------- */
double evaluateSymmetry(double *p, void *prm)
{
    double rot=p[1];
	double tilt=p[2];

	// Compute symmetry axis
	Matrix2D<double> Euler;
	Euler_angles2matrix(rot, tilt, 0, Euler);
	Matrix1D<double> sym_axis(3);
	Euler.getRow(2, sym_axis);
	sym_axis.selfTranspose();

	// Symmetrize along this axis
	volume_sym() = volume();
	for (int n = 1; n < rot_sym; n++)
	{
		Matrix2D<double> sym_matrix;
		sym_matrix = rotation3DMatrix(360.0 / rot_sym * n, sym_axis);
		if (useSplines)
			applyGeometryBSpline(volume_aux(), sym_matrix, volume(),
				3, IS_NOT_INV, DONT_WRAP);
		else
			applyGeometry(volume_aux(), sym_matrix, volume(), IS_NOT_INV,
				DONT_WRAP);
		volume_sym() += volume_aux();
	}

	// Measure correlation
	return -correlation_index(volume(), volume_sym(),
		&mask_prm.get_binary_mask3D());
}
