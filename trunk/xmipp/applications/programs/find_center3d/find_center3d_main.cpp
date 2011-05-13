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

#include <data/image.h>
#include <data/args.h>
#include <data/mask.h>
#include <data/filters.h>
#include <data/geometry.h>
#include <data/program.h>

#include <cstdio>

// Prototype
double evaluateSymmetryWrapper(double *p, void *prm);

class ProgVolumeFindSymmetry: public XmippProgram
{
public:
    int rot_sym;
    bool useSplines;
    FileName fn_input, fn_output;
    double   rot0,  rotF,  step_rot;
    double   tilt0, tiltF, step_tilt;
    bool     local;
    Mask mask_prm;

    // Define parameters
    void defineParams()
    {
        addUsageLine("Finds a symmetry rotational axis.");
        addUsageLine("+It is important that the volume is correctly centered");
        addSeeAlsoLine("volume_center");
        addParamsLine(" -i <volumeFile>               : Volume to process");
        addParamsLine("[-o+ <file=\"\">]              : Metadata with the orientation of the symmetry axis");
        addParamsLine("--rot_sym <n>                  : Order of the rotational axis");
        addParamsLine("[--rot  <rot0=0>  <rotF=355> <step=5>]: Limits for rotational angle search");
        addParamsLine("[--tilt <tilt0=0> <tiltF=90> <step=5>]: Limits for tilt angle search");
        addParamsLine("[--local <rot0> <tilt0>]       : Perform a local search around this angle");
        addParamsLine("[--useSplines+]                 : Use cubic B-Splines for the interpolations");
        mask_prm.defineParams(this,INT_MASK,NULL,"Restrict the comparison to the mask area.");
    }

    // Read parameters
    void readParams()
    {
        fn_input = getParam("-i");
        fn_output = getParam("-o");
        rot_sym = getIntParam("--rot_sym");
        useSplines = checkParam("--useSplines");
        local = checkParam("--local");
        if (local)
        {
            rot0=getDoubleParam("--local",0);
            tilt0=getDoubleParam("--local",1);
        }
        else
        {
            rot0=getDoubleParam("--rot",0);
            rotF=getDoubleParam("--rot",1);
            step_rot=getDoubleParam("--rot",2);
            tilt0=getDoubleParam("--tilt",0);
            tiltF=getDoubleParam("--tilt",1);
            step_tilt=getDoubleParam("--tilt",2);
        }
        mask_prm.allowed_data_types = INT_MASK;
        if (checkParam("--mask"))
            mask_prm.readParams(this);
    }

    void run()
    {
        // Read input volume
        volume.read(fn_input);
        volume().setXmippOrigin();
        mask_prm.generate_mask(volume());

        // Look for the rotational symmetry axis
        double best_corr = 0, best_rot = rot0 - step_rot, best_tilt = tilt0 - step_tilt;
        if (!local)
        {
            int maxsteps = FLOOR((rotF - rot0) / step_rot+1) *
                           FLOOR((tiltF - tilt0) / step_tilt +1);
            if (verbose>0)
            {
                std::cerr << "Searching symmetry axis ...\n";
                init_progress_bar(maxsteps);
            }
            int i = 0;
            Matrix1D<double> p(2);
            for (double rot = rot0; rot <= rotF; rot += step_rot)
                for (double tilt = tilt0; tilt <= tiltF; tilt += step_tilt)
                {
                    p(0)=rot;
                    p(1)=tilt;
                    double corr=-evaluateSymmetry(MATRIX1D_ARRAY(p)-1);
                    if (corr > best_corr)
                    {
                        best_corr = corr;
                        best_rot = rot;
                        best_tilt = tilt;
                        if (verbose!=0)
                            std::cout << "rot=" << rot << " tilt=" << tilt
                            << " corr=" << corr << std::endl;
                    }

                    // progress bar
                    if ((i++) % XMIPP_MAX(maxsteps / 60, 1) == 0 && verbose!=0)
                        progress_bar(i);
                }
            if (verbose!=0)
                progress_bar(maxsteps);
        }
        else
        {
            Matrix1D<double> p(2), steps(2);
            p(0)=rot0;
            p(1)=tilt0;
            double fitness;
            int iter;
            steps.initConstant(1);
            powellOptimizer(p,1,2,&evaluateSymmetryWrapper,this,0.01,
                            fitness,iter,steps,true);
            best_rot=p(0);
            best_tilt=p(1);
        }
        Euler_angles2matrix(best_rot, best_tilt, 0, Euler);
        Euler.getRow(2, sym_axis);
        if (verbose!=0)
        	std::cout << "Symmetry axis (rot,tilt)= " << best_rot << " "
        	<< best_tilt << " --> " << sym_axis << std::endl;
        if (fn_output!="")
        {
        	std::vector<double> direction;
        	direction.push_back(XX(sym_axis));
        	direction.push_back(YY(sym_axis));
        	direction.push_back(ZZ(sym_axis));

        	MetaData MD;
        	size_t id=MD.addObject();
        	MD.setValue(MDL_ANGLEROT,best_rot,id);
        	MD.setValue(MDL_ANGLETILT,best_tilt,id);
        	MD.setValue(MDL_DIRECTION,direction,id);
        	MD.write(fn_output);
        }
    }

    Image<double> volume;
    MultidimArray<double> volume_sym, volume_aux;
    Matrix2D<double> Euler, sym_matrix;
    Matrix1D<double> sym_axis;

    /* Evaluate symmetry ------------------------------------------------------- */
    double evaluateSymmetry(double *p)
    {
        double rot=p[1];
        double tilt=p[2];

        // Compute symmetry axis
        Euler_angles2matrix(rot, tilt, 0, Euler);
        Euler.getRow(2, sym_axis);
        sym_axis.selfTranspose();

        // Symmetrize along this axis
        volume_sym = volume();
        for (int n = 1; n < rot_sym; n++)
        {
            rotation3DMatrix(360.0 / rot_sym * n, sym_axis, sym_matrix);
            if (useSplines)
                applyGeometry(BSPLINE3, volume_aux, volume(), sym_matrix,
                              IS_NOT_INV, DONT_WRAP);
            else
                applyGeometry(LINEAR, volume_aux, volume(), sym_matrix,
                              IS_NOT_INV, DONT_WRAP);
            volume_sym += volume_aux;
        }

        // Measure correlation
        return -correlationIndex(volume(), volume_sym, &mask_prm.get_binary_mask());
    }
};

double evaluateSymmetryWrapper(double *p, void *prm)
{
	return ((ProgVolumeFindSymmetry*)prm)->evaluateSymmetry(p);
}

int main(int argc, char **argv)
{
    ProgVolumeFindSymmetry prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
