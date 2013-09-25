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

#include <data/args.h>
#include <data/mask.h>
#include <data/filters.h>
#include <data/geometry.h>
#include <data/symmetries.h>
#include <data/xmipp_program.h>
#include <data/xmipp_threads.h>

#include <cstdio>

// Prototypes
void globalThreadEvaluateSymmetry(ThreadArgument &thArg);
double evaluateSymmetryWrapper(double *p, void *prm);

class ProgVolumeFindSymmetry: public XmippProgram
{
public:
    int rot_sym;
    bool useSplines;
    FileName fn_input, fn_output;
    double   rot0,  rotF,  step_rot;
    double   tilt0, tiltF, step_tilt;
    double   z0, zF, step_z;
    bool     local;
    bool     helical;
    Mask mask_prm;
    int numberOfThreads;

    // Define parameters
    void defineParams()
    {
        addUsageLine("Find a symmetry rotational axis.");
        addUsageLine("+The output is of the form");
        addUsageLine("+Symmetry axis (rot,tilt)= 10 20 -->    0.33682   0.059391    0.93969",true);
        addUsageLine("+The angles represent the rot and tilt angles of the symmetry axis ");
        addUsageLine("+see [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Conventions#Euler_Angles][the note on Euler angles]] ");
        addUsageLine("+convention in Xmipp). The rest of the numbers is the axis itself in X, Y, Z ");
        addUsageLine("+coordinates. In this way, the symmetry axis is a line that passes through the ");
        addUsageLine("+center of the volume and whose direction is this vector.");
        addUsageLine("+It is important that the volume is correctly centered.");
        addSeeAlsoLine("volume_center, transform_geometry");
        addParamsLine(" -i <volumeFile>                : Volume to process");
        addParamsLine("[-o+ <file=\"\">]               : Metadata with the orientation of the symmetry axis");
        addParamsLine("                                : In helical mode, there is a file called output.xmp ");
        addParamsLine("                                : with the correlation map (vertical axis is the rotation, ");
        addParamsLine("                                : horizontal axis is the translation)");
        addParamsLine("[--thr <N=1>]                   : Number of threads");
        addParamsLine("--sym <mode>                    : Symmetry mode");
        addParamsLine("    where <mode>");
        addParamsLine("           rot <n>              : Order of the rotational axis");
        addParamsLine("           helical              : Helical symmetry.");
        addParamsLine("==Locate rotational axis==");
        addParamsLine("[--rot  <rot0=0>  <rotF=355> <step=5>]: Limits for rotational angle search");
        addParamsLine("[--tilt <tilt0=0> <tiltF=90> <step=5>]: Limits for tilt angle search");
        addParamsLine("[--localRot <rot0> <tilt0>]     : Perform a local search around this angle");
        addParamsLine("[--useSplines+]                 : Use cubic B-Splines for the interpolations");
        addParamsLine("==Locate helical parameters==");
        addParamsLine("[-z <z0=1> <zF=10> <zstep=0.5>]          : Search space for the shift in Z");
        addParamsLine("[--rotHelical <rot0=0> <rotF=357> <step=3>]: Search space for rotation around Z");
        addParamsLine("[--localHelical <z> <rot> <rot0=0>]  : Perform a local search around this angle and shift");
        mask_prm.defineParams(this,INT_MASK,NULL,"Restrict the comparison to the mask area.",true);
        addExampleLine("A typical application for a rotational symmetry axis is ",false);
        addExampleLine("xmipp_volume_center -i volume.vol");
        addExampleLine("xmipp_volume_find_symmetry -i volume.vol --sym rot 3");
        addExampleLine("Presume the symmetry axis is in rot=20, tilt=10. To align vertically the axis use",false);
        addExampleLine("xmipp_transform_geometry -i volume.vol --rotate euler 20 10 0");
        addExampleLine("For locating the helical parameters use",false);
        addExampleLine("xmipp_volume_find_symmetry -i volume --sym helical -z 0 6 1 --mask circular -32 --thr 2 -o parameters.xmd");
    }

    // Read parameters
    void readParams()
    {
        fn_input = getParam("-i");
        fn_output = getParam("-o");
        String mode;
        mode=getParam("--sym");
        if (mode=="helical")
            helical=true;
        else
        {
            helical=false;
            rot_sym = getIntParam("--sym",1);
        }
        useSplines = checkParam("--useSplines");
        if (helical)
        {
            local=checkParam("--localHelical");
            if (local)
            {
                z0=getDoubleParam("--localHelical",0);
                rot0=getDoubleParam("--localHelical",1);
            }
            else
            {
                z0=getDoubleParam("-z",0);
                zF=getDoubleParam("-z",1);
                step_z=getDoubleParam("-z",2);
                rot0=getDoubleParam("--rotHelical",0);
                rotF=getDoubleParam("--rotHelical",1);
                step_rot=getDoubleParam("--rotHelical",2);
            }
        }
        else
        {
            local=checkParam("--localRot");
            if (local)
            {
                rot0=getDoubleParam("--localRot",0);
                tilt0=getDoubleParam("--localRot",1);
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
        }
        mask_prm.allowed_data_types = INT_MASK;
        if (checkParam("--mask"))
            mask_prm.readParams(this);
        numberOfThreads=getIntParam("--thr");
    }

    void threadEvaluateSymmetry(int thrId)
    {
        Matrix1D<double> p(2);
        DIRECT_A1D_ELEM(vbest_corr,thrId) = 0;
        DIRECT_A1D_ELEM(vbest_rot,thrId)  = 0;
        if (helical)
            DIRECT_A1D_ELEM(vbest_z,thrId) = 0;
        else
            DIRECT_A1D_ELEM(vbest_tilt,thrId) = 0;
        size_t first, last;
        if (thrId==0)
            init_progress_bar(rotVector.size());
        while (td->getTasks(first, last))
        {
            for (size_t i=first; i<=last; i++)
            {
                XX(p)=rotVector[i];
                if (helical)
                    YY(p)=zVector[i];
                else
                    YY(p)=tiltVector[i];
                double corr=-evaluateSymmetry(MATRIX1D_ARRAY(p)-1);
                if (helical)
                    DIRECT_MULTIDIM_ELEM(helicalCorrelation(),i)=corr;
                if (corr > DIRECT_A1D_ELEM(vbest_corr,thrId))
                {
                    DIRECT_A1D_ELEM(vbest_corr,thrId) = corr;
                    DIRECT_A1D_ELEM(vbest_rot,thrId) = XX(p);
                    if (helical)
                        DIRECT_A1D_ELEM(vbest_z,thrId) = YY(p);
                    else
                        DIRECT_A1D_ELEM(vbest_tilt,thrId) = YY(p);
                }
            }
            if (thrId==0)
                progress_bar(last);
        }
        if (thrId==0)
            progress_bar(rotVector.size());
    }

    void run()
    {
    	// Read input volume
        volume.read(fn_input);
        volume().setXmippOrigin();
        mask_prm.generate_mask(volume());
        double best_corr, best_rot, best_tilt, best_z;
        td=NULL;

        if (!helical)
        {
            // Look for the rotational symmetry axis
            if (!local)
            {
                if (verbose>0)
                    std::cerr << "Searching symmetry axis ...\n";
                for (double rot = rot0; rot <= rotF; rot += step_rot)
                    for (double tilt = tilt0; tilt <= tiltF; tilt += step_tilt)
                    {
                        rotVector.push_back(rot);
                        tiltVector.push_back(tilt);
                    }
                vbest_corr.resizeNoCopy(numberOfThreads);
                vbest_corr.initConstant(-1e38);
                vbest_rot.initZeros(numberOfThreads);
                vbest_tilt.initZeros(numberOfThreads);
                td = new ThreadTaskDistributor(rotVector.size(), 5);
                ThreadManager thMgr(numberOfThreads,this);
                thMgr.run(globalThreadEvaluateSymmetry);
                best_corr=-1e38;
                FOR_ALL_ELEMENTS_IN_ARRAY1D(vbest_corr)
                if (vbest_corr(i)>best_corr)
                {
                    best_corr=vbest_corr(i);
                    best_rot=vbest_rot(i);
                    best_tilt=vbest_tilt(i);
                }
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
            Matrix2D<double> Euler;
            Matrix1D<double> sym_axis;
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
                MD.setValue(MDL_ANGLE_ROT,best_rot,id);
                MD.setValue(MDL_ANGLE_TILT,best_tilt,id);
                MD.setValue(MDL_DIRECTION,direction,id);
                MD.write(fn_output);
            }
        }
        else
        {
            // If helical
            if (!local)
            {
                int ydim=0, xdim=0;
				for (double rot = rot0; rot <= rotF; rot += step_rot)
                {
                    ydim++;
					for (double z = z0; z <= zF; z += step_z)
					{
						if (ydim==1)
							xdim++;
						rotVector.push_back(rot);
						zVector.push_back(z);
					}
                }
                vbest_corr.resizeNoCopy(numberOfThreads);
                vbest_corr.initConstant(-1e38);
                vbest_rot.initZeros(numberOfThreads);
                vbest_z.initZeros(numberOfThreads);
                helicalCorrelation().initZeros(ydim,xdim);
                td = new ThreadTaskDistributor(rotVector.size(), 5);
                ThreadManager thMgr(numberOfThreads,this);
                thMgr.run(globalThreadEvaluateSymmetry);
                best_corr=-1e38;
                FOR_ALL_ELEMENTS_IN_ARRAY1D(vbest_corr)
                if (vbest_corr(i)>best_corr)
                {
                    best_corr=vbest_corr(i);
                    best_rot=vbest_rot(i);
                    best_z=vbest_z(i);
                }
            }
            else
            {
                Matrix1D<double> p(2), steps(2);
                p(0)=rot0;
                p(1)=z0;
                double fitness;
                int iter;
                steps.initConstant(1);
                powellOptimizer(p,1,2,&evaluateSymmetryWrapper,this,0.01,
                                fitness,iter,steps,true);
                best_rot=p(0);
                best_z=p(1);

            }
            if (verbose!=0)
                std::cout << "Symmetry parameters (rot,z)= " << best_rot << " " << best_z << std::endl;
            if (fn_output!="")
            {
                MetaData MD;
                size_t id=MD.addObject();
                MD.setValue(MDL_ANGLE_ROT,best_rot,id);
                MD.setValue(MDL_SHIFT_Z,best_z,id);
                MD.write(fn_output);
                if (!local)
                    helicalCorrelation.write(fn_output.removeAllExtensions()+".xmp");
            }
        }
        delete td;
    }

    Image<double> volume;
    std::vector<double> rotVector, tiltVector, zVector;
    ThreadTaskDistributor * td;
    MultidimArray<double> vbest_corr, vbest_rot, vbest_tilt, vbest_z;
    Image<double> helicalCorrelation;

    /* Evaluate symmetry ------------------------------------------------------- */
    double evaluateSymmetry(double *p)
    {
        MultidimArray<double> volume_sym, volume_aux;
        const MultidimArray<double> &mVolume=volume();
        Matrix2D<double> sym_matrix;
        if (!helical)
        {
            Matrix2D<double> Euler;
            Matrix1D<double> sym_axis;

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
                    applyGeometry(BSPLINE3, volume_aux, mVolume, sym_matrix,
                                  IS_NOT_INV, DONT_WRAP);
                else
                    applyGeometry(LINEAR, volume_aux, mVolume, sym_matrix,
                                  IS_NOT_INV, DONT_WRAP);
                volume_sym += volume_aux;
            }
            return -correlationIndex(mVolume, volume_sym, &mask_prm.get_binary_mask());
        }
        else
        {
            double rotHelical=DEG2RAD(p[1]);
            double zHelical=p[2];
            if (zHelical<0 || zHelical>ZSIZE(mVolume)*0.4)
            	return 1e38;
            symmetry_Helical(volume_sym, mVolume, zHelical, rotHelical, 0, &mask_prm.get_binary_mask());
			double corr=correlationIndex(mVolume, volume_sym, &mask_prm.get_binary_mask());
//#define DEBUG
#ifdef DEBUG

                Image<double> save;
                save()=mVolume;
                save.write("PPPvol.vol");
                save()=volume_sym;
                save.write("PPPvolsym.vol");
                typeCast(mask_prm.get_binary_mask(),save());
                save.write("PPPmask.vol");
                std::cout << p[1] << " " << p[2] << " " << p[3] << " -> " << corr << std::endl;
                char c;
                std::cin >> c;
#endif
			return -corr;
        }
    }
};

void globalThreadEvaluateSymmetry(ThreadArgument &thArg)
{
    ((ProgVolumeFindSymmetry*)thArg.workClass)->threadEvaluateSymmetry(thArg.thread_id);
}

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
