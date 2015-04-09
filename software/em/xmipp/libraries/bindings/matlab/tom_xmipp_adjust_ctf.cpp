/*=================================================================
 *
 * tom_xmipp_adjust_ctf is a wrapper to xmipp_adjust_ctf
 *
 * The calling syntax is:
 *
 *		psd = tom_xmipp_adjust_ctf_wrapper(image,min_freq,max_freq,ctfmodelSize,Ca,Tm,voltage,Cs,DeltafU,f1,f2,enhanced_weight)
 *
 * Electron Tomography toolbox of the
 * Max-Planck-Institute for Biochemistry
 * Dept. Molecular Structural Biology
 * 82152 Martinsried, Germany
 * http://www.biochem.mpg.de
 *
 * and
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 *
 * created: 27/08/2007
 * by: Andreas Korinek & Carlos Oscar Sorzano
 *
 *=================================================================*/

/*xmipp includes */
#include <reconstruction/ctf_estimate_from_psd.h>
#include "tom_xmipp_helpers.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{

    ProgCTFEstimateFromPSD adjustParams;
    adjustParams.show_optimization=false;
    getMatrix2D(prhs[0],adjustParams.ctftomodel());
    adjustParams.min_freq=(double) mxGetScalar(prhs[1]);
    adjustParams.max_freq=(double) mxGetScalar(prhs[2]);
    adjustParams.defocus_range=30000;
    adjustParams.ctfmodelSize=(int) mxGetScalar(prhs[3]);
    adjustParams.initial_ctfmodel.Ca=(double) mxGetScalar(prhs[4]);

    adjustParams.initial_ctfmodel.enable_CTF = 
       adjustParams.initial_ctfmodel.enable_CTFnoise = true;
    adjustParams.Tm =
       adjustParams.initial_ctfmodel.Tm = (double) mxGetScalar(prhs[5]);
    adjustParams.initial_ctfmodel.kV = (double) mxGetScalar(prhs[6]);
    adjustParams.initial_ctfmodel.Cs = (double) mxGetScalar(prhs[7]);
    adjustParams.initial_ctfmodel.DeltafU = 
       adjustParams.initial_ctfmodel.DeltafV =
       (double) mxGetScalar(prhs[8]);
    
    adjustParams.f1=(double) mxGetScalar(prhs[9]);
    adjustParams.f2=(double) mxGetScalar(prhs[10]);
    adjustParams.enhanced_weight=(double) mxGetScalar(prhs[11]);

    CTFDescription ctfmodel;
    try 
    {
        ROUT_Adjust_CTF(adjustParams,ctfmodel,false);
    }
    catch (XmippError Xe)
    {
        mexErrMsgTxt(Xe.msg.c_str());
    }
    
    const char *field_names[] = {"DeltafU","DeltafV","AzimuthalAngle",
       "kV","K","Cs","Ca","espr","ispr","alpha","DeltaF","DeltaR","Q0",
       "base_line","sqrt_K","sqU","sqV","sqrt_angle","gaussian_K",
       "sigmaU","sigmaV","gaussian_angle","cU","cV","gaussian_K2",
       "sigmaU2","sigmaV2","gaussian_angle2","cU2","cV2",
       "CTFmodelhalf","CTFmodelquadrant","zeros","Tm"};

    mwSize dims[2] = {1, 1};
    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);
    
    mxArray *field1 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field1) = ctfmodel.DeltafU;
    mxSetField(plhs[0],0,"DeltafU",field1);

    mxArray *field2 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field2) = ctfmodel.DeltafV;
    mxSetField(plhs[0],0,"DeltafV",field2);

    mxArray *field3 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field3) = ctfmodel.azimuthal_angle;
    mxSetField(plhs[0],0,"AzimuthalAngle",field3);

    mxArray *field4 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field4) = ctfmodel.kV;
    mxSetField(plhs[0],0,"kV",field4);
    
    mxArray *field5 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field5) = ctfmodel.K;
    mxSetField(plhs[0],0,"K",field5);

    mxArray *field6 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field6) = ctfmodel.Cs;
    mxSetField(plhs[0],0,"Cs",field6);

    mxArray *field7 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field7) = ctfmodel.Ca;
    mxSetField(plhs[0],0,"Ca",field7);

    mxArray *field8 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field8) = ctfmodel.espr;
    mxSetField(plhs[0],0,"espr",field8);

    mxArray *field9 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field9) = ctfmodel.ispr;
    mxSetField(plhs[0],0,"ispr",field9);
    
    mxArray *field10 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field10) = ctfmodel.alpha;
    mxSetField(plhs[0],0,"alpha",field10);

    mxArray *field11 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field11) = ctfmodel.DeltaF;
    mxSetField(plhs[0],0,"DeltaF",field11);

    mxArray *field12 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field12) = ctfmodel.DeltaR;
    mxSetField(plhs[0],0,"DeltaR",field12);

    mxArray *field13 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field13) = ctfmodel.Q0;
    mxSetField(plhs[0],0,"Q0",field13);

    mxArray *field14 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field14) = ctfmodel.base_line;
    mxSetField(plhs[0],0,"base_line",field14);

    mxArray *field15 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field15) = ctfmodel.sqrt_K;
    mxSetField(plhs[0],0,"sqrt_K",field15);

    mxArray *field16 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field16) = ctfmodel.sqU;
    mxSetField(plhs[0],0,"sqU",field16);

    mxArray *field17 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field17) = ctfmodel.sqV;
    mxSetField(plhs[0],0,"sqV",field17);

    mxArray *field18 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field18) = ctfmodel.sqrt_angle;
    mxSetField(plhs[0],0,"sqrt_angle",field18);

    mxArray *field19 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field19) = ctfmodel.gaussian_K;
    mxSetField(plhs[0],0,"gaussian_K",field19);

    mxArray *field20 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field20) = ctfmodel.sigmaU;
    mxSetField(plhs[0],0,"sigmaU",field20);

    mxArray *field21 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field21) = ctfmodel.sigmaV;
    mxSetField(plhs[0],0,"sigmaV",field21);

    mxArray *field22 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field22) = ctfmodel.gaussian_angle;
    mxSetField(plhs[0],0,"gaussian_angle",field22);

    mxArray *field23 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field23) = ctfmodel.cU;
    mxSetField(plhs[0],0,"cU",field23);

    mxArray *field24 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field24) = ctfmodel.cV;
    mxSetField(plhs[0],0,"cV",field24);

    mxArray *field25 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field25) = ctfmodel.gaussian_K2;
    mxSetField(plhs[0],0,"gaussian_K2",field25);

    mxArray *field26 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field26) = ctfmodel.sigmaU2;
    mxSetField(plhs[0],0,"sigmaU2",field26);

    mxArray *field27 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field27) = ctfmodel.sigmaV2;
    mxSetField(plhs[0],0,"sigmaV2",field27);
    
    mxArray *field28 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field28) = ctfmodel.gaussian_angle2;
    mxSetField(plhs[0],0,"gaussian_angle2",field28);

    mxArray *field29 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field29) = ctfmodel.cU2;
    mxSetField(plhs[0],0,"cU2",field29);

    mxArray *field30 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field30) = ctfmodel.cV2;
    mxSetField(plhs[0],0,"cV2",field30);
    
    
    // Get the ctfmodels
    if (adjustParams.ctfmodelSize!=0)
    {
        MultidimArray<double> half;
        adjustParams.generate_model_halfplane(
        adjustParams.ctfmodelSize,adjustParams.ctfmodelSize,half);
        mxArray *field31;
        setMatrix2D(half,field31);
        mxSetField(plhs[0],0,"CTFmodelhalf",field31);
        
        MultidimArray<double> quadrant;
        adjustParams.generate_model_quadrant(
        adjustParams.ctfmodelSize,adjustParams.ctfmodelSize,quadrant);
        mxArray *field32;
        setMatrix2D(quadrant,field32);
        mxSetField(plhs[0],0,"CTFmodelquadrant",field32);
        
        MultidimArray<double> zeros(10,100,2);
        for (int n=0; n<ZSIZE(zeros); n++) {
            for (int iu=0; iu<YSIZE(zeros); iu++) {
                Matrix1D<double> u(2);
                VECTOR_R2(u,cos(iu*2*PI/YSIZE(zeros)),sin(iu*2*PI/YSIZE(zeros)));
                Matrix1D<double> contfreq(2), digfreq(2);
                ctfmodel.lookFor(n+1,u,contfreq);
                contfreq2digfreq(contfreq,digfreq,ctfmodel.Tm);
                digfreq*=adjustParams.ctfmodelSize;
                digfreq+=adjustParams.ctfmodelSize/2+1;
                zeros(n,iu,0)=XX(digfreq);
                zeros(n,iu,1)=YY(digfreq);
            }
        }
        mxArray *field33;
        setMatrix3D(zeros,field33);
        mxSetField(plhs[0],0,"zeros",field33);   
    }
    
    mxArray *field34 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field34) = adjustParams.Tm;
    mxSetField(plhs[0],0,"Tm",field34);
    
}	
