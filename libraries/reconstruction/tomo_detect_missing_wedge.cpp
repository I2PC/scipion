/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2009)
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

#include "tomo_detect_missing_wedge.h"
#include <data/args.h>
#include <data/xmipp_fftw.h>
#include <data/numerical_tools.h>

// Evaluate plane ----------------------------------------------------------
double evaluatePlane(double rot, double tilt,
                     const MultidimArray<double> *V, const MultidimArray<double> *Vmag,
                     double maxFreq, double planeWidth, int direction,
                     MultidimArray<double> *Vdraw=NULL,
                     bool setPos=false, double rotPos=0, double tiltPos=0)
{
    if (rot<0 || rot>360 || tilt<-90 || tilt>90)
        return 0;

    Matrix2D<double> E, Einv;
    Euler_angles2matrix(rot,tilt,0,E);
    Einv=E.transpose();

    if (setPos)
    {
        Matrix2D<double> Epos;
        Euler_angles2matrix(rotPos,tiltPos,0,Epos);
        double angle=acos(E(2,0)*Epos(2,0)+E(2,1)*Epos(2,1)+E(2,2)*Epos(2,2));
        angle=RAD2DEG(angle);
        if (ABS(angle)<20 || ABS(180-angle)<20)
            return 0;
    }

    size_t N=XMIPP_MAX(XSIZE(*Vmag),YSIZE(*Vmag)/2);
    N=XMIPP_MAX(N,ZSIZE(*Vmag)/2);
    double df=0.5/N;
    Matrix1D<double> freq(3), freqp(3);
    Matrix1D<int> idx(3);
    double sumNeg=0, sumPos=0;
    int Nneg=0, Npos=0;
    double maxFreq2=maxFreq*maxFreq;
    int iPlaneWidth=(int)ceil(planeWidth);
    for (double ix=0; ix<=N; ix++)
    {
        XX(freq)=ix*df;
        double fx2=XX(freq)*XX(freq);
        for (double iy=-N; iy<=N; iy++)
        {
            YY(freq)=iy*df;
            double fx2fy2=fx2+YY(freq)*YY(freq);
            if (fx2fy2>maxFreq2)
                continue;
            for (int iz=-iPlaneWidth; iz<=iPlaneWidth; iz++)
            {
                if (iz==0 || ix==0 || iy==0)
                    continue;

                // Frequency in the coordinate system of the plane
                ZZ(freq)=iz*df;

                // Frequency in the coordinate system of the volume
                SPEED_UP_temps012;
                M3x3_BY_V3x1(freqp,Einv,freq);
                bool inverted=false;
                if (XX(freqp)<0)
                {
                    XX(freqp)=-XX(freqp);
                    YY(freqp)=-YY(freqp);
                    ZZ(freqp)=-ZZ(freqp);
                    inverted=true;
                }

                // Get the corresponding index
                DIGFREQ2FFT_IDX(ZZ(freqp), ZSIZE(*V), ZZ(idx));
                DIGFREQ2FFT_IDX(YY(freqp), YSIZE(*V), YY(idx));
                DIGFREQ2FFT_IDX(XX(freqp), XSIZE(*V), XX(idx));
                if (XX(idx) < STARTINGX(*Vmag) || XX(idx) > FINISHINGX(*Vmag) ||
                    YY(idx) < STARTINGY(*Vmag) || YY(idx) > FINISHINGY(*Vmag) ||
                    ZZ(idx) < STARTINGZ(*Vmag) || ZZ(idx) > FINISHINGZ(*Vmag))
                    continue;

                // Make the corresponding sums
                bool negativeSum;
                if (direction==1)
                    negativeSum=iz<0;
                else
                    negativeSum=iz>0;
                double val=A3D_ELEM(*Vmag,ZZ(idx),YY(idx),XX(idx));
                if (negativeSum ^ inverted) // XOR
                {
                    sumNeg+=val;
                    Nneg++;
                    if (Vdraw!=NULL)
                        (*Vdraw)(idx)=2*direction*val;
                }
                else
                {
                    sumPos+=val;
                    Npos++;
                    if (Vdraw!=NULL)
                        (*Vdraw)(idx)=1.0/2.0*direction*val;
                }
            }
        }
    }
    if (fabs(Nneg-Npos)/(0.5*(Nneg+Npos))>0.5)
        // If there is a difference of more than 50%
        return 1e38;
    if (Nneg!=0)
        sumNeg/=Nneg;
    else
        return 1e38;
    if (Npos!=0)
        sumPos/=Npos;
    else
        return 1e38;

    return -(sumPos-sumNeg);
}

// Look for plane ----------------------------------------------------------
class WedgeSolver: public DESolver
{
public:
    WedgeSolver(int dim,int pop,
                const MultidimArray<double> *_V, const MultidimArray<double> *_Vmag,
                double _maxFreq, double _planeWidth, int _direction):
            DESolver(dim,pop)
    {
        count=0;
        V=_V;
        Vmag=_Vmag;
        maxFreq=_maxFreq;
        planeWidth=_planeWidth;
        direction=_direction;
        setPos=false;
    }

    ~WedgeSolver()
    {
        V = NULL;
        Vmag = NULL;
    }

    double EnergyFunction(double trial[],bool &bAtSolution)
    {
        double result=evaluatePlane(trial[0],trial[1],V,Vmag,
                                    maxFreq, planeWidth, direction, NULL, setPos, rotPos, tiltPos);
        if (count++ % (5*nPop) == 0)
            std::cout << "Evaluations= " << count/nPop
            << " energy= "     << Energy()
            << " rot= "        << Solution()[0]
            << " tilt= "       << Solution()[1]
            << std::endl;
        bAtSolution=false;
        return result;
    }
public:
    int count;
    const MultidimArray<double> *V;
    const MultidimArray<double> *Vmag;
    double maxFreq;
    double planeWidth;
    int direction;
    bool setPos;
    double rotPos;
    double tiltPos;
};

//// Sjors wrapperFitness also exists in nma_alignment
double wrapperFitnessDetectMissingWedge(double *p, void* extraArgs)
{
    bool dummy;
    WedgeSolver *wegde_solver=(WedgeSolver *) extraArgs;
    return wegde_solver->EnergyFunction(p+1,dummy);
}

void lookForPlane(const MultidimArray<double> *V, const MultidimArray<double> *Vmag,
                  double maxFreq, double planeWidth, int direction,
                  double &rot, double &tilt,
                  bool setPos=false, double rotPos=0, double tiltPos=0)
{
    // Optimize with DE
    int length=2;
    int Npop=50;

    double min_allowed[] = {0,-90};
    double max_allowed[] = {360,90};

    WedgeSolver * solver = new WedgeSolver(length,length*Npop,V,Vmag,maxFreq,planeWidth,direction);
    solver->Setup( min_allowed, max_allowed, stBest2Bin, 0.5, 0.8);

    if (setPos)
    {
        solver->setPos=true;
        solver->rotPos=rotPos;
        solver->tiltPos=tiltPos;
    }

    solver->Solve(50);

    double * bestSolution = solver->Solution();

    // Refine with Powell
    Matrix1D<double> x(2), steps(2);
    double fitness;
    int iter;
    steps.initConstant(1);
    x(0)=bestSolution[0];
    x(1)=bestSolution[1];
    powellOptimizer(x,1,2,&wrapperFitnessDetectMissingWedge,solver,0.01,fitness,iter,steps,true);
    rot=x(0);
    tilt=x(1);

    delete solver;

    solver = NULL;
}

// Draw wedge --------------------------------------------------------------
void drawWedge(double rotPos, double tiltPos, double rotNeg, double tiltNeg,
               const MultidimArray<double> *V,
               const MultidimArray<double> *Vmag, MultidimArray<double> *Vdraw)
{
    Matrix2D<double> Epos, Eneg;
    Euler_angles2matrix(rotPos,tiltPos,0,Epos);
    Euler_angles2matrix(rotNeg,tiltNeg,0,Eneg);

    Matrix1D<double> freq(3), freqPos, freqNeg;
    Matrix1D<int> idx(3);
    Vdraw->initZeros(*Vmag);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(*Vdraw)
    {
        // Frequency in the coordinate system of the volume
        VECTOR_R3(idx,j,i,k);
        FFT_idx2digfreq(*V,idx,freq);

        // Frequency in the coordinate system of the plane
        freqPos=Epos*freq;
        freqNeg=Eneg*freq;
        if (ZZ(freqPos)<0 || ZZ(freqNeg)>0)
            (*Vdraw)(k,i,j)+=1;
    }
}

// Read from command line --------------------------------------------------
void ProgDetectMissingWedge::readParams()
{
    fn_vol  = getParam("-i");
    maxFreq = getDoubleParam("--maxFreq");
    planeWidth = getDoubleParam("--width");
    saveMarks = checkParam("--saveMarks");
    saveMask = checkParam("--saveMask");
}

// Produce side info -------------------------------------------------------
void ProgDetectMissingWedge::produceSideInfo()
{
    V.read(fn_vol);
    FourierTransformer transformer;
    MultidimArray< std::complex<double> > Vfft;
    transformer.FourierTransform(MULTIDIM_ARRAY(V),Vfft,false);
    Vmag = new MultidimArray<double>();
    Vmag->resize(Vfft);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(*Vmag)
    (*Vmag)(k,i,j)=20*log10(abs(Vfft(k,i,j)));
}

// Show --------------------------------------------------------------------
void ProgDetectMissingWedge::show() const
{
    std::cout << "Detecting a missing wedge\n";
    std::cout << "Input volume:   " << fn_vol     << std::endl
    << "Maximum Freq.:  " << maxFreq    << std::endl
    << "Plane width:    " << planeWidth << std::endl
    << "saveMarks:      " << saveMarks  << std::endl
    << "saveMask:       " << saveMask   << std::endl
    ;
}

// Usage -------------------------------------------------------------------
void ProgDetectMissingWedge::defineParams()
{
    addUsageLine("Detect the orientation of the missing wedge in a tomogram. ");
    addUsageLine("+For doing so it fits a couple of planes along which there is a maximum ");
    addUsageLine("+variation between the energy of the Fourier transform on its left and ");
    addUsageLine("+on its right. The missing wedge is coded with four angles (two for each ");
    addUsageLine("+plane). You may also produce a mask with 1 where the missing wedge is, ");
    addUsageLine("+and 0 where the data has been actually measured. Finally, you can also ");
    addUsageLine("+produce a marked magnitude volume (i.e., the magnitude of the Fourier ");
    addUsageLine("+transform where the position of the two planes have been marked). ");
    addParamsLine("  -i <file>           : Input tomogram");
    addParamsLine(" [--maxFreq <f=0.25>] : Maximum frequency to fit the plane (normalized to 0.5)");
    addParamsLine(" [--width <w=2>]      : Width of the probe plane");
    addParamsLine(" [--saveMarks]        : Save the magnitude of the FFT with");
    addParamsLine("                      : marks showing the two planes");
    addParamsLine(" [--saveMask]         : Save a mask for the FFT of this tomogram");
    addParamsLine("                      : 1=Missing wedge, 0=non missing wedge");
    addExampleLine("xmipp_tomo_detect_missing_wedge -i tomogram.vol");
}

// Run ---------------------------------------------------------------------
void ProgDetectMissingWedge::run()
{
    produceSideInfo();
    FileName fn_root=V.name().withoutExtension();

    MultidimArray<double> * mdaV = &MULTIDIM_ARRAY(V);

    // Detect one of the planes
    lookForPlane(mdaV, Vmag, maxFreq, planeWidth, 1, rotPos, tiltPos);

    Image<double> * Vdraw = new Image<double>();

    if (saveMarks)
    {
        (*Vdraw)()=(*Vmag);

        evaluatePlane(rotPos, tiltPos, mdaV, Vmag, maxFreq, planeWidth,
                      1, &(*Vdraw)());
    }

    // Detect the other plane
    lookForPlane(mdaV, Vmag, maxFreq, planeWidth, -1, rotNeg, tiltNeg, true, rotPos, tiltPos);

    if (saveMarks)
    {
        evaluatePlane(rotNeg, tiltNeg, mdaV, Vmag, maxFreq, planeWidth,
                      -1, &(*Vdraw)());
        Vdraw->write(fn_root+"_marks.vol");
    }

    if (saveMask)
    {
        Vdraw->clear();
        drawWedge(rotPos, tiltPos, rotNeg, tiltNeg, mdaV, Vmag, &(*Vdraw)());
        Vdraw->write(fn_root+"_mask.vol");
    }

    std::cout << "Plane1: " << rotPos << " " << tiltPos << std::endl;
    std::cout << "Plane2: " << rotNeg << " " << tiltNeg << std::endl;

    delete Vdraw;
    delete Vmag;
}
