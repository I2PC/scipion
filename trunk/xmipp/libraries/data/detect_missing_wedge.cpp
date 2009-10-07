/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2009)
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

#include "detect_missing_wedge.h"
#include "args.h"
#include "fftw.h"
#include "de_solver.h"

// Evaluate plane ----------------------------------------------------------
double evaluatePlane(double rot, double tilt, 
    const Matrix3D<double> &V, const Matrix3D<double> &Vmag,
    double maxFreq, double planeWidth, int direction,
    Matrix3D<double> *Vdraw=NULL,
    bool setPos=false, double rotPos=0, double tiltPos=0)
{
    if (rot<0 || rot>360 || tilt<-90 || tilt>90)
        return 0;

    Matrix2D<double> E, Einv;
    Euler_angles2matrix(rot,tilt,0,E);
    Einv=E.transpose();
    
    if (setPos) {
        Matrix2D<double> Epos;
        Euler_angles2matrix(rotPos,tiltPos,0,Epos);
        double angle=acos(E(2,0)*Epos(2,0)+E(2,1)*Epos(2,1)+E(2,2)*Epos(2,2));
        angle=RAD2DEG(angle);
        if (ABS(angle)<20 || ABS(180-angle)<20) return 0;
    }
    
    int N=XMIPP_MAX(XSIZE(Vmag),YSIZE(Vmag)/2);
    N=XMIPP_MAX(N,ZSIZE(Vmag)/2);
    double df=0.5/N;
    Matrix1D<double> freq(3);
    Matrix1D<int> idx(3);
    double sumNeg=0, sumPos=0;
    int Nneg=0, Npos=0;
    for (double ix=0; ix<=N; ix++)
        for (double iy=-N; iy<=N; iy++)
        {
            VECTOR_R3(freq,ix*df,iy*df,0);
            if (freq.module()>maxFreq) continue;
            for (int iz=-planeWidth; iz<=planeWidth; iz++)
            {
                if (iz==0 || ix==0 || iy==0) continue;
                
                // Frequency in the coordinate system of the plane
                VECTOR_R3(freq,ix*df,iy*df,iz*df);

                // Frequency in the coordinate system of the volume
                freq=Einv*freq;
                bool inverted=false;
                if (XX(freq)<0)
                {
                    XX(freq)=-XX(freq);
                    YY(freq)=-YY(freq);
                    ZZ(freq)=-ZZ(freq);
                    inverted=true;
                }

                // Get the corresponding index
                digfreq2FFT_idx(V,freq,idx);
                if (Vmag.outside(ZZ(idx),YY(idx),XX(idx)))
                    continue;

                // Make the corresponding sums
                bool negativeSum;
                if (direction==1) negativeSum=iz<0;
                else              negativeSum=iz>0;
                if (negativeSum ^ inverted) // XOR
                {
                    sumNeg+=Vmag(idx);
                    Nneg++;
                    if (Vdraw!=NULL)
                        (*Vdraw)(idx)=2*direction*Vmag(idx);
                }
                else
                {
                    sumPos+=Vmag(idx);
                    Npos++;
                    if (Vdraw!=NULL)
                        (*Vdraw)(idx)=1.0/2.0*direction*Vmag(idx);
                }
            }
        }
    if (ABS(Nneg-Npos)/(0.5*(Nneg+Npos))>0.5)
        // If there is a difference of more than 50%
        return 1e38;
    if (Nneg!=0) sumNeg/=Nneg;
    else return 1e38;
    if (Npos!=0) sumPos/=Npos;
    else return 1e38;
    
    return -(sumPos-sumNeg);
}

// Look for plane ----------------------------------------------------------
class WedgeSolver: public DESolver {
public:
    WedgeSolver(int dim,int pop,
        const Matrix3D<double> &_V, const Matrix3D<double> &_Vmag,
        double _maxFreq, double _planeWidth, int _direction):
        DESolver(dim,pop) {
            count=0;
            V=&_V;
            Vmag=&_Vmag;
            maxFreq=_maxFreq;
            planeWidth=_planeWidth;
            direction=_direction;
            setPos=false;
        }
    double EnergyFunction(double trial[],bool &bAtSolution) {
        double result=evaluatePlane(trial[0],trial[1],*V,*Vmag,
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
   const Matrix3D<double> *V;
   const Matrix3D<double> *Vmag;
   double maxFreq;
   double planeWidth;
   int direction;
   bool setPos;
   double rotPos;
   double tiltPos;
};

double wrapperFitness(double *p, void* extraArgs)
{
    bool dummy;
    WedgeSolver *wegde_solver=(WedgeSolver *) extraArgs;
    return wegde_solver->EnergyFunction(p+1,dummy);
}

void lookForPlane(const Matrix3D<double> &V, const Matrix3D<double> &Vmag,
    double maxFreq, double planeWidth, int direction,
    double &rot, double &tilt,
    bool setPos=false, double rotPos=0, double tiltPos=0)
{
    // Optimize with DE
    int length=2;
    int Npop=50;
    WedgeSolver solver(length,length*Npop,V,Vmag,maxFreq,planeWidth,direction);
    Matrix1D<double> min_allowed(2), max_allowed(2);
    min_allowed(0)=0;   min_allowed(1)=-90;
    max_allowed(0)=360; max_allowed(1)= 90;
    solver.Setup(MULTIDIM_ARRAY(min_allowed),
                 MULTIDIM_ARRAY(max_allowed), stBest2Bin, 0.5, 0.8);
    if (setPos)
    {
        solver.setPos=true;
        solver.rotPos=rotPos;
        solver.tiltPos=tiltPos;
    }
    solver.Solve(50);
    double* bestSolution=solver.Solution();

    // Refine with Powell
    Matrix1D<double> x(2), steps(2);
    double fitness;
    int iter;
    steps.initConstant(1);
    x(0)=bestSolution[0];
    x(1)=bestSolution[1];
    powellOptimizer(x,1,2,&wrapperFitness,&solver,0.01,fitness,
        iter,steps,true);
    rot=x(0);
    tilt=x(1);
}

// Draw wedge --------------------------------------------------------------
void drawWedge(double rotPos, double tiltPos, double rotNeg, double tiltNeg,
    const Matrix3D<double> &V,
    const Matrix3D<double> &Vmag, Matrix3D<double> &Vdraw)
{
    Matrix2D<double> Epos, Eneg;
    Euler_angles2matrix(rotPos,tiltPos,0,Epos);
    Euler_angles2matrix(rotNeg,tiltNeg,0,Eneg);
    
    Matrix1D<double> freq(3), freqPos, freqNeg;
    Matrix1D<int> idx(3);
    Vdraw.initZeros(Vmag);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vdraw)
    {
        // Frequency in the coordinate system of the volume
        VECTOR_R3(idx,j,i,k);
        FFT_idx2digfreq(V,idx,freq);

        // Frequency in the coordinate system of the plane
        freqPos=Epos*freq;
        freqNeg=Eneg*freq;
        if (ZZ(freqPos)<0 || ZZ(freqNeg)>0)
            Vdraw(k,i,j)+=1;
    }
}

// Remove wedge ------------------------------------------------------------
void MissingWedge::removeWedge(Matrix3D<double> &V) const
{
    Matrix2D<double> Epos, Eneg;
    Euler_angles2matrix(rotPos,tiltPos,0,Epos);
    Euler_angles2matrix(rotNeg,tiltNeg,0,Eneg);
    
    Matrix1D<double> freq(3), freqPos, freqNeg;
    Matrix1D<int> idx(3);

    XmippFftw transformer;
    Matrix3D< std::complex<double> > Vfft;
    transformer.FourierTransform(V,Vfft,false);

    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vfft)
    {
        // Frequency in the coordinate system of the volume
        VECTOR_R3(idx,j,i,k);
        FFT_idx2digfreq(V,idx,freq);

        // Frequency in the coordinate system of the plane
        freqPos=Epos*freq;
        freqNeg=Eneg*freq;
        if (ZZ(freqPos)<0 || ZZ(freqNeg)>0)
            Vfft(k,i,j)=0;
    }
    transformer.inverseFourierTransform();
}

// Read from command line --------------------------------------------------
void DetectMissingWedge_parameters::read(int argc, char **argv)
{
    fn_vol  = getParameter(argc, argv, "-i");
    maxFreq = textToFloat(getParameter(argc, argv, "-maxFreq", "2"));
    planeWidth = textToFloat(getParameter(argc, argv, "-width", "2"));
    saveMarks = checkParameter(argc, argv, "-saveMarks");
    saveMask = checkParameter(argc, argv, "-saveMask");
}

// Produce side info -------------------------------------------------------
void DetectMissingWedge_parameters::produceSideInfo()
{
    V.read(fn_vol);
    XmippFftw transformer;
    Matrix3D< std::complex<double> > Vfft;
    transformer.FourierTransform(V(),Vfft,false);
    Vmag.resize(Vfft);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vmag)
        Vmag(k,i,j)=20*log10(abs(Vfft(k,i,j)));
}

// Show --------------------------------------------------------------------
void DetectMissingWedge_parameters::show() const
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
void DetectMissingWedge_parameters::usage() const
{
    std::cerr << "detect_missing_wedge\n"
              << "  -i <volume>        : Input tomogram\n"
              << " [-maxFreq <f=0.25>] : Normalized to 0.5\n"
              << " [-width <w=2>]      : Width of the probe plane\n"
              << " [-saveMarks]        : Save the magnitude of the FFT with\n"
              << "                       marks showing the two planes\n"
              << " [-saveMask]         : Save a mask for the FFT of this tomogram\n"
              << "                       1=Missing wedge, 0=non missing wedge\n"
    ;
}

// Run ---------------------------------------------------------------------
void DetectMissingWedge_parameters::run()
{
    VolumeXmipp Vdraw;
    FileName fn_root=V.name().without_extension();

    // Detect one of the planes
    lookForPlane(V(), Vmag, maxFreq, planeWidth, 1, rotPos, tiltPos);
    if (saveMarks)
    {
        Vdraw()=Vmag;
        evaluatePlane(rotPos, tiltPos, V(), Vmag, maxFreq, planeWidth,
            1, &(Vdraw()));
    }

    // Detect the other plane
    lookForPlane(V(), Vmag, maxFreq, planeWidth, -1, rotNeg, tiltNeg,
        true, rotPos, tiltPos);
    if (saveMarks)
    {
        evaluatePlane(rotNeg, tiltNeg, V(), Vmag, maxFreq, planeWidth,
            -1, &(Vdraw()));
        Vdraw.write(fn_root+"_marks.vol");
    }

    if (saveMask)
    {
        Vdraw.clear();
        drawWedge(rotPos, tiltPos, rotNeg, tiltNeg, V(), Vmag, Vdraw());
        Vdraw.write(fn_root+"_mask.vol");
    }

    std::cout << "Plane1: " << rotPos << " " << tiltPos << std::endl;
    std::cout << "Plane2: " << rotNeg << " " << tiltNeg << std::endl;
}
