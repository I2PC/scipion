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

#include "convert_vol2pseudo.h"
#include "fourier_filter.h"
#include <algorithm>
#include <stdio.h>

/* Pseudo atoms ------------------------------------------------------------ */
PseudoAtom::PseudoAtom()
{
    location.initZeros(3);
    intensity=0;
}

bool operator <(const PseudoAtom &a, const PseudoAtom &b)
{
    return a.intensity<b.intensity;
}

std::ostream& operator << (std::ostream &o, const PseudoAtom &a)
{
    o << a.location.transpose() << " " << a.intensity;
    return o;
}

/* I/O --------------------------------------------------------------------- */
void Prog_Convert_Vol2Pseudo::read(int argc, char **argv)
{
    fnVol = getParameter(argc,argv,"-i");
    fnOut = getParameter(argc,argv,"-o","");
    useMask = checkParameter(argc, argv, "-mask");
    if (useMask) mask_prm.read(argc, argv);
    sigma = textToFloat(getParameter(argc,argv,"-sigma","1.5"));
    targetError = textToFloat(getParameter(argc,argv,"-targetError","0.02"));
    stop = textToFloat(getParameter(argc,argv,"-stop","0.001"));
    initialSeeds = textToInteger(getParameter(argc,argv,"-initialSeeds","300"));
    growSeeds = textToFloat(getParameter(argc,argv,"-growSeeds","30"));
    allowMovement = !checkParameter(argc,argv,"-dontAllowMovement");
    allowIntensity = !checkParameter(argc,argv,"-dontAllowIntensity");
    intensityColumn = getParameter(argc,argv,"-intensityColumn","occupancy");
}

void Prog_Convert_Vol2Pseudo::show() const
{
    std::cout << "Input volume:   " << fnVol           << std::endl
              << "Output volume:  " << fnOut           << std::endl
              << "Sigma:          " << sigma           << std::endl
              << "Initial seeds:  " << initialSeeds    << std::endl
              << "Grow seeds:     " << growSeeds       << std::endl
              << "Target error:   " << targetError     << std::endl
              << "Stop:           " << stop            << std::endl
              << "AllowMovement:  " << allowMovement   << std::endl
              << "AllowIntensity: " << allowIntensity  << std::endl
              << "Intensity Col:  " << intensityColumn << std::endl
    ;    
    if (useMask) mask_prm.show();
    else std::cout << "No mask\n";
}

void Prog_Convert_Vol2Pseudo::usage() const
{
    std::cout << "Approximation algorithm:\n"
              << "   -i <volume>                     : Input volume\n"
              << "  [-o <rootname>]                  : Output rootname\n"
              << "  [-sigma <s=1.5>]                 : Sigma of gaussians\n"
              << "  [-initialSeeds <N=300>]          : Initial number of gaussians\n"
              << "  [-growSeeds <%=30>]              : Percentage of growth\n"
              << "  [-stop <p=0.001>]                : Stop criterion (0<p<1) for inner iterations\n"
              << "  [-targetError <e=0.02>]          : Finish when the average representation\n"
              << "                                     error is below this threshold\n"
              << "  [-dontAllowMovement]             : Don't allow Gaussians to move\n"
              << "  [-dontAllowIntensity]            : Don't allow Gaussians to change intensity\n"
              << "  [-intensityColumn <s=occupancy>] : Where to write the intensity in the PDB file\n"
              << "                                     Valid values: occupancy, Bfactor\n"
    ;
    mask_prm.usage();
}

void Prog_Convert_Vol2Pseudo::produceSideInfo()
{
    if (intensityColumn!="occupancy" && intensityColumn!="Bfactor")
        REPORT_ERROR(1,(std::string)"Unknown column: "+intensityColumn);

    Vin.read(fnVol);
    Vin().setXmippOrigin();
    
    Vcurrent().initZeros(Vin());
    mask_prm.generate_3Dmask(Vin());
    
    sigma3=3*sigma;
    gaussianTable.resize(CEIL(sigma3*sqrt(3)*1000));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(gaussianTable)
        gaussianTable(i)=gaussian1D(i/1000.0,sigma);
    
    energyOriginal=0;
    double N=0;
    double minval=1e38, maxval=-1e38;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vin())
    {
        if (useMask && mask_prm.imask3D(k,i,j)==0) continue;
        double v=Vin(k,i,j);
        energyOriginal+=v*v;
        minval=XMIPP_MIN(minval,v);
        maxval=XMIPP_MAX(maxval,v);
        N++;
    }
    energyOriginal/=N;
    
    histogram1D hist;
    if (useMask)
        compute_hist_within_binary_mask(mask_prm.imask3D, Vin(), hist,
            minval, maxval, 200);
    else
        compute_hist(Vin(), hist, minval, maxval, 200);
    percentil1=hist.percentil(1);
    if (percentil1<=0)
        percentil1=maxval/500;
    range=hist.percentil(99)-percentil1;
    smallAtom=(range*targetError)/2;
}

//#define DEBUG
void Prog_Convert_Vol2Pseudo::placeSeeds(int Nseeds)
{
    // Convolve the difference with the Gaussian to know
    // where it would be better to put a Gaussian
    FourierMask Filter;
    Filter.FilterBand=LOWPASS;
    Filter.FilterShape=REALGAUSSIAN;
    Filter.w1=sigma;
    Filter.generate_mask(Vin());
    Filter.do_generate_3dmask=false;
 
    Matrix3D<double> Vdiff=Vin();
    Vdiff-=Vcurrent();
    Filter.apply_mask_Space(Vdiff);
    
    // Place all seeds
    int rmax=3*sigma;
    for (int n=0; n<Nseeds; n++)
    {
        // Look for the maximum error
        bool first=true;
        int kmax, imax, jmax;
        double maxVal;
        FOR_ALL_ELEMENTS_IN_MATRIX3D(Vdiff)
        {
            if (useMask && mask_prm.imask3D(k,i,j)==0) continue;
            if (first || Vdiff(k,i,j)>maxVal)
            {
                kmax=k;
                imax=i;
                jmax=j;
                maxVal=Vdiff(k,i,j);
                first=false;
            }
        }
        
        // Keep this as an atom
        PseudoAtom a;
        a.location(0)=kmax;
        a.location(1)=imax;
        a.location(2)=jmax;
        if (allowIntensity) a.intensity=maxVal;
        else 
        {
            if (maxVal<smallAtom) break;
            a.intensity=smallAtom;
        }
        atoms.push_back(a);
        
        // Remove this density from the difference
        drawGaussian(kmax,imax,jmax,Vdiff,-a.intensity);
        
        #ifdef DEBUG
            std::cout << "New atom: " << a << std::endl;
            VolumeXmipp save;
            save()=Vdiff; save.write("PPPDiff.vol");
            std::cout << "Press any key\n";
            char c; std::cin >> c;
        #endif
    }
} 
#undef DEBUG

/* Remove seeds ------------------------------------------------------------ */
void Prog_Convert_Vol2Pseudo::removeSeeds(int Nseeds)
{
    int fromNegative=ROUND(Nseeds*0.5);
    int fromSmall=Nseeds-fromNegative;

    if (allowIntensity)
    {
        // Remove too small atoms
        std::sort(atoms.begin(),atoms.end());
        atoms.erase(atoms.begin(),atoms.begin()+fromSmall);
    }
    else
    {
        fromNegative=Nseeds;
        fromSmall=0;
    }        
    
    // Remove atoms from regions in which the error is too negative
    Matrix3D<double> Vdiff=Vin();
    Vdiff-=Vcurrent();
    int alreadyRemoved=0;
    double vmin=Vdiff.computeMin();
    if (vmin<0)
    {
        for (double v=vmin+vmin/20; v<0; v-=vmin/20)
        {
            int oldListSize;
            do {
                oldListSize=atoms.size();

                // Search for a point within a negative region
                bool found=false;
                int kneg, ineg, jneg;
                for (int k=STARTINGZ(Vdiff); k<=FINISHINGZ(Vdiff) && !found; k++)
                    for (int i=STARTINGY(Vdiff); i<=FINISHINGY(Vdiff) && !found; i++)
                        for (int j=STARTINGX(Vdiff); j<=FINISHINGX(Vdiff) && !found; j++)
                        {
                            if (useMask && mask_prm.imask3D(k,i,j)==0) continue;
                            if (Vdiff(k,i,j)<v)
                            {
                                kneg=k;
                                ineg=i;
                                jneg=j;
                                Vdiff(k,i,j)=0;
                                found=true;
                            }
                        }

                // If found such a point, search for a nearby atom
                if (found)
                {
                    // Search for atom
                    int nmax=atoms.size();
                    for (int n=0; n<nmax; n++)
                    {
                        double r=
                            (kneg-atoms[n].location(0))*(kneg-atoms[n].location(0))+
                            (ineg-atoms[n].location(1))*(ineg-atoms[n].location(1))+
                            (jneg-atoms[n].location(2))*(jneg-atoms[n].location(2));                        r=sqrt(r);
                        if (r<sigma3)
                        {
                            drawGaussian(atoms[n].location(0),
                                atoms[n].location(1),
                                atoms[n].location(2),
                                Vdiff,
                                atoms[n].intensity);
                            atoms.erase(atoms.begin()+n);
                            alreadyRemoved++;
                            break;
                        }
                    }
                }
            } while (oldListSize>atoms.size() && alreadyRemoved<fromNegative);
            if (alreadyRemoved==fromNegative) break;
        }
    }
}

/* Draw approximation ------------------------------------------------------ */
void Prog_Convert_Vol2Pseudo::drawApproximation()
{
    Vcurrent().initZeros(Vin());
    int nmax=atoms.size();
    for (int n=0; n<nmax; n++)
        drawGaussian(atoms[n].location(0),atoms[n].location(1),
            atoms[n].location(2),Vcurrent(),atoms[n].intensity);

    energyDiff=0;
    double N=0;
    percentageDiff=0;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vcurrent())
    {
        if (useMask && mask_prm.imask3D(k,i,j)==0) continue;
        double Vinv=Vin(k,i,j);
        if (Vinv<=0) continue;
        double vdiff=Vinv-Vcurrent(k,i,j);
        double vperc=ABS(vdiff);
        energyDiff+=vdiff*vdiff;
        percentageDiff+=vperc;
        N++;
    }
    energyDiff/=N;
    percentageDiff/=(N*range);
}

/* Gaussian operations ----------------------------------------------------- */
double Prog_Convert_Vol2Pseudo::computeAverage(int k, int i, int j,
    Matrix3D<double> &V)
{
    int k0=XMIPP_MAX(STARTINGZ(V),k-sigma3);
    int i0=XMIPP_MAX(STARTINGY(V),i-sigma3);
    int j0=XMIPP_MAX(STARTINGX(V),j-sigma3);
    int kF=XMIPP_MIN(FINISHINGZ(V),k+sigma3);
    int iF=XMIPP_MIN(FINISHINGY(V),i+sigma3);
    int jF=XMIPP_MIN(FINISHINGX(V),j+sigma3);
    double sum=0;
    for (int kk=k0; kk<=kF; kk++)
        for (int ii=i0; ii<=iF; ii++)
            for (int jj=j0; jj<=jF; jj++)
                sum+=V(kk,ii,jj);
    return sum/((kF-k0+1)*(iF-i0+1)*(jF-j0+1));
}

void Prog_Convert_Vol2Pseudo::drawGaussian(double k, double i, double j,
    Matrix3D<double> &V, double intensity)
{
    int k0=CEIL(XMIPP_MAX(STARTINGZ(V),k-sigma3));
    int i0=CEIL(XMIPP_MAX(STARTINGY(V),i-sigma3));
    int j0=CEIL(XMIPP_MAX(STARTINGX(V),j-sigma3));
    int kF=FLOOR(XMIPP_MIN(FINISHINGZ(V),k+sigma3));
    int iF=FLOOR(XMIPP_MIN(FINISHINGY(V),i+sigma3));
    int jF=FLOOR(XMIPP_MIN(FINISHINGX(V),j+sigma3));
    for (int kk=k0; kk<=kF; kk++)
    {
        double diffkk2=(kk-k)*(kk-k);
        for (int ii=i0; ii<=iF; ii++)
        {
            double diffiikk2=(ii-i)*(ii-i)+diffkk2;
            for (int jj=j0; jj<=jF; jj++)
            {
                double r=sqrt(diffiikk2+(jj-j)*(jj-j));
                V(kk,ii,jj)+=intensity*
                    DIRECT_VEC_ELEM(gaussianTable,ROUND(r*1000));
            }
        }
    }
}

void Prog_Convert_Vol2Pseudo::extractRegion(int idxGaussian,
    Matrix3D<double> &region, bool extended) const
{
    double k=atoms[idxGaussian].location(0);
    double i=atoms[idxGaussian].location(1);
    double j=atoms[idxGaussian].location(2);

    double sigma3ToUse=sigma3;
    if (extended)
        sigma3ToUse+=1.5;

    int k0=CEIL(XMIPP_MAX(STARTINGZ(Vcurrent()),k-sigma3ToUse));
    int i0=CEIL(XMIPP_MAX(STARTINGY(Vcurrent()),i-sigma3ToUse));
    int j0=CEIL(XMIPP_MAX(STARTINGX(Vcurrent()),j-sigma3ToUse));
    int kF=FLOOR(XMIPP_MIN(FINISHINGZ(Vcurrent()),k+sigma3ToUse));
    int iF=FLOOR(XMIPP_MIN(FINISHINGY(Vcurrent()),i+sigma3ToUse));
    int jF=FLOOR(XMIPP_MIN(FINISHINGX(Vcurrent()),j+sigma3ToUse));
    
    region.resize(kF-k0+1,iF-i0+1,jF-j0+1);
    STARTINGZ(region)=k0;
    STARTINGY(region)=i0;
    STARTINGX(region)=j0;
    for (int k=k0; k<=kF; k++)
        for (int i=i0; i<=iF; i++)
            for (int j=j0; j<=jF; j++)
                region(k,i,j)=Vcurrent(k,i,j);
}

double Prog_Convert_Vol2Pseudo::evaluateRegion(const Matrix3D<double> &region)
    const 
{
    double avgDiff=0;
    double N=0;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(region)
    {
        double Vinv=Vin(k,i,j);
        if (Vinv<=0) continue;
        if (useMask && mask_prm.imask3D(k,i,j)==0) continue;
        double vdiff=region(k,i,j)-Vinv;
        double vperc=ABS(vdiff);
        avgDiff+=vperc;
        N++;
    }
    return avgDiff/(N*range);
}

void Prog_Convert_Vol2Pseudo::insertRegion(const Matrix3D<double> &region)
{
    FOR_ALL_ELEMENTS_IN_MATRIX3D(region)
        Vcurrent(k,i,j)=VOL_ELEM(region,k,i,j);
}

/* Optimize ---------------------------------------------------------------- */
void Prog_Convert_Vol2Pseudo::optimizeCurrentAtoms()
{
    if (!allowIntensity && !allowMovement) return;
    bool finished=false;
    int iter=0;
    Matrix3D<double> region, regionBackup;
    do
    {
        double oldError=percentageDiff;
        int Nintensity=0;
        int Nmovement=0;
        int nmax=atoms.size();
        for (int n=0; n<nmax; n++)
        {
            extractRegion(n,region,true);
            double currentRegionEval=evaluateRegion(region);
            drawGaussian(atoms[n].location(0), atoms[n].location(1),
                atoms[n].location(2),region,-atoms[n].intensity);
            regionBackup=region;

            // Change intensity
            if (allowIntensity)
            {
                // Try with a Gaussian that is of different intensity
                double tryCoeffs[5]={0, 0.9, 0.99, 1.01, 1.1};
                double bestRed=0;
                int bestT=-1;
                for (int t=0; t<5; t++)
                {
                    region=regionBackup;
                    drawGaussian(atoms[n].location(0), atoms[n].location(1),
                        atoms[n].location(2),region,
                        tryCoeffs[t]*atoms[n].intensity);
                    double trialRegionEval=evaluateRegion(region);
                    double reduction=trialRegionEval-currentRegionEval;
                    if (reduction<bestRed)
                    {
                        bestRed=reduction;
                        bestT=t;
                    }
                }
                if (bestT!=-1)
                {
                    atoms[n].intensity*=tryCoeffs[bestT];
                    region=regionBackup;
                    drawGaussian(atoms[n].location(0), atoms[n].location(1),
                        atoms[n].location(2),region,atoms[n].intensity);
                    insertRegion(region);
                    currentRegionEval=evaluateRegion(region);
                    drawGaussian(atoms[n].location(0), atoms[n].location(1),
                        atoms[n].location(2),region,-atoms[n].intensity);
                    regionBackup=region;
                    Nintensity++;
                }
            }
            
            // Change location
            if (allowMovement && atoms[n].intensity>0)
            {
                double tryX[6]={-0.45,0.5, 0.0 ,0.0, 0.0 ,0.0};
                double tryY[6]={ 0.0 ,0.0,-0.45,0.5, 0.0 ,0.0};
                double tryZ[6]={ 0.0 ,0.0, 0.0 ,0.0,-0.45,0.5};
                double bestRed=0;
                int bestT=-1;
                for (int t=0; t<6; t++)
                {
                    region=regionBackup;
                    drawGaussian(atoms[n].location(0)+tryZ[t],
                        atoms[n].location(1)+tryY[t],
                        atoms[n].location(2)+tryX[t],
                        region,atoms[n].intensity);
                    double trialRegionEval=evaluateRegion(region);
                    double reduction=trialRegionEval-currentRegionEval;
                    if (reduction<bestRed)
                    {
                        bestRed=reduction;
                        bestT=t;
                    }
                }
                if (bestT!=-1)
                {
                    atoms[n].location(0)+=tryZ[bestT];
                    atoms[n].location(1)+=tryY[bestT];
                    atoms[n].location(2)+=tryX[bestT];
                    region=regionBackup;
                    drawGaussian(atoms[n].location(0), atoms[n].location(1),
                        atoms[n].location(2),region,atoms[n].intensity);
                    insertRegion(region);
                    currentRegionEval=evaluateRegion(region);
                    Nmovement++;
                }
            }
        }
        
        // Remove all the removed atoms
        for (int n=nmax-1; n>=0; n--)
            if (atoms[n].intensity==0)
                atoms.erase(atoms.begin()+n);

        drawApproximation();
        std::cout << "Iteration " << iter << " error= " << percentageDiff
                  << " Natoms= " << atoms.size()
                  << " Intensity= " << Nintensity
                  << " Location= " << Nmovement
                  << std::endl;
        
        if (iter>0)
            if ((oldError-percentageDiff)/oldError<stop) finished=true;
        iter++;
    } while (!finished);
}

/* Write ------------------------------------------------------------------- */
void Prog_Convert_Vol2Pseudo::writeResults()
{
    if (fnOut!="")
        Vcurrent.write(fnOut+".vol");

    // Compute the histogram of intensities
    Matrix1D<double> intensities;
    intensities.initZeros(atoms.size());
    FOR_ALL_ELEMENTS_IN_MATRIX1D(intensities)
        intensities(i)=atoms[i].intensity;
    histogram1D hist;
    compute_hist(intensities, hist, 100);
    hist.write(fnOut+".hist");

    // Save the difference
    VolumeXmipp Vdiff;
    Vdiff()=Vin()-Vcurrent();
    if (useMask && XSIZE(mask_prm.imask3D)!=0)
        FOR_ALL_ELEMENTS_IN_MATRIX3D(Vdiff())
            if (!mask_prm.imask3D(k,i,j))
                Vdiff(k,i,j)=0;
    Vdiff.write(fnOut+"_rawDiff.vol");
    
    Vdiff()/=range;
    Vdiff.write(fnOut+"_relativeDiff.vol");
    
    // Write the PDB
    double minIntensity=intensities.computeMin();
    double maxIntensity=intensities.computeMax();
    double a=1.000/(maxIntensity-minIntensity);
    
    FILE *fhOut=NULL;
    fhOut=fopen((fnOut+".pdb").c_str(),"w");
    if (!fhOut)
        REPORT_ERROR(1,(std::string)"Cannot open "+fnOut+".pdb for output");
    int nmax=atoms.size();
    int col=1;
    if (intensityColumn=="Bfactor") col=2;
    for (int n=0; n<nmax; n++)
    {
        double intensity=1.0;
        if (allowIntensity)
            intensity=ROUND(100*a*(atoms[n].intensity-minIntensity))/100.0;
        if (col==1)
            fprintf(fhOut,
                "ATOM  %5d DENS DENS    1    %8.3f%8.3f%8.3f%6.2f     1      DENS\n",
                n,
                (float)atoms[n].location(0),(float)atoms[n].location(1),
                (float)atoms[n].location(2),(float)intensity);
        else
            fprintf(fhOut,
                "ATOM  %5d DENS DENS    1    %8.3f%8.3f%8.3f     1%6.2f      DENS\n",
                n,
                (float)atoms[n].location(0),(float)atoms[n].location(1),
                (float)atoms[n].location(2),(float)intensity);
    }
    fclose(fhOut);
}

/* Run --------------------------------------------------------------------- */
void Prog_Convert_Vol2Pseudo::run()
{
    int iter=0;
    do
    {
        // Place seeds
        if (iter==0) placeSeeds(initialSeeds);
        else
        {
            double Natoms=atoms.size();
            removeSeeds(FLOOR(Natoms*(growSeeds/2)/100));
            placeSeeds(FLOOR(Natoms*growSeeds/100));
        }
        drawApproximation();
        if (iter==0)
            std::cout << "Initial error with " << atoms.size()
                      << " pseudo-atoms " << percentageDiff << std::endl;
        
        // Optimize seeds until convergence
        optimizeCurrentAtoms();
        std::cout << "Error with " << atoms.size() << " pseudo-atoms "
                  << percentageDiff << std::endl;
        writeResults();
        iter++;
    } while (percentageDiff>targetError);
    writeResults();
}
