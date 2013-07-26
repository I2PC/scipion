/***************************************************************************
 *
 * Authors:     Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
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

#include "pdb.h"
#include "fstream"
#include "args.h"
#include "matrix2d.h"
#include "mask.h"
#include "xmipp_fftw.h"
#include "integration.h"
#include "numerical_tools.h"

/* Atom charge ------------------------------------------------------------- */
int atomCharge(const std::string &atom)
{
    switch (atom[0])
    {
    case 'H':
        return 1;
        break;
    case 'C':
        return 6;
        break;
    case 'N':
        return 7;
        break;
    case 'O':
        return 8;
        break;
    case 'P':
        return 15;
        break;
    case 'S':
        return 16;
        break;
    case 'F': // Iron
        return 26;
        break;
    default:
        return 0;
    }
}

/* Atom radius ------------------------------------------------------------- */
double atomRadius(const std::string &atom)
{
    switch (atom[0])
    {
    case 'H':
        return 0.25;
        break;
    case 'C':
        return 0.70;
        break;
    case 'N':
        return 0.65;
        break;
    case 'O':
        return 0.60;
        break;
    case 'P':
        return 1.00;
        break;
    case 'S':
        return 1.00;
        break;
    case 'F': // Iron
        return 1.40;
        break;
    default:
        return 0;
    }
}

/* Compute geometry -------------------------------------------------------- */
void computePDBgeometry(const std::string &fnPDB,
                        Matrix1D<double> &centerOfMass,
                        Matrix1D<double> &limit0, Matrix1D<double> &limitF,
                        const std::string &intensityColumn)
{
    // Initialization
    centerOfMass.initZeros(3);
    limit0.initZeros(3);
    limitF.initZeros(3);
    limit0.initConstant(1e30);
    limitF.initConstant(-1e30);
    double total_mass = 0;

    // Open the file
    std::ifstream fh_pdb;
    fh_pdb.open(fnPDB.c_str());
    if (!fh_pdb)
        REPORT_ERROR(ERR_IO_NOTEXIST, fnPDB);

    // Process all lines of the file
    int col=1;
    if (intensityColumn=="Bfactor")
        col=2;
    while (!fh_pdb.eof())
    {
        // Read a ATOM line
        std::string line;
        getline(fh_pdb, line);
        if (line == "")
            continue;
        std::string kind = line.substr(0,4);
        if (kind != "ATOM" && kind!="HETA")
            continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        std::string atom_type = line.substr(13,2);
        double x = textToFloat(line.substr(30,8));
        double y = textToFloat(line.substr(38,8));
        double z = textToFloat(line.substr(46,8));

        // Update center of mass and limits
        if (x < XX(limit0))
            XX(limit0) = x;
        else if (x > XX(limitF))
            XX(limitF) = x;
        if (y < YY(limit0))
            YY(limit0) = y;
        else if (y > YY(limitF))
            YY(limitF) = y;
        if (z < ZZ(limit0))
            ZZ(limit0) = z;
        else if (z > ZZ(limitF))
            ZZ(limitF) = z;
        double weight;
        if (atom_type=="EN")
        {
            if      (col==1)
                weight=textToFloat(line.substr(54,6));
            else if (col==2)
                weight=textToFloat(line.substr(60,6));
        }
        else
        {
            if (kind=="HETA")
                continue;
            weight=(double) atomCharge(atom_type);
        }
        total_mass += weight;
        XX(centerOfMass) += weight * x;
        YY(centerOfMass) += weight * y;
        ZZ(centerOfMass) += weight * z;
    }

    // Finish calculations
    centerOfMass /= total_mass;

    // Close file
    fh_pdb.close();
}

/* Apply geometry ---------------------------------------------------------- */
void applyGeometryToPDBFile(const std::string &fn_in, const std::string &fn_out,
                   const Matrix2D<double> &A, bool centerPDB,
                   const std::string &intensityColumn)
{
    Matrix1D<double> centerOfMass, limit0, limitF;
    if (centerPDB)
    {
        computePDBgeometry(fn_in, centerOfMass,limit0, limitF,
                           intensityColumn);
        limit0 -= centerOfMass;
        limitF -= centerOfMass;
    }

    // Open files
    std::ifstream fh_in;
    fh_in.open(fn_in.c_str());
    if (!fh_in)
        REPORT_ERROR(ERR_IO_NOTEXIST, fn_in);
    std::ofstream fh_out;
    fh_out.open(fn_out.c_str());
    if (!fh_out)
        REPORT_ERROR(ERR_IO_NOWRITE, fn_out);

    // Process all lines of the file
    while (!fh_in.eof())
    {
        // Read an ATOM line
        std::string line;
        getline(fh_in, line);
        if (line == "")
        {
            fh_out << "\n";
            continue;
        }
        std::string kind = line.substr(0,4);
        if (kind != "ATOM" && kind != "HETA")
        {
            fh_out << line << std::endl;
            continue;
        }

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        double x = textToFloat(line.substr(30,8));
        double y = textToFloat(line.substr(38,8));
        double z = textToFloat(line.substr(46,8));

        Matrix1D<double> v(4);
        if (centerPDB)
        {
            VECTOR_R3(v,x-XX(centerOfMass),
                      y-YY(centerOfMass),z-ZZ(centerOfMass));
        }
        else
        {
            VECTOR_R3(v,x,y,z);
        }
        v(3)=1;
        v=A*v;

        char aux[15];
        sprintf(aux,"%8.3f",XX(v));
        line.replace(30,8,aux);
        sprintf(aux,"%8.3f",YY(v));
        line.replace(38,8,aux);
        sprintf(aux,"%8.3f",ZZ(v));
        line.replace(46,8,aux);

        fh_out << line << std::endl;
    }

    // Close files
    fh_in.close();
    fh_out.close();
}

/* Read phantom from PDB --------------------------------------------------- */
void PDBPhantom::read(const FileName &fnPDB)
{
    // Open file
    std::ifstream fh_in;
    fh_in.open(fnPDB.c_str());
    if (!fh_in)
        REPORT_ERROR(ERR_IO_NOTEXIST, fnPDB);

    // Process all lines of the file
    std::string line, kind;
    Atom atom;
    while (!fh_in.eof())
    {
        // Read an ATOM line
        getline(fh_in, line);
        if (line == "")
        {
            continue;
        }
        kind = line.substr(0,4);
        if (kind != "ATOM" && kind != "HETA")
            continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        atom.atomType = line[13];
        atom.x = textToFloat(line.substr(30,8));
        atom.y = textToFloat(line.substr(38,8));
        atom.z = textToFloat(line.substr(46,8));
        atomList.push_back(atom);
    }

    // Close files
    fh_in.close();
}

/* Shift ------------------------------------------------------------------- */
void PDBPhantom::shift(double x, double y, double z)
{
    int imax=atomList.size();
    for (int i=0; i<imax; i++)
    {
        atomList[i].x+=x;
        atomList[i].y+=y;
        atomList[i].z+=z;
    }
}

/* Read phantom from PDB --------------------------------------------------- */
void PDBRichPhantom::read(const FileName &fnPDB)
{
    // Open file
    std::ifstream fh_in;
    fh_in.open(fnPDB.c_str());
    if (!fh_in)
        REPORT_ERROR(ERR_IO_NOTEXIST, fnPDB);

    // Process all lines of the file
    std::string line, kind;
    RichAtom atom;
    while (!fh_in.eof())
    {
        // Read an ATOM line
        getline(fh_in, line);
        if (line == "")
        {
            continue;
        }
        kind = line.substr(0,4);
        if (kind == "ATOM" || kind == "HETA")
        {
			// Extract atom type and position
			// Typical line:
			// ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
			atom.name=line.substr(12,4);
			atom.atomType = line[13];
			atom.altloc=line[16];
			atom.resname=line.substr(17,3);
			atom.chainid=line[21];
			atom.resseq = textToInteger(line.substr(22,4));
			atom.icode = line[26];
			atom.x = textToFloat(line.substr(30,8));
			atom.y = textToFloat(line.substr(38,8));
			atom.z = textToFloat(line.substr(46,8));
			atom.occupancy = textToFloat(line.substr(54,6));
			atom.bfactor = textToFloat(line.substr(60,6));
			atomList.push_back(atom);
        } else if (kind == "REMA")
        	remarks.push_back(line);
    }

    // Close files
    fh_in.close();
}

/* Write phantom to PDB --------------------------------------------------- */
void PDBRichPhantom::write(const FileName &fnPDB)
{
    FILE* fh_out=fopen(fnPDB.c_str(),"w");
    if (!fh_out)
        REPORT_ERROR(ERR_IO_NOWRITE, fnPDB);
    size_t imax=remarks.size();
    for (size_t i=0; i<imax; ++i)
    	fprintf(fh_out,"%s\n",remarks[i].c_str());
    imax=atomList.size();
    for (size_t i=0; i<imax; ++i)
    {
    	const RichAtom &atom=atomList[i];
    	fprintf (fh_out,"ATOM  %5lu %4s%c%-4s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
    			(unsigned long int)i,atom.name.c_str(),
    			atom.altloc,atom.resname.c_str(),atom.chainid,
    			atom.resseq,atom.icode,atom.x,atom.y,atom.z,atom.occupancy,atom.bfactor,
    			atom.name.c_str());
    }
    fclose(fh_out);
}

/* Atom descriptors -------------------------------------------------------- */
void atomDescriptors(const std::string &atom, Matrix1D<double> &descriptors)
{
    descriptors.initZeros(11);
    if (atom=="H")
    {
        descriptors( 0)= 1;     // Z
        descriptors( 1)= 0.0088; // a1
        descriptors( 2)= 0.0449; // a2
        descriptors( 3)= 0.1481; // a3
        descriptors( 4)= 0.2356; // a4
        descriptors( 5)= 0.0914; // a5
        descriptors( 6)= 0.1152; // b1
        descriptors( 7)= 1.0867; // b2
        descriptors( 8)= 4.9755; // b3
        descriptors( 9)=16.5591; // b4
        descriptors(10)=43.2743; // b5
    }
    else if (atom=="C")
    {
        descriptors( 0)= 6;     // Z
        descriptors( 1)= 0.0489; // a1
        descriptors( 2)= 0.2091; // a2
        descriptors( 3)= 0.7537; // a3
        descriptors( 4)= 1.1420; // a4
        descriptors( 5)= 0.3555; // a5
        descriptors( 6)= 0.1140; // b1
        descriptors( 7)= 1.0825; // b2
        descriptors( 8)= 5.4281; // b3
        descriptors( 9)=17.8811; // b4
        descriptors(10)=51.1341; // b5
    }
    else if (atom=="N")
    {
        descriptors( 0)= 7;     // Z
        descriptors( 1)= 0.0267; // a1
        descriptors( 2)= 0.1328; // a2
        descriptors( 3)= 0.5301; // a3
        descriptors( 4)= 1.1020; // a4
        descriptors( 5)= 0.4215; // a5
        descriptors( 6)= 0.0541; // b1
        descriptors( 7)= 0.5165; // b2
        descriptors( 8)= 2.8207; // b3
        descriptors( 9)=10.6297; // b4
        descriptors(10)=34.3764; // b5
    }
    else if (atom=="O")
    {
        descriptors( 0)= 8;     // Z
        descriptors( 1)= 0.0365; // a1
        descriptors( 2)= 0.1729; // a2
        descriptors( 3)= 0.5805; // a3
        descriptors( 4)= 0.8814; // a4
        descriptors( 5)= 0.3121; // a5
        descriptors( 6)= 0.0652; // b1
        descriptors( 7)= 0.6184; // b2
        descriptors( 8)= 2.9449; // b3
        descriptors( 9)= 9.6298; // b4
        descriptors(10)=28.2194; // b5
    }
    else if (atom=="P")
    {
        descriptors( 0)=15;     // Z
        descriptors( 1)= 0.1005; // a1
        descriptors( 2)= 0.4615; // a2
        descriptors( 3)= 1.0663; // a3
        descriptors( 4)= 2.5854; // a4
        descriptors( 5)= 1.2725; // a5
        descriptors( 6)= 0.0977; // b1
        descriptors( 7)= 0.9084; // b2
        descriptors( 8)= 4.9654; // b3
        descriptors( 9)=18.5471; // b4
        descriptors(10)=54.3648; // b5
    }
    else if (atom=="S")
    {
        descriptors( 0)=16;     // Z
        descriptors( 1)= 0.0915; // a1
        descriptors( 2)= 0.4312; // a2
        descriptors( 3)= 1.0847; // a3
        descriptors( 4)= 2.4671; // a4
        descriptors( 5)= 1.0852; // a5
        descriptors( 6)= 0.0838; // b1
        descriptors( 7)= 0.7788; // b2
        descriptors( 8)= 4.3462; // b3
        descriptors( 9)=15.5846; // b4
        descriptors(10)=44.6365; // b5
    }
    else if (atom=="Fe")
    {
        descriptors( 0)=26;     // Z
        descriptors( 1)= 0.1929; // a1
        descriptors( 2)= 0.8239; // a2
        descriptors( 3)= 1.8689; // a3
        descriptors( 4)= 2.3694; // a4
        descriptors( 5)= 1.9060; // a5
        descriptors( 6)= 0.1087; // b1
        descriptors( 7)= 1.0806; // b2
        descriptors( 8)= 4.7637; // b3
        descriptors( 9)=22.8500; // b4
        descriptors(10)=76.7309; // b5
    }
    else
        REPORT_ERROR(ERR_VALUE_INCORRECT,(std::string)"atomDescriptors: Unknown atom "+atom);
}

/* Electron form factor in Fourier ----------------------------------------- */
double electronFormFactorFourier(double f, const Matrix1D<double> &descriptors)
{
    double retval=0;
    for (int i=1; i<=5; i++)
    {
        double ai=VEC_ELEM(descriptors,i);
        double bi=VEC_ELEM(descriptors,i+5);
        retval+=ai*exp(-bi*f*f);
    }
    return retval;
}

/* Electron form factor in real space -------------------------------------- */
/* We know the transform pair

   sqrt(pi/b)*exp(-x^2/(4*b)) <----> exp(-b*W^2)
   
   We also know that 
   
   X(f)=sum_i ai exp(-bi*f^2)
   
   Therefore, using W=2*pi*f
   
   X(W)=sum_i ai exp(-bi/(2*pi)^2*W^2)
   
   Thus, the actual b for the inverse Fourier transform is bi/(2*pi)^2
   And we have to divide by 2*pi to account for the Jacobian of the
   transformation.
*/
double electronFormFactorRealSpace(double r,
                                   const Matrix1D<double> &descriptors)
{
    double retval=0;
    for (int i=1; i<=5; i++)
    {
        double ai=descriptors(i);
        double bi=descriptors(i+5);
        double b=bi/(4*PI*PI);
        retval+=ai*sqrt(PI/b)*exp(-r*r/(4*b));
    }
    retval/=2*PI;
    return retval;
}

/* Computation of the low pass filter -------------------------------------- */
// Returns the impulse response of the lowpass filter
void hlpf(MultidimArray<double> &f, int M, double T, const std::string &filterType,
          MultidimArray<double> &filter, double reductionFactor=0.8,
          double ripple=0.01, double deltaw=1.0/8.0)
{
    filter.initZeros(XSIZE(f));
    filter.setXmippOrigin();

    int Nmax=(int)CEIL(M/2.0);
    if (filterType=="SimpleAveraging")
    {
        FOR_ALL_ELEMENTS_IN_ARRAY1D(filter)
        if (ABS(i)<=Nmax)
            filter(i)=1.0/(2*Nmax+1);
    }
    else if (filterType=="SincKaiser")
    {
        SincKaiserMask(filter,reductionFactor*PI/M,ripple,deltaw);
        filter/=filter.sum();
        if (FINISHINGX(f)>FINISHINGX(filter))
            filter.selfWindow(STARTINGX(f),FINISHINGX(f));
        else
            f.selfWindow(STARTINGX(filter),FINISHINGX(filter));
    }
}

/* Convolution between f and the hlpf -------------------------------------- */
void fhlpf(const MultidimArray<double> &f, const MultidimArray<double> &filter,
           int M, MultidimArray<double> &convolution)
{
    // Expand the two input signals
    int Nmax=FINISHINGX(filter);
    MultidimArray<double> auxF, auxFilter;
    auxF=f;
    auxFilter=filter;
    auxF.selfWindow(STARTINGX(f)-Nmax,FINISHINGX(f)+Nmax);
    auxFilter.selfWindow(STARTINGX(filter)-Nmax,FINISHINGX(filter)+Nmax);

    // Convolve in Fourier
    MultidimArray< std::complex<double> > F, Filter;
    FourierTransform(auxF,F);
    FourierTransform(auxFilter,Filter);
    F*=Filter;

    // Correction for the double phase factor and the double
    // amplitude factor
    const double K1=2*PI*(STARTINGX(auxFilter)-1);
    const double K2=XSIZE(auxFilter);
    std::complex<double> aux;
    double * ptrAux=(double*)&aux;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(F)
    {
        double w;
        FFT_IDX2DIGFREQ(i,XSIZE(F),w);
        double arg=w*K1;
        sincos(arg,ptrAux+1,ptrAux);
        *ptrAux*=K2;
        *(ptrAux+1)*=K2;
        A1D_ELEM(F,i)*=aux;
    }
    InverseFourierTransform(F,convolution);
    convolution.setXmippOrigin();
}

/* Optimization of the low pass filter to fit a given atom ----------------- */
Matrix1D<double> globalHlpfPrm(3);
MultidimArray<double> globalf;
int globalM;
double globalT;
std::string globalAtom;

double Hlpf_fitness(double *p, void *prm)
{
    double reductionFactor=p[1];
    double ripple=p[2];
    double deltaw=p[3];

    if (reductionFactor<0.7 || reductionFactor>1.3)
        return 1e38;
    if (ripple<0 || ripple>0.2)
        return 1e38;
    if (deltaw<0 || deltaw>0.2)
        return 1e38;

    // Construct the filter with the current parameters
    MultidimArray<double> filter, auxf;
    auxf=globalf;
    hlpf(auxf, globalM, globalT, "SincKaiser", filter, reductionFactor,
         ripple, deltaw);

    // Convolve the filter with the atomic profile
    MultidimArray<double> fhlpfFinelySampled;
    fhlpf(auxf, filter, globalM, fhlpfFinelySampled);

    // Coarsely sample
    double Rmax=FINISHINGX(fhlpfFinelySampled)*globalT;
    int imax=CEIL(Rmax/(globalM*globalT));
    MultidimArray<double> fhlpfCoarselySampled(2*imax+1);
    MultidimArray<double> splineCoeffsfhlpfFinelySampled;
    produceSplineCoefficients(BSPLINE3,splineCoeffsfhlpfFinelySampled,fhlpfFinelySampled);
    fhlpfCoarselySampled.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY1D(fhlpfCoarselySampled)
    {
        double r=i*(globalM*globalT)/globalT;
        fhlpfCoarselySampled(i)=
            splineCoeffsfhlpfFinelySampled.interpolatedElementBSpline1D(r,3);
    }

    // Build the frequency response of the convolved and coarsely sampled
    // atom
    MultidimArray<double> aux, FfilterMag, freq;
    MultidimArray< std::complex<double> > Ffilter;
    aux=fhlpfCoarselySampled;
    aux.selfWindow(-10*FINISHINGX(aux),10*FINISHINGX(aux));
    FourierTransform(aux,Ffilter);
    FFT_magnitude(Ffilter,FfilterMag);
    freq.initZeros(XSIZE(Ffilter));

    FOR_ALL_ELEMENTS_IN_ARRAY1D(FfilterMag)
    FFT_IDX2DIGFREQ(i,XSIZE(FfilterMag),freq(i));
    freq/=globalM*globalT;
    double amplitudeFactor=fhlpfFinelySampled.sum()/
                           fhlpfCoarselySampled.sum();

    // Compute the error in representation
    double error=0;
    Matrix1D<double> descriptors;
    atomDescriptors(globalAtom, descriptors);
    double iglobalT=1.0/globalT;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(FfilterMag)
    if (A1D_ELEM(freq,i)>=0)
    {
        double f1=log10(A1D_ELEM(FfilterMag,i)*XSIZE(FfilterMag)*amplitudeFactor);
        double f2=log10(iglobalT*
                           electronFormFactorFourier(A1D_ELEM(freq,i),descriptors));
        double diff=20*(f1-f2);
        error+=diff*diff;
    }

    return error/XSIZE(FfilterMag);
}

/** Optimize the low pass filter.
    The optimization is so that the Fourier response of the coarsely
    downsampled and convolved atom profile resembles as much as possible
    the ideal atomic response up to the maximum frequency provided by the
    Nyquist frequency associated to M*T.
    
    f is the electron scattering factor in real space sampled at a sampling
    rate T. M is the downsampling factor. atom is the name of the atom
    being optimized. filter is an output parameter with the optimal impulse
    response sampled at a sampling rate T. bestPrm(0)=reduction function
    of the cutoff frequency, bestPrm(1)=ripple of the Kaiser selfWindow,
    bestPrm(2)=deltaw of the Kaiser selfWindow.
*/
void optimizeHlpf(MultidimArray<double> &f, int M, double T, const std::string &atom,
		MultidimArray<double> &filter, Matrix1D<double> &bestPrm)
{
    globalHlpfPrm(0)=1.0;     // reduction factor
    globalHlpfPrm(1)=0.01;    // ripple
    globalHlpfPrm(2)=1.0/8.0; // deltaw
    globalf=f;
    globalM=M;
    globalT=T;
    globalAtom=atom;
    double fitness;
    int iter;
    Matrix1D<double> steps(3);
    steps.initConstant(1);
    powellOptimizer(globalHlpfPrm, 1, 3,
                    &Hlpf_fitness, NULL, 0.05, fitness, iter, steps, false);
    bestPrm=globalHlpfPrm;
    hlpf(f, M, T, "SincKaiser", filter, bestPrm(0), bestPrm(1), bestPrm(2));
}

/* Atom radial profile ----------------------------------------------------- */
void atomRadialProfile(int M, double T, const std::string &atom,
		MultidimArray<double> &profile)
{
    // Compute the electron form factor in real space
    double largestb1=76.7309/(4*PI*PI);
    double Rmax=4*sqrt(2*largestb1);
    int imax=(int)CEIL(Rmax/T);
    Matrix1D<double> descriptors;
    atomDescriptors(atom, descriptors);
    MultidimArray<double> f(2*imax+1);
    f.setXmippOrigin();
    for (int i=-imax; i<=imax; i++)
    {
        double r=i*T;
        f(i)=electronFormFactorRealSpace(r,descriptors);
    }

    // Compute the optimal filter
    MultidimArray<double> filter;
    Matrix1D<double> bestPrm;
    optimizeHlpf(f, M, T, atom, filter, bestPrm);

    // Perform the convolution
    fhlpf(f, filter, M, profile);

    // Remove zero values
    int ileft=STARTINGX(profile);
    if (fabs(profile(ileft))<1e-3)
        for (ileft=STARTINGX(profile)+1; ileft<=0; ileft++)
            if (fabs(profile(ileft))>1e-3)
                break;
    int iright=FINISHINGX(profile);
    if (fabs(profile(iright))<1e-3)
        for (iright=FINISHINGX(profile)-1; iright>=0; iright--)
            if (fabs(profile(iright))>1e-3)
                break;
    profile.selfWindow(ileft,iright);
}

/** Atom projection profile ------------------------------------------------ */
class AtomValueFunc: public doubleFunction
{
public:
    int M;
    double r0_2, z;
    const MultidimArray<double> *profileCoefficients;
    virtual double operator()()
    {
        double r=M*sqrt(r0_2+z*z);
        if (ABS(r)>FINISHINGX(*profileCoefficients))
            return 0;
        return profileCoefficients->interpolatedElementBSpline1D(r,3);
    }
};

#define INTEGRATION 2
void atomProjectionRadialProfile(int M,
                                 const MultidimArray<double> &profileCoefficients,
                                 MultidimArray<double> &projectionProfile)
{
    AtomValueFunc atomValue;
    atomValue.profileCoefficients=&profileCoefficients;
    atomValue.M=M;
    double radius=(double)FINISHINGX(profileCoefficients)/M;
    double r2=radius*radius;

    projectionProfile.initZeros(profileCoefficients);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(projectionProfile)
    {
        double r0=(double)i/M;
        atomValue.r0_2=r0*r0;
        if (atomValue.r0_2>r2)
            continue; // Because of numerical instabilities

        double maxZ=sqrt(r2-atomValue.r0_2);

#if INTEGRATION==1

        double dz=1/24.0;
        double integral=0;
        for (atomValue.z=-maxZ; atomValue.z<=maxZ; atomValue.z+=dz)
        {
            projectionProfile(i)+=atomValue();
        }
        projectionProfile(i)*=dz;
#else

        Romberg Rom(atomValue, atomValue.z,-maxZ,maxZ);
        projectionProfile(i) = Rom();
#endif

    }
}

/** Atom interpolations ---------------------------------------------------- */
void AtomInterpolator::setup(int m, double hights, bool computeProjection)
{
    M=m;
    highTs=hights;
    if (volumeProfileCoefficients.size()==7)
    	return;
    addAtom("H",computeProjection);
    addAtom("C",computeProjection);
    addAtom("N",computeProjection);
    addAtom("O",computeProjection);
    addAtom("P",computeProjection);
    addAtom("S",computeProjection);
    addAtom("Fe",computeProjection);
}

void AtomInterpolator::addAtom(const std::string &atom, bool computeProjection)
{
    MultidimArray<double> profile;
    MultidimArray<double> splineCoeffs;

    // Atomic profile
    atomRadialProfile(M, highTs, atom, profile);
    produceSplineCoefficients(BSPLINE3,splineCoeffs,profile);
    volumeProfileCoefficients.push_back(splineCoeffs);

    // Radius
    radii.push_back((double)FINISHINGX(profile)/M);

    // Projection profile
    if (computeProjection)
    {
        atomProjectionRadialProfile(M, splineCoeffs, profile);
        produceSplineCoefficients(BSPLINE3,splineCoeffs,profile);
        projectionProfileCoefficients.push_back(splineCoeffs);
    }
}

/** PDB projection --------------------------------------------------------- */
// Taken from phantom.cpp, Feature::project_to
void projectAtom(const Atom &atom, Projection &P,
                 const Matrix2D<double> &VP, const Matrix2D<double> &PV,
                 const AtomInterpolator &interpolator)
{
#define SUBSAMPLING 2                  // for every measure 2x2 line
    // integrals will be taken to
    // avoid numerical errors
#define SUBSTEP 1/(SUBSAMPLING*2.0)

    Matrix1D<double> origin(3);
    Matrix1D<double> direction;
    VP.getRow(2, direction);
    direction.selfTranspose();
    Matrix1D<double> corner1(3), corner2(3);
    Matrix1D<double> act(3);
    SPEED_UP_temps012;

    // Find center of the feature in the projection plane ...................
    // Step 1). Project the center to the plane, the result is in the
    //          universal coord system
    Matrix1D<double> Center(3);
    VECTOR_R3(Center, atom.x, atom.y, atom.z);
    double max_distance=interpolator.atomRadius(atom.atomType);
    M3x3_BY_V3x1(origin, VP, Center);

    //#define DEBUG_LITTLE
#ifdef DEBUG_LITTLE

    std::cout << "Actual atom\n"        << atom.atomType << " ("
    << atom.x << "," << atom.y << "," << atom.z << ")\n";
    std::cout << "Center              " << Center.transpose() << std::endl;
    std::cout << "VP matrix\n"          << VP << std::endl;
    std::cout << "P.direction         " << P.direction.transpose() << std::endl;
    std::cout << "direction           " << direction.transpose() << std::endl;
    std::cout << "P.euler matrix      " << P.euler << std::endl;
    std::cout << "max_distance        " << max_distance << std::endl;
    std::cout << "origin              " << origin.transpose() << std::endl;
#endif

    // Find limits for projection ...........................................
    // Choose corners for the projection of this feature. It is supposed
    // to have at the worst case a projection of size max_distance
    VECTOR_R3(corner1, max_distance, max_distance, max_distance);
    VECTOR_R3(corner2, -max_distance, -max_distance, -max_distance);
#ifdef DEBUG_LITTLE

    std::cout << "Corner1 : " << corner1.transpose() << std::endl
    << "Corner2 : " << corner2.transpose() << std::endl;
#endif

    box_enclosing(corner1, corner2, VP, corner1, corner2);
#ifdef DEBUG_LITTLE

    std::cout << "Corner1 moves to : " << corner1.transpose() << std::endl
    << "Corner2 moves to : " << corner2.transpose() << std::endl;
#endif

    V3_PLUS_V3(corner1, origin, corner1);
    V3_PLUS_V3(corner2, origin, corner2);
#ifdef DEBUG_LITTLE

    std::cout << "Corner1 finally is : " << corner1.transpose() << std::endl
    << "Corner2 finally is : " << corner2.transpose() << std::endl;
#endif
    // Discard not necessary components
    corner1.resize(2);
    corner2.resize(2);

    // Clip to image size
    sortTwoVectors(corner1, corner2);
    XX(corner1) = CLIP(ROUND(XX(corner1)), STARTINGX(P()), FINISHINGX(P()));
    YY(corner1) = CLIP(ROUND(YY(corner1)), STARTINGY(P()), FINISHINGY(P()));
    XX(corner2) = CLIP(ROUND(XX(corner2)), STARTINGX(P()), FINISHINGX(P()));
    YY(corner2) = CLIP(ROUND(YY(corner2)), STARTINGY(P()), FINISHINGY(P()));

#ifdef DEBUG_LITTLE

    std::cout << "corner1      " << corner1.transpose() << std::endl;
    std::cout << "corner2      " << corner2.transpose() << std::endl;
    std::cout.flush();
#endif

    // Check if there is something to project
    if (XX(corner1) == XX(corner2))
        return;
    if (YY(corner1) == YY(corner2))
        return;

    // Study the projection for each point in the projection plane ..........
    // (u,v) are in the deformed projection plane (if any deformation)
    for (int v = (int)YY(corner1); v <= (int)YY(corner2); v++)
        for (int u = (int)XX(corner1); u <= (int)XX(corner2); u++)
        {
            double length = 0;
            //#define DEBUG_EVEN_MORE
#ifdef DEBUG_EVEN_MORE

            std::cout << "Studying point (" << u << "," << v << ")\n";
            std::cout.flush();
#endif

            // Perform subsampling ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
            double u0 = u - (int)(SUBSAMPLING / 2.0) * SUBSTEP;
            double v0 = v - (int)(SUBSAMPLING / 2.0) * SUBSTEP;
            double actv = v0;
            for (int subv = 0; subv < SUBSAMPLING; subv++)
            {
                double actu = u0;
                for (int subu = 0; subu < SUBSAMPLING; subu++)
                {
                    // Compute the coordinates of point (subu,subv) which is
                    // within the plane in the universal coordinate system
                    XX(act) = actu;
                    YY(act) = actv;
                    ZZ(act) = 0;
                    M3x3_BY_V3x1(act, PV, act);

                    // Compute the intersection of a ray which passes through
                    // this point and its direction is perpendicular to the
                    // projection plane
                    double r=point_line_distance_3D(Center, act, direction);
                    double possible_length=interpolator.
                                           projectionAtDistance(atom.atomType,r);
                    if (possible_length > 0)
                        length += possible_length;

#ifdef DEBUG_EVEN_MORE

                    std::cout << "Averaging at (" << actu << "," << actv << ")\n";
                    std::cout << "   which in univ. coords is " << act.transpose() << std::endl;
                    std::cout << "   r=" << r << std::endl;
                    std::cout << "   intersection there " << possible_length << std::endl;
                    std::cout.flush();
#endif
                    // Prepare for next iteration
                    actu += SUBSTEP * 2.0;
                }
                actv += SUBSTEP * 2.0;
            }
            length /= (SUBSAMPLING * SUBSAMPLING);
            //#define DEBUG
#ifdef DEBUG

            std::cout << "Final value added at position (" << u << "," << v << ")="
            << length << std::endl;
            std::cout.flush();
#endif

            // Add at the correspondant pixel the found intersection ,,,,,,,,,,
            IMGPIXEL(P, v, u) += length;
        }
}
#undef DEBUG
#undef DEBUG_LITTLE
#undef DEBUG_EVEN_MORE

void projectPDB(const PDBPhantom &phantomPDB,
                const AtomInterpolator &interpolator, Projection &proj,
                int Ydim, int Xdim, double rot, double tilt, double psi)
{
    // Initialise projection
    proj().initZeros(Ydim, Xdim);
    proj().setXmippOrigin();
    proj.setAngles(rot, tilt, psi);

    // Compute volume to Projection matrix
    Matrix2D<double> VP = proj.euler;
    Matrix2D<double> PV = VP.inv();

    // Project all elements
    for (size_t i = 0; i < phantomPDB.getNumberOfAtoms(); i++)
    {
        try
        {
            projectAtom(phantomPDB.getAtom(i), proj, VP, PV, interpolator);
        }
        catch (XmippError XE) {}
    }
}

void distanceHistogramPDB(const PDBPhantom &phantomPDB, size_t Nnearest, int Nbins, Histogram1D &hist)
{
    // Compute the histogram of distances
	const std::vector<Atom> &atoms=phantomPDB.atomList;
    int Natoms=atoms.size();
    MultidimArray<double> NnearestDistances;
    NnearestDistances.resize((Natoms-1)*Nnearest);
    for (int i=0; i<Natoms; i++)
    {
        std::vector<double> NnearestToThisAtom;
        const Atom& atom_i=atoms[i];
        for (int j=i+1; j<Natoms; j++)
        {
            const Atom& atom_j=atoms[j];
            double diffx=atom_i.x-atom_j.x;
            double diffy=atom_i.y-atom_j.y;
            double diffz=atom_i.z-atom_j.z;
            double dist=sqrt(diffx*diffx+diffy*diffy+diffz*diffz);
        	//std::cout << "Analyzing " << i << " and " << j << " -> d=" << dist << std::endl;
            size_t nearestSoFar=NnearestToThisAtom.size();
            if (nearestSoFar==0)
            {
                NnearestToThisAtom.push_back(dist);
            	//std::cout << "Pushing d" << std::endl;
            }
            else
            {
                size_t idx=0;
                while (idx<nearestSoFar && NnearestToThisAtom[idx]<dist)
                    idx++;
                if (idx<nearestSoFar)
                {
                    NnearestToThisAtom.insert(NnearestToThisAtom.begin()+idx,1,dist);
                    if (NnearestToThisAtom.size()>Nnearest)
                        NnearestToThisAtom.erase(NnearestToThisAtom.begin()+Nnearest);
                }
                if (idx==nearestSoFar && nearestSoFar<Nnearest)
                {
                    NnearestToThisAtom.push_back(dist);
                	//std::cout << "Pushing d" << std::endl;
                }
            }
        }
		if (i<Natoms-1)
			for (size_t k=0; k<Nnearest; k++)
				NnearestDistances(i*Nnearest+k)=NnearestToThisAtom[k];
    }
    compute_hist(NnearestDistances, hist, 0, NnearestDistances.computeMax(), Nbins);
}
