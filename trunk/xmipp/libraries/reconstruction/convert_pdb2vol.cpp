/***************************************************************************
 *
 * Authors:
 *
 * Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include "convert_pdb2vol.h"
#include "fourier_filter.h"
#include "blobs.h"

#include <data/args.h>
#include <data/fft.h>

#include <fstream>

/* Empty constructor ------------------------------------------------------- */
Prog_PDBPhantom_Parameters::Prog_PDBPhantom_Parameters()
{
    blob.radius = 2;   // Blob radius in voxels
    blob.order  = 2;   // Order of the Bessel function
    blob.alpha  = 3.6; // Smoothness parameter
    output_dim = -1;
    fn_pdb = "";
    Ts = 1;
    highTs = 1.0/12.0;
    useBlobs=false;

    // Periodic table for the blobs
    periodicTable.resize(7, 2);
    periodicTable(0, 0) = 0.25;
    periodicTable(0, 1) = 1;  // Hydrogen
    periodicTable(1, 0) = 0.70;
    periodicTable(1, 1) = 6; // Carbon
    periodicTable(2, 0) = 0.65;
    periodicTable(2, 1) = 7; // Nitrogen
    periodicTable(3, 0) = 0.60;
    periodicTable(3, 1) = 8; // Oxygen
    periodicTable(4, 0) = 1.00;
    periodicTable(4, 1) = 15; // Phosphorus
    periodicTable(5, 0) = 1.00;
    periodicTable(5, 1) = 16; // Sulfur
    periodicTable(6, 0) = 1.40;
    periodicTable(6, 1) = 26; // Iron
    // Correct the atom weights by the blob weight
    for (int i = 0; i < YSIZE(periodicTable); i++)
    {
        periodicTable(i, 1) /=
            basvolume(periodicTable(i, 0) / highTs, blob.alpha, blob.order, 3);
    }
}

/* Produce Side Info ------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::produceSideInfoAtom(
    const string &atom) {
    Matrix1D<double> profile;
    Matrix1D<double> *splineCoeffs=NULL;
    atomRadialProfile(M, highTs, atom, profile);
    splineCoeffs=new Matrix1D<double>; 
    profile.produceSplineCoefficients(*splineCoeffs,3);
    atomProfiles.push_back(splineCoeffs);
}

void Prog_PDBPhantom_Parameters::produceSideInfo() {
    if (!useBlobs) {
	// Compute the downsampling factor
	M=(int)ROUND(Ts/highTs);

	// Atom profiles for the electron scattering method
	produceSideInfoAtom("H");
	produceSideInfoAtom("C");
	produceSideInfoAtom("N");
	produceSideInfoAtom("O");
	produceSideInfoAtom("P");
	produceSideInfoAtom("S");
	produceSideInfoAtom("Fe");
    }
}

/* Atom description ------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::atomBlobDescription(const string &_element,
        double &weight, double &radius) const
{
    int idx = -1;
    weight = radius = 0;
    switch (_element[0])
    {
    case 'H':
        idx = 0;
        break;
    case 'C':
        idx = 1;
        break;
    case 'N':
        idx = 2;
        break;
    case 'O':
        idx = 3;
        break;
    case 'P':
        idx = 4;
        break;
    case 'S':
        idx = 5;
        break;
    case 'F':
        idx = 6;
        break;
    default:
        std::cout << "Unknown :" << _element << std::endl;
        return;
    }
    radius = periodicTable(idx, 0);
    weight = periodicTable(idx, 1);
}

/* Read parameters --------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::read(int argc, char **argv)
{
    fn_pdb = getParameter(argc, argv, "-i");
    fn_out = getParameter(argc, argv, "-o", "");
    if (fn_out == "") fn_out = fn_pdb.without_extension();
    Ts = textToFloat(getParameter(argc, argv, "-sampling_rate", "1"));
    highTs = textToFloat(getParameter(argc, argv, "-high_sampling_rate", "0.08333333"));
    output_dim = textToInteger(getParameter(argc, argv, "-size", "-1"));
    useBlobs = checkParameter(argc, argv, "-blobs");
}

/* Usage ------------------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::usage()
{
    std::cerr << "convert_pdb2vol\n"
	      << "   -i <pdb file>		       : File to process\n"
	      << "  [-o <fn_root>]		       : Root name for output\n"
	      << "  [-sampling_rate <Ts=1>]	       : Sampling rate (Angstroms/pixel)\n"
	      << "  [-high_sampling_rate <highTs=1/12>]: Sampling rate before downsampling\n"
	      << "  [-size <output_dim>]	       : Final size in pixels (must be a power of 2)\n"
	      << "  [-blobs]                           : Use blobs instead of scattering factors\n"
	      << "\n"
	      << "Example of use: Sample at 1.6A and limit the frequency to 10A\n"
	      << "   xmipp_convert_pdb2vol -i 1o7d.pdb -sampling_rate 1.6\n"
	      << "   xmipp_fourier_filter -i 1o7d.vol -o 1o7d_filtered.vol -low_pass 10 -sampling 1.6 -fourier_mask raised_cosine 0.1\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::show()
{
    std::cout << "PDB file:           " << fn_pdb     << std::endl
              << "Sampling rate:      " << Ts         << std::endl
              << "High sampling rate: " << highTs     << std::endl
              << "Size:               " << output_dim << std::endl
	      << "Use blobs:          " << useBlobs   << std::endl
    ;
}

/* Compute protein geometry ------------------------------------------------ */
void Prog_PDBPhantom_Parameters::compute_protein_geometry()
{
    // Initialization
    centerOfMass.initZeros(3);
    Matrix1D<double> limit0(3), limitF(3);
    limit0.init_constant(1e30);
    limitF.init_constant(-1e30);
    double total_mass = 0;

    // Open the file
    std::ifstream fh_pdb;
    fh_pdb.open(fn_pdb.c_str());
    if (!fh_pdb)
        REPORT_ERROR(1, (string)"Prog_PDBPhantom_Parameters::protein_geometry:"
                     "Cannot open " + fn_pdb + " for reading");

    // Process all lines of the file�
    while (!fh_pdb.eof())
    {
        // Read a ATOM line
        string line;
        getline(fh_pdb, line);
        if (line == "") continue;
        string kind = firstToken(line);
        if (kind != "ATOM") continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        string dummy = nextToken();
        string atom_type = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        double x = textToFloat(nextToken());
        double y = textToFloat(nextToken());
        double z = textToFloat(nextToken());

        // Update center of mass and limits
        if (x < XX(limit0)) XX(limit0) = x;
        else if (x > XX(limitF)) XX(limitF) = x;
        if (y < YY(limit0)) YY(limit0) = y;
        else if (y > YY(limitF)) YY(limitF) = y;
        if (z < ZZ(limit0)) ZZ(limit0) = z;
        else if (z > ZZ(limitF)) ZZ(limitF) = z;
        double weight, radius;
        atomBlobDescription(atom_type, weight, radius);
        total_mass += weight;
        XX(centerOfMass) += weight * x;
        YY(centerOfMass) += weight * y;
        ZZ(centerOfMass) += weight * z;
    }

    // Finish calculations
    centerOfMass /= total_mass;
    limit0 -= centerOfMass;
    limitF -= centerOfMass;
    limit.resize(3);
    XX(limit) = MAX(ABS(XX(limit0)), ABS(XX(limitF)));
    YY(limit) = MAX(ABS(YY(limit0)), ABS(YY(limitF)));
    ZZ(limit) = MAX(ABS(ZZ(limit0)), ABS(ZZ(limitF)));

    // Update output size if necessary
    if (output_dim == -1)
    {
        int max_dim = MAX(CEIL(ZZ(limit) * 2 / Ts) + 5, CEIL(YY(limit) * 2 / Ts) + 5);
        max_dim = MAX(max_dim, CEIL(XX(limit) * 2 / Ts) + 5);
        output_dim = (int)NEXT_POWER_OF_2(max_dim);
        std::cout << "Setting output_dim to " << output_dim << std::endl;
    }

    // Close file
    fh_pdb.close();
}

/* Create protein at a high sampling rate ---------------------------------- */
void Prog_PDBPhantom_Parameters::create_protein_at_high_sampling_rate()
{
    // Create an empty volume to hold the protein
    Vhigh().initZeros((int)NEXT_POWER_OF_2(output_dim / highTs),
                      (int)NEXT_POWER_OF_2(output_dim / highTs),
                      (int)NEXT_POWER_OF_2(output_dim / highTs));
    Vhigh().setXmippOrigin();
    std::cout << "The highly sampled volume is of size " << XSIZE(Vhigh())
              << std::endl;

    // Fill the volume with the different atoms
    std::ifstream fh_pdb;
    fh_pdb.open(fn_pdb.c_str());
    if (!fh_pdb)
        REPORT_ERROR(1, (string)"Prog_PDBPhantom_Parameters::protein_geometry:"
                     "Cannot open " + fn_pdb + " for reading");

    // Process all lines of the file
    while (!fh_pdb.eof())
    {
        // Read an ATOM line
        string line;
        getline(fh_pdb, line);
        if (line == "") continue;
        string kind = firstToken(line);
        if (kind != "ATOM") continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        string dummy = nextToken();
        string atom_type = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        double x = textToFloat(nextToken());
        double y = textToFloat(nextToken());
        double z = textToFloat(nextToken());

        // Correct position
        Matrix1D<double> r(3);
        VECTOR_R3(r, x, y, z);
        r -= centerOfMass;
        r /= highTs;

        // Characterize atom
        double weight, radius;
        atomBlobDescription(atom_type, weight, radius);
        blob.radius = radius;

        // Find the part of the volume that must be updated
        int k0 = MAX(FLOOR(ZZ(r) - radius), STARTINGZ(Vhigh()));
        int kF = MIN(CEIL(ZZ(r) + radius), FINISHINGZ(Vhigh()));
        int i0 = MAX(FLOOR(YY(r) - radius), STARTINGY(Vhigh()));
        int iF = MIN(CEIL(YY(r) + radius), FINISHINGY(Vhigh()));
        int j0 = MAX(FLOOR(XX(r) - radius), STARTINGX(Vhigh()));
        int jF = MIN(CEIL(XX(r) + radius), FINISHINGX(Vhigh()));

        // Fill the volume with this atom
        for (int k = k0; k <= kF; k++)
            for (int i = i0; i <= iF; i++)
                for (int j = j0; j <= jF; j++)
                {
                    Matrix1D<double> rdiff(3);
                    VECTOR_R3(rdiff, XX(r) - j, YY(r) - i, ZZ(r) - k);
                    Vhigh(k, i, j) += weight * blob_val(rdiff.module(), blob);
                }
    }

    // Close file
    fh_pdb.close();
}

/* Create protein at a high sampling rate ---------------------------------- */
void Prog_PDBPhantom_Parameters::create_protein_at_low_sampling_rate()
{
    // Compute the integer downsapling factor
    int M = FLOOR(Ts / highTs);
    double current_Ts = highTs;

    // Use Bsplines pyramid if possible
    int levels = FLOOR(log10((double)M) / log10(2.0) + XMIPP_EQUAL_ACCURACY);
    Vhigh().pyramidReduce(Vlow(), levels);
    current_Ts *= pow(2.0, levels);
    Vhigh.clear();

    // Now scale using Bsplines
    int new_output_dim = CEIL(XSIZE(Vlow()) * current_Ts / Ts);
    Vlow().scaleToSizeBSpline(3, new_output_dim, new_output_dim, new_output_dim,
                                 Vhigh());
    Vlow() = Vhigh();
    Vlow().setXmippOrigin();

    // Return to the desired size
    Vlow().window(FIRST_XMIPP_INDEX(output_dim), FIRST_XMIPP_INDEX(output_dim),
                  FIRST_XMIPP_INDEX(output_dim), LAST_XMIPP_INDEX(output_dim),
                  LAST_XMIPP_INDEX(output_dim), LAST_XMIPP_INDEX(output_dim));
}

/* Blob properties --------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::blob_properties() const
{
    std::ofstream fh_out;
    fh_out.open((fn_out + "_Fourier_profile.txt").c_str());
    if (!fh_out)
        REPORT_ERROR(1, (string)"Prog_PDBPhantom_Parameters::blob_properties:"
                     " Cannot open " + fn_out + "_Fourier_profile.txt for output");
    fh_out << "# Freq(1/A) 10*log10(|Blob(f)|^2) Ts=" << highTs << endl;
    for (double w = 0; w < 1.0 / (2*highTs); w += 1.0 / (highTs * 500))
    {
        double H = kaiser_Fourier_value(w * highTs, periodicTable(0, 0) / highTs,
                                        blob.alpha, blob.order);
        fh_out << w << " " << 10*log10(H*H) << endl;
    }
    fh_out.close();
}

/* Atom descriptors -------------------------------------------------------- */
void atomDescriptors(const std::string &atom, Matrix1D<double> &descriptors) {
    descriptors.initZeros(11);
    if (atom=="H") {
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
    } else if (atom=="C") {
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
    } else if (atom=="N") {
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
    } else if (atom=="O") {
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
    } else if (atom=="P") {
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
    } else if (atom=="S") {
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
    } else if (atom=="Fe") {
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
    } else
        REPORT_ERROR(1,(std::string)"atomDescriptors: Unknown atom "+atom);
}

/* Electron form factor in Fourier ----------------------------------------- */
double electronFormFactorFourier(double f, const Matrix1D<double> &descriptors) {
    double retval=0;
    for (int i=1; i<=5; i++) {
        double ai=descriptors(i);
	double bi=descriptors(i+5);
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
    const Matrix1D<double> &descriptors) {
    double retval=0;
    for (int i=1; i<=5; i++) {
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
void hlpf(Matrix1D<double> &f, int M, double T, const string &filterType,
   Matrix1D<double> &filter, double reductionFactor=0.8,
   double ripple=0.01, double deltaw=1.0/8.0) {
   filter.initZeros(XSIZE(f));
   filter.setXmippOrigin();
   
   int Nmax=(int)CEIL(M/2.0);
   if (filterType=="SimpleAveraging") {
      FOR_ALL_ELEMENTS_IN_MATRIX1D(filter)
        if (ABS(i)<=Nmax) filter(i)=1.0/(2*Nmax+1);
   } else if (filterType=="SincKaiser") {
      SincKaiserMask(filter,reductionFactor*PI/M,ripple,deltaw);
      filter/=filter.sum();
      if (FINISHINGX(f)>FINISHINGX(filter))
          filter.window(STARTINGX(f),FINISHINGX(f));
      else
          f.window(STARTINGX(filter),FINISHINGX(filter));
   }
}

/* Convolution between f and the hlpf -------------------------------------- */
void fhlpf(const Matrix1D<double> &f, const Matrix1D<double> &filter,
   int M, Matrix1D<double> &convolution) {
   // Expand the two input signals
   int Nmax=FINISHINGX(filter);
   Matrix1D<double> auxF, auxFilter;
   auxF=f;
   auxFilter=filter;
   auxF.window(STARTINGX(f)-Nmax,FINISHINGX(f)+Nmax);
   auxFilter.window(STARTINGX(filter)-Nmax,FINISHINGX(filter)+Nmax);
   
   // Convolve in Fourier
   Matrix1D< complex<double> > F, Filter;
   FourierTransform(auxF,F);
   FourierTransform(auxFilter,Filter);
   F*=Filter;
   
   // Correction for the double phase factor and the double
   // amplitude factor
   FOR_ALL_ELEMENTS_IN_MATRIX1D(F) {
      double w; FFT_IDX2DIGFREQ(i,XSIZE(F),w);
      w*=2*PI;
      F(i)*=complex<double>(cos(w*(STARTINGX(auxFilter)-1)),
                            sin(w*(STARTINGX(auxFilter)-1)));
      F(i)*=XSIZE(auxFilter);
   }
   InverseFourierTransform(F,convolution);
   convolution.setXmippOrigin();
}

/* Optimization of the low pass filter to fit a given atom ----------------- */
Matrix1D<double> globalHlpfPrm(3);
Matrix1D<double> globalf;
int globalM;
double globalT;
string globalAtom;

double Hlpf_fitness(double *p) {
    double reductionFactor=p[1];
    double ripple=p[2];
    double deltaw=p[3];
    
    if (reductionFactor<0.7 || reductionFactor>1.3) return 1e38;
    if (ripple<0 || ripple>0.2) return 1e38;
    if (deltaw<0 || deltaw>0.2) return 1e38;
    
    // Construct the filter with the current parameters
    Matrix1D<double> filter, auxf;
    auxf=globalf;
    hlpf(auxf, globalM, globalT, "SincKaiser", filter, reductionFactor,
       ripple, deltaw);
    
    // Convolve the filter with the atomic profile
    Matrix1D<double> fhlpfFinelySampled;
    fhlpf(auxf, filter, globalM, fhlpfFinelySampled);

    // Coarsely sample
    double Rmax=FINISHINGX(fhlpfFinelySampled)*globalT;
    int imax=CEIL(Rmax/(globalM*globalT));
    Matrix1D<double> fhlpfCoarselySampled(2*imax+1);
    Matrix1D<double> splineCoeffsfhlpfFinelySampled;
    fhlpfFinelySampled.produceSplineCoefficients(
       splineCoeffsfhlpfFinelySampled,3);
    fhlpfCoarselySampled.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX1D(fhlpfCoarselySampled) {
       double r=i*(globalM*globalT)/globalT;
       fhlpfCoarselySampled(i)=
          splineCoeffsfhlpfFinelySampled.interpolatedElementBSpline(r,3);
    }

    // Build the frequency response of the convolved and coarsely sampled
    // atom
    Matrix1D<double> aux, FfilterMag, freq;
    Matrix1D< complex<double> > Ffilter;
    aux=fhlpfCoarselySampled;
    aux.window(-10*FINISHINGX(aux),10*FINISHINGX(aux));
    FourierTransform(aux,Ffilter);
    FFT_magnitude(Ffilter,FfilterMag);
    freq.initZeros(XSIZE(Ffilter));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FfilterMag)
       FFT_IDX2DIGFREQ(i,XSIZE(FfilterMag),freq(i));
    freq/=globalM*globalT;
    double amplitudeFactor=fhlpfFinelySampled.sum()/
                           fhlpfCoarselySampled.sum();

    // Compute the error in representation
    double error=0;
    Matrix1D<double> descriptors; atomDescriptors(globalAtom, descriptors);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FfilterMag)
        if (freq(i)>=0) {
            double f1=20*log10(FfilterMag(i)*XSIZE(FfilterMag)*amplitudeFactor);
	    double f2=20*log10(1/globalT*
	                 electronFormFactorFourier(freq(i),descriptors));
            error+=(f1-f2)*(f1-f2);
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
    of the cutoff frequency, bestPrm(1)=ripple of the Kaiser window,
    bestPrm(2)=deltaw of the Kaiser window.
*/
void optimizeHlpf(Matrix1D<double> &f, int M, double T, const string &atom,
    Matrix1D<double> &filter, Matrix1D<double> &bestPrm) {
    globalHlpfPrm(0)=1.0;     // reduction factor
    globalHlpfPrm(1)=0.01;    // ripple
    globalHlpfPrm(2)=1.0/8.0; // deltaw
    globalf=f;
    globalM=M;
    globalT=T;
    globalAtom=atom;
    double fitness;
    int iter;
    Matrix1D<double> steps(3); steps.init_constant(1);
    powellOptimizer(globalHlpfPrm, 1, 3,
                      &Hlpf_fitness, 0.05, fitness, iter, steps, false);
    bestPrm=globalHlpfPrm;
    hlpf(f, M, T, "SincKaiser", filter, bestPrm(0), bestPrm(1), bestPrm(2));
}

/* Atom radial profile ----------------------------------------------------- */
void atomRadialProfile(int M, double T, const string &atom,
    Matrix1D<double> &profile) {
    // Compute the electron form factor in real space
    double largestb1=76.7309/(4*PI*PI);
    double Rmax=4*sqrt(2*largestb1);
    int imax=(int)CEIL(Rmax/T);
    Matrix1D<double> descriptors;  atomDescriptors(atom, descriptors);
    Matrix1D<double> f(2*imax+1);
    f.setXmippOrigin();
    for (int i=-imax; i<=imax; i++) {
	double r=i*T;
	f(i)=electronFormFactorRealSpace(r,descriptors);
    }

    // Compute the optimal filter
    Matrix1D<double> filter, bestPrm;
    optimizeHlpf(f, M, T, atom, filter, bestPrm);
    
    // Perform the convolution
    fhlpf(f, filter, M, profile);
}

/* Create protein using scattering profiles -------------------------------- */
void Prog_PDBPhantom_Parameters::create_protein_using_scattering_profiles() {
    // Create an empty volume to hold the protein
    Vlow().initZeros(output_dim,output_dim,output_dim);
    Vlow().setXmippOrigin();

    // Fill the volume with the different atoms
    std::ifstream fh_pdb;
    fh_pdb.open(fn_pdb.c_str());
    if (!fh_pdb)
        REPORT_ERROR(1, (string)"Prog_PDBPhantom_Parameters::create_protein_using_scattering_profiles:"
                     "Cannot open " + fn_pdb + " for reading");

    // Process all lines of the file
    while (!fh_pdb.eof())
    {
        // Read an ATOM line
        string line;
        getline(fh_pdb, line);
        if (line == "") continue;
        string kind = firstToken(line);
        if (kind != "ATOM") continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        string dummy = nextToken();
        string atom_type = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        double x = textToFloat(nextToken());
        double y = textToFloat(nextToken());
        double z = textToFloat(nextToken());

        // Correct position
        Matrix1D<double> r(3);
        VECTOR_R3(r, x, y, z);
        r -= centerOfMass;
        r /= Ts;

        // Characterize atom
	int idx;
	switch (atom_type[0]) {
	   case 'H': idx=0; break;
	   case 'C': idx=1; break;
	   case 'N': idx=2; break;
	   case 'O': idx=3; break;
	   case 'P': idx=4; break;
	   case 'S': idx=5; break;
	   case 'F': idx=6; break;
	   default: continue;
	}
	double radius=(double)FINISHINGX(*(atomProfiles[idx]))/M;

        // Find the part of the volume that must be updated
        int k0 = MAX(FLOOR(ZZ(r) - radius), STARTINGZ(Vlow()));
        int kF = MIN(CEIL(ZZ(r) + radius), FINISHINGZ(Vlow()));
        int i0 = MAX(FLOOR(YY(r) - radius), STARTINGY(Vlow()));
        int iF = MIN(CEIL(YY(r) + radius), FINISHINGY(Vlow()));
        int j0 = MAX(FLOOR(XX(r) - radius), STARTINGX(Vlow()));
        int jF = MIN(CEIL(XX(r) + radius), FINISHINGX(Vlow()));

        // Fill the volume with this atom
        for (int k = k0; k <= kF; k++)
            for (int i = i0; i <= iF; i++)
                for (int j = j0; j <= jF; j++)
                {
                    Matrix1D<double> rdiff(3);
                    VECTOR_R3(rdiff, XX(r) - j, YY(r) - i, ZZ(r) - k);
		    double rdiffModule=rdiff.module();
		    if (rdiffModule<radius) 
                	Vlow(k, i, j) += (*(atomProfiles[idx])).
			   interpolatedElementBSpline(rdiffModule*M,3);
                }
    }

    // Close file
    fh_pdb.close();
}

/* Run --------------------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::run()
{
    produceSideInfo();
    compute_protein_geometry();
    std::cout << "Center of mass: " << centerOfMass.transpose() << std::endl
              << "Limits: " << limit.transpose() << std::endl;
    if (useBlobs) {
	create_protein_at_high_sampling_rate();
	create_protein_at_low_sampling_rate();
	blob_properties();
    } else {
    	create_protein_using_scattering_profiles();
    }
    Vlow.write(fn_out + ".vol");
}
