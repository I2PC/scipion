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
    highTs = 0.1;

    periodic_table.resize(7, 2);
    periodic_table(0, 0) = 0.53;
    periodic_table(0, 1) = 1;  // Hydrogen
    periodic_table(1, 0) = 0.67;
    periodic_table(1, 1) = 12; // Carbon
    periodic_table(2, 0) = 0.56;
    periodic_table(2, 1) = 14; // Nitrogen
    periodic_table(3, 0) = 0.48;
    periodic_table(3, 1) = 16; // Oxygen
    periodic_table(4, 0) = 0.98;
    periodic_table(4, 1) = 31; // Phosphorus
    periodic_table(5, 0) = 0.88;
    periodic_table(5, 1) = 32; // Sulfur
    periodic_table(6, 0) = 1.56;
    periodic_table(6, 1) = 55.85; // Iron
    // Correct the atom weights by the blob weight
    for (int i = 0; i < YSIZE(periodic_table); i++)
    {
        periodic_table(i, 1) /=
            basvolume(periodic_table(i, 0) / highTs, blob.alpha, blob.order, 3);
    }
}

/* Atom description ------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::atom_description(const string &_element,
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
    default:
        std::cout << "Unknown :" << _element << std::endl;
        return;
    }
    radius = periodic_table(idx, 0);
    weight = periodic_table(idx, 1);
}

/* Read parameters --------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::read(int argc, char **argv)
{
    fn_pdb = getParameter(argc, argv, "-i");
    fn_out = getParameter(argc, argv, "-o", "");
    if (fn_out == "") fn_out = fn_pdb.without_extension();
    Ts = AtoF(getParameter(argc, argv, "-sampling_rate", "1"));
    highTs = AtoF(getParameter(argc, argv, "-high_sampling_rate", "0.1"));
    output_dim = AtoI(getParameter(argc, argv, "-output_dim", "-1"));
}

/* Usage ------------------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::usage()
{
    std::cerr << "convert_pdb2vol\n"
	      << "   -i <pdb file>		      : File to process\n"
	      << "  [-o <fn_root>]		      : Root name for output\n"
	      << "  [-sampling_rate <Ts=1>]	      : Sampling rate (Angstroms/pixel)\n"
	      << "  [-high_sampling_rate <highTs=0.1>]: Sampling rate before downsampling\n"
	      << "  [-size <output_dim>]	      : Final size in pixels (must be a power of 2)\n"
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
        double x = AtoF(nextToken());
        double y = AtoF(nextToken());
        double z = AtoF(nextToken());

        // Update center of mass and limits
        if (x < XX(limit0)) XX(limit0) = x;
        else if (x > XX(limitF)) XX(limitF) = x;
        if (y < YY(limit0)) YY(limit0) = y;
        else if (y > YY(limitF)) YY(limitF) = y;
        if (z < ZZ(limit0)) ZZ(limit0) = z;
        else if (z > ZZ(limitF)) ZZ(limitF) = z;
        double weight, radius;
        atom_description(atom_type, weight, radius);
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
        double x = AtoF(nextToken());
        double y = AtoF(nextToken());
        double z = AtoF(nextToken());

        // Correct position
        Matrix1D<double> r(3);
        VECTOR_R3(r, x, y, z);
        r -= centerOfMass;
        r /= highTs;

        // Characterize atom
        double weight, radius;
        atom_description(atom_type, weight, radius);
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

#ifdef NEVER_DEFINED
    // Filter Vlow before downsampling
    Vlow.write("PPPlow_before_filtering.vol");
    Matrix3D<complex<double> > FFTVlow;
    FourierTransform(Vlow(), FFTVlow);
    STARTINGX(FFTVlow) = STARTINGY(FFTVlow) = STARTINGZ(FFTVlow) = 0;
    Matrix1D<double> dig_freq(3);
    double freq_c = current_Ts / (2 * Ts);
    std::cout << "Filtering at " << freq_c << std::endl;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(FFTVlow)
    {
        FFT_IDX2DIGFREQ(i, XSIZE(FFTVlow), XX(dig_freq));
        FFT_IDX2DIGFREQ(j, YSIZE(FFTVlow), YY(dig_freq));
        FFT_IDX2DIGFREQ(k, ZSIZE(FFTVlow), ZZ(dig_freq));
        double freq_module = dig_freq.module();
        if (freq_module > freq_c) FFTVlow(k, i, j) = 0;
    }
    InverseFourierTransform(FFTVlow, Vlow());
    Vlow().setXmippOrigin();
    Vlow.write("PPPlow_after_filtering.vol");
#endif

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
        double H = kaiser_Fourier_value(w * highTs, periodic_table(0, 0) / highTs,
                                        blob.alpha, blob.order);
        fh_out << w << " " << 10*log10(H*H) << endl;
    }
    fh_out.close();
}

/* Run --------------------------------------------------------------------- */
void Prog_PDBPhantom_Parameters::run()
{
    compute_protein_geometry();
    std::cout << "Center of mass: " << centerOfMass.transpose() << std::endl
              << "Limits: " << limit.transpose() << std::endl;
    create_protein_at_high_sampling_rate();
    create_protein_at_low_sampling_rate();
    blob_properties();
    Vlow.write(fn_out + ".vol");
}
