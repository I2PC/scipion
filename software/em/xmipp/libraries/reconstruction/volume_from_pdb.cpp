/***************************************************************************
 *
 * Authors:
 *
 * Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "volume_from_pdb.h"

#include <data/args.h>

#include <fstream>

/* Empty constructor ------------------------------------------------------- */
ProgPdbConverter::ProgPdbConverter()
{
    blob.radius = 2;   // Blob radius in voxels
    blob.order  = 2;   // Order of the Bessel function
    blob.alpha  = 3.6; // Smoothness parameter
    output_dim = -1;
    fn_pdb = "";
    Ts = 1;
    highTs = 1.0/12.0;
    useBlobs=false;
    usePoorGaussian=false;
    useFixedGaussian=false;
    doCenter=false;

    // Periodic table for the blobs
    periodicTable.resize(7, 2);
    periodicTable(0, 0) = atomRadius("H");
    periodicTable(0, 1) = atomCharge("H");
    periodicTable(1, 0) = atomRadius("C");
    periodicTable(1, 1) = atomCharge("C");
    periodicTable(2, 0) = atomRadius("N");
    periodicTable(2, 1) = atomCharge("N");
    periodicTable(3, 0) = atomRadius("O");
    periodicTable(3, 1) = atomCharge("O");
    periodicTable(4, 0) = atomRadius("P");
    periodicTable(4, 1) = atomCharge("P");
    periodicTable(5, 0) = atomRadius("S");
    periodicTable(5, 1) = atomCharge("S");
    periodicTable(6, 0) = atomRadius("Fe");
    periodicTable(6, 1) = atomCharge("Fe");

    // Correct the atom weights by the blob weight
    for (size_t i = 0; i < MAT_YSIZE(periodicTable); i++)
    {
        periodicTable(i, 1) /=
            basvolume(periodicTable(i, 0) / highTs, blob.alpha, blob.order, 3);
    }
}

/* Produce Side Info ------------------------------------------------------- */
void ProgPdbConverter::produceSideInfo()
{
    if (useFixedGaussian && sigmaGaussian<0)
    {
        // Check if it is a pseudodensity volume
        std::ifstream fh_pdb;
        fh_pdb.open(fn_pdb.c_str());
        if (!fh_pdb)
            REPORT_ERROR(ERR_IO_NOTEXIST, fn_pdb);
        while (!fh_pdb.eof())
        {
            // Read an ATOM line
            std::string line;
            getline(fh_pdb, line);
            if (line == "")
                continue;
            std::string kind = line.substr(0,6);
            if (kind!="REMARK")
                continue;
            std::vector< std::string > results;
            splitString(line," ",results);
            if (results[1]=="xmipp_convert_vol2pseudo")
                useFixedGaussian=true;
            if (useFixedGaussian && results[1]=="fixedGaussian")
                sigmaGaussian=textToFloat(results[2]);
            if (useFixedGaussian && results[1]=="intensityColumn")
                intensityColumn=results[2];
        }
        fh_pdb.close();
    }

    if (!useBlobs && !usePoorGaussian && !useFixedGaussian)
    {
        // Compute the downsampling factor
        M=(int)ROUND(Ts/highTs);

        // Atom profiles for the electron scattering method
        atomProfiles.setup(M,highTs,false);
    }
}

/* Atom description ------------------------------------------------------- */
void ProgPdbConverter::atomBlobDescription(
    const std::string &_element, double &weight, double &radius) const
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
    	if (verbose>0)
    		std::cout << "Unknown :" << _element << std::endl;
        return;
    }
    radius = periodicTable(idx, 0);
    weight = periodicTable(idx, 1);
}
/* Usage ------------------------------------------------------------------- */
void ProgPdbConverter::defineParams()
{
    addUsageLine("Covert a PDB file to a volume.");
    addExampleLine("Sample at 1.6A and limit the frequency to 10A",false);
    addExampleLine("   xmipp_volume_from_pdb -i 1o7d.pdb --sampling 1.6");
    addExampleLine("   xmipp_transform_filter -i 1o7d.vol -o 1o7dFiltered.vol --fourier low_pass 10 raised_cosine 0.05 --sampling 1.6");

    addParamsLine("   -i <pdb_file>                     : File to process");
    addParamsLine("  [-o <fn_root>]                     : Root name for output");
    addParamsLine("  [--sampling <Ts=1>]                : Sampling rate (Angstroms/pixel)");
    addParamsLine("  [--high_sampling_rate <highTs=0.08333333>]: Sampling rate before downsampling");
    addParamsLine("  [--size <output_dim=-1>]               : Final size in pixels (must be a power of 2, if blobs are used)");
    addParamsLine("  [--centerPDB]                       : Center PDB with the center of mass");
    addParamsLine("  [--blobs]                           : Use blobs instead of scattering factors");
    addParamsLine("  [--poor_Gaussian]                   : Use a simple Gaussian adapted to each atom");
    addParamsLine("  [--fixed_Gaussian <std=-1>]         : Use a fixed Gausian for each atom with");
    addParamsLine("                                     :  this standard deviation");
    addParamsLine("                                     :  If not given, the standard deviation is taken from the PDB file");
    addParamsLine("  [--intensityColumn <intensity_type=occupancy>]   : Where to write the intensity in the PDB file");
    addParamsLine("     where <intensity_type> occupancy Bfactor     : Valid values: occupancy, Bfactor");
}
/* Read parameters --------------------------------------------------------- */
void ProgPdbConverter::readParams()
{
    fn_pdb = getParam("-i");
    fn_out = checkParam("-o") ? getParam("-o") : fn_pdb.withoutExtension();
    Ts = getDoubleParam("--sampling");
    highTs = getDoubleParam("--high_sampling_rate");
    output_dim = getIntParam("--size");
    useBlobs = checkParam("--blobs");
    usePoorGaussian = checkParam("--poor_Gaussian");
    useFixedGaussian = checkParam("--fixed_Gaussian");
    if (useFixedGaussian)
        sigmaGaussian = getDoubleParam("--fixed_Gaussian");
    doCenter = checkParam("--centerPDB");
    intensityColumn = getParam("--intensityColumn");
}

/* Show -------------------------------------------------------------------- */
void ProgPdbConverter::show()
{
    if (verbose==0)
        return;
    std::cout << "PDB file:           " << fn_pdb           << std::endl
    << "Sampling rate:      " << Ts               << std::endl
    << "High sampling rate: " << highTs           << std::endl
    << "Size:               " << output_dim       << std::endl
    << "Center PDB:         " << doCenter         << std::endl
    << "Use blobs:          " << useBlobs         << std::endl
    << "Use poor Gaussian:  " << usePoorGaussian  << std::endl
    << "Use fixed Gaussian: " << useFixedGaussian << std::endl
    ;
    if (useFixedGaussian)
        std::cout << "Intensity Col:      " << intensityColumn  << std::endl
        << "Sigma:              " << sigmaGaussian  << std::endl;
}

/* Compute protein geometry ------------------------------------------------ */
void ProgPdbConverter::computeProteinGeometry()
{
    Matrix1D<double> limit0(3), limitF(3);
    computePDBgeometry(fn_pdb, centerOfMass, limit0, limitF, intensityColumn);
    if (doCenter)
    {
        limit0-=centerOfMass;
        limitF-=centerOfMass;
    }
    limit.resize(3);
    XX(limit) = XMIPP_MAX(ABS(XX(limit0)), ABS(XX(limitF)));
    YY(limit) = XMIPP_MAX(ABS(YY(limit0)), ABS(YY(limitF)));
    ZZ(limit) = XMIPP_MAX(ABS(ZZ(limit0)), ABS(ZZ(limitF)));

    // Update output size if necessary
    if (output_dim == -1)
    {
        int max_dim = XMIPP_MAX(CEIL(ZZ(limit) * 2 / Ts) + 5, CEIL(YY(limit) * 2 / Ts) + 5);
        max_dim = XMIPP_MAX(max_dim, CEIL(XX(limit) * 2 / Ts) + 5);
        if (useBlobs)
            output_dim = (int)NEXT_POWER_OF_2(max_dim);
        else
            output_dim = max_dim+10;
    }
}

/* Create protein at a high sampling rate ---------------------------------- */
void ProgPdbConverter::createProteinAtHighSamplingRate()
{
    // Create an empty volume to hold the protein
    int finalDim;
    if (highTs!=Ts)
        finalDim=(int)NEXT_POWER_OF_2(output_dim / (highTs/Ts));
    else
        finalDim=output_dim;
    Vhigh().initZeros(finalDim,finalDim,finalDim);
    Vhigh().setXmippOrigin();
    if (verbose)
    	std::cout << "The highly sampled volume is of size " << XSIZE(Vhigh())
    	<< std::endl;

    // Fill the volume with the different atoms
    std::ifstream fh_pdb;
    fh_pdb.open(fn_pdb.c_str());
    if (!fh_pdb)
        REPORT_ERROR(ERR_IO_NOTEXIST, fn_pdb);

    // Process all lines of the file
    int col=1;
    if (intensityColumn=="Bfactor")
        col=2;
    while (!fh_pdb.eof())
    {
        // Read an ATOM line
        std::string line;
        getline(fh_pdb, line);
        if (line == "")
            continue;
        std::string kind = line.substr(0,4);
        if (kind != "ATOM" && kind !="HETA")
            continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        std::string atom_type = line.substr(13,2);
        double x = textToFloat(line.substr(30,8));
        double y = textToFloat(line.substr(38,8));
        double z = textToFloat(line.substr(46,8));

        // Correct position
        Matrix1D<double> r(3);
        VECTOR_R3(r, x, y, z);
        if (doCenter)
            r -= centerOfMass;
        r /= highTs;

        // Characterize atom
        double weight, radius;
        if (!useFixedGaussian)
        {
            if (atom_type=="HETA")
                continue;
            atomBlobDescription(atom_type, weight, radius);
        }
        else
        {
            radius=4.5*sigmaGaussian;
            if (col==1)
                weight=textToFloat(line.substr(54,6));
            else
                weight=textToFloat(line.substr(60,6));
        }
        blob.radius = radius;
        if (usePoorGaussian)
            radius=XMIPP_MAX(radius/Ts,4.5);
        double GaussianSigma2=(radius/(3*sqrt(2.0)));
        if (useFixedGaussian)
            GaussianSigma2=sigmaGaussian;
        GaussianSigma2*=GaussianSigma2;
        double GaussianNormalization = 1.0/pow(2*PI*GaussianSigma2,1.5);

        // Find the part of the volume that must be updated
        int k0 = XMIPP_MAX(FLOOR(ZZ(r) - radius), STARTINGZ(Vhigh()));
        int kF = XMIPP_MIN(CEIL(ZZ(r) + radius), FINISHINGZ(Vhigh()));
        int i0 = XMIPP_MAX(FLOOR(YY(r) - radius), STARTINGY(Vhigh()));
        int iF = XMIPP_MIN(CEIL(YY(r) + radius), FINISHINGY(Vhigh()));
        int j0 = XMIPP_MAX(FLOOR(XX(r) - radius), STARTINGX(Vhigh()));
        int jF = XMIPP_MIN(CEIL(XX(r) + radius), FINISHINGX(Vhigh()));

        // Fill the volume with this atom
        Matrix1D<double> rdiff(3);
        for (int k = k0; k <= kF; k++)
            for (int i = i0; i <= iF; i++)
                for (int j = j0; j <= jF; j++)
                {
                    VECTOR_R3(rdiff, XX(r) - j, YY(r) - i, ZZ(r) - k);
                    rdiff*=highTs;
                    if (useBlobs)
                        Vhigh(k, i, j) += weight * blob_val(rdiff.module(), blob);
                    else if (usePoorGaussian || useFixedGaussian)
                        Vhigh(k, i, j) += weight *
                                          exp(-rdiff.module()*rdiff.module()/(2*GaussianSigma2))*
                                          GaussianNormalization;
                }
    }

    // Close file
    fh_pdb.close();
}

/* Create protein at a low sampling rate ----------------------------------- */
void ProgPdbConverter::createProteinAtLowSamplingRate()
{
    // Compute the integer downsapling factor
    int M = FLOOR(Ts / highTs);
    double current_Ts = highTs;

    // Use Bsplines pyramid if possible
    int levels = FLOOR(log10((double)M) / log10(2.0) + XMIPP_EQUAL_ACCURACY);
    pyramidReduce(BSPLINE3, Vlow(), Vhigh(), levels);
    current_Ts *= pow(2.0, levels);
    Vhigh.clear();

    // Now scale using Bsplines
    int new_output_dim = CEIL(XSIZE(Vlow()) * current_Ts / Ts);
    scaleToSize(BSPLINE3, Vhigh(), Vlow(),
                new_output_dim, new_output_dim, new_output_dim);
    Vlow() = Vhigh();
    Vlow().setXmippOrigin();

    // Return to the desired size
    Vlow().selfWindow(FIRST_XMIPP_INDEX(output_dim), FIRST_XMIPP_INDEX(output_dim),
                  FIRST_XMIPP_INDEX(output_dim), LAST_XMIPP_INDEX(output_dim),
                  LAST_XMIPP_INDEX(output_dim), LAST_XMIPP_INDEX(output_dim));
}

/* Blob properties --------------------------------------------------------- */
void ProgPdbConverter::blobProperties() const
{
    std::ofstream fh_out;
    fh_out.open((fn_out + "_Fourier_profile.txt").c_str());
    if (!fh_out)
        REPORT_ERROR(ERR_IO_NOWRITE, fn_out);
    fh_out << "# Freq(1/A) 10*log10(|Blob(f)|^2) Ts=" << highTs << std::endl;
    for (double w = 0; w < 1.0 / (2*highTs); w += 1.0 / (highTs * 500))
    {
        double H = kaiser_Fourier_value(w * highTs, periodicTable(0, 0) / highTs,
                                        blob.alpha, blob.order);
        fh_out << w << " " << 10*log10(H*H) << std::endl;
    }
    fh_out.close();
}

/* Create protein using scattering profiles -------------------------------- */
void ProgPdbConverter::createProteinUsingScatteringProfiles()
{
    // Create an empty volume to hold the protein
    Vlow().initZeros(output_dim,output_dim,output_dim);
    Vlow().setXmippOrigin();

    // Fill the volume with the different atoms
    std::ifstream fh_pdb;
    fh_pdb.open(fn_pdb.c_str());
    if (!fh_pdb)
        REPORT_ERROR(ERR_IO_NOTEXIST, fn_pdb);

    // Process all lines of the file
    std::string line, kind, atom_type;
    double iTs=1.0/Ts;
    Matrix1D<double> r(3), rdiff(3);
    while (!fh_pdb.eof())
    {
        // Read an ATOM line
        getline(fh_pdb, line);
        if (line == "")
            continue;
        kind =line.substr(0,4);
        if (kind != "ATOM")
            continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        atom_type = line.substr(13,2);
        char atom_type0=atom_type[0];
        double x = textToFloat(line.substr(30,8));
        double y = textToFloat(line.substr(38,8));
        double z = textToFloat(line.substr(46,8));

        // Correct position
        VECTOR_R3(r, x, y, z);
        if (doCenter)
            r -= centerOfMass;
        r *= iTs;

        // Characterize atom
        try
        {
            double radius=atomProfiles.atomRadius(atom_type[0]);
            double radius2=radius*radius;

            // Find the part of the volume that must be updated
            const MultidimArray<double> &mVlow=Vlow();
            int k0 = XMIPP_MAX(FLOOR(ZZ(r) - radius), STARTINGZ(mVlow));
            int kF = XMIPP_MIN(CEIL(ZZ(r) + radius), FINISHINGZ(mVlow));
            int i0 = XMIPP_MAX(FLOOR(YY(r) - radius), STARTINGY(mVlow));
            int iF = XMIPP_MIN(CEIL(YY(r) + radius), FINISHINGY(mVlow));
            int j0 = XMIPP_MAX(FLOOR(XX(r) - radius), STARTINGX(mVlow));
            int jF = XMIPP_MIN(CEIL(XX(r) + radius), FINISHINGX(mVlow));

            // Fill the volume with this atom
            for (int k = k0; k <= kF; k++)
            {
                double zdiff=ZZ(r) - k;
                double zdiff2=zdiff*zdiff;
                for (int i = i0; i <= iF; i++)
                {
                    double ydiff=YY(r) - i;
                    double zydiff2=zdiff2+ydiff*ydiff;
                    for (int j = j0; j <= jF; j++)
                    {
                        double xdiff=XX(r) - j;
                        double rdiffModule2=zydiff2+xdiff*xdiff;
                        if (rdiffModule2<radius2)
                        {
                            double rdiffModule=sqrt(rdiffModule2);
                            A3D_ELEM(mVlow,k, i, j) += atomProfiles.volumeAtDistance(
                                                 atom_type0,rdiffModule);
                        }
                    }
                }
            }
        }
        catch (XmippError XE)
        {
        	if (verbose)
        		std::cerr << "Ignoring atom of type *" << atom_type << "*" << std::endl;
        }
    }

    // Close file
    fh_pdb.close();
}

/* Run --------------------------------------------------------------------- */
void ProgPdbConverter::run()
{
    produceSideInfo();
    show();
    computeProteinGeometry();
    if (useBlobs)
    {
        createProteinAtHighSamplingRate();
        createProteinAtLowSamplingRate();
        blobProperties();
    }
    else if (usePoorGaussian || useFixedGaussian)
    {
        highTs=Ts;
        createProteinAtHighSamplingRate();
        Vlow=Vhigh;
        Vhigh.clear();
    }
    else
    {
        createProteinUsingScatteringProfiles();
    }
    if (fn_out!="")
        Vlow.write(fn_out + ".vol");
}
