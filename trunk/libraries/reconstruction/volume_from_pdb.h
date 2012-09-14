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
#ifndef _PROG_VOLUME_FROM_PDB_HH
#  define _PROG_VOLUME_FROM_PDB_HH

#include <data/blobs.h>
#include <data/pdb.h>
#include <data/xmipp_program.h>

/**@defgroup PDBPhantom convert_pdb2vol (PDB Phantom program)
   @ingroup ReconsLibrary */
//@{
/* PDB Phantom Program Parameters ------------------------------------------ */
/** Parameter class for the PDB Phantom program */
class ProgPdbConverter: public XmippProgram
{
public:
    /** Sampling rate */
    double Ts;

    /** PDB file */
    FileName fn_pdb;

    /** Output fileroot */
    FileName fn_out;

    /** Blob */
    struct blobtype blob;

    /** Final size in pixels */
    int output_dim;
    
    /** Use blobs instead of scattering factors */
    bool useBlobs;

    /** Use poor Gaussian instead of scattering factors */
    bool usePoorGaussian;
    
    /** Use fixed Gaussian instead of scattering factors */
    bool useFixedGaussian;
    
    /** Center the PDB */
    bool doCenter;

    /** Fixed Gaussian standard deviation */
    double sigmaGaussian;

    /// Column for the intensity (if any). Only valid for fixed_gaussians
    std::string intensityColumn;
public:
    /** Empty constructor */
    ProgPdbConverter();

    /** Params definitions */
    void defineParams();
    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void readParams();

    /** Produce side information.
        Produce the atomic profiles. */
    void produceSideInfo();

    /** Show parameters. */
    void show();

    /** Run. */
    void run();
public:
    /* Downsampling factor */
    int M;

    /* Sampling rate used for highly sampled volumes */
    double highTs;

    /* periodic_table(i,0)=radius
       periodic_table(i,1)=atomic weight */
    Matrix2D<double> periodicTable;

    /* Atom interpolator. */
    AtomInterpolator atomProfiles;

    // Protein geometry
    Matrix1D<double> centerOfMass, limit;

    /* Volume at a high sampling rate */
    Image<double> Vhigh;

    /* Volume at a low sampling rate */
    Image<double> Vlow;

    /* Blob properties at the high sampling rate */
    void blobProperties() const;

    /* Atom weight and radius */
    void atomBlobDescription(const std::string &_element,
        double &weight, double &radius) const;

    /* Protein geometry */
    void computeProteinGeometry();

    /* Create protein at a high sampling rate */
    void createProteinAtHighSamplingRate();

    /* Create protein at a low sampling rate */
    void createProteinAtLowSamplingRate();

    /* Create protein using scattering profiles */
    void createProteinUsingScatteringProfiles();
};
//@}
#endif
