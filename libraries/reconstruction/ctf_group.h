/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#ifndef _CTF_GROUP_HH
#define _CTF_GROUP_HH

#include <data/args.h>
#include <data/filters.h>
#include "fourier_filter.h"
#include <data/ctf.h>

#include <data/xmipp_program.h>

/**@defgroup CTFGroup ctf_group
   @ingroup ReconsLibrary */
//@{
/// CTFGroup program class
class ProgCtfGroup: public XmippProgram
{
public:
    /// Filenames
    //FileName fn_sel,
    FileName fn_ctfdat, fn_root, fn_split,format;

    /// Maximum allowed error
    double max_error;

    /// Resolution for maximum allowed error (in dig. freq.)
    double resol_error;

    /// Resolution for maximum allowed error (in pixels)
    int iresol_error;

    // Pixel size  (Angstrom)
    double pixel_size;

    /// Flag for phase-flipped data
    bool phase_flipped;

    // Size of the images and of the output CTFs
    size_t dim, xpaddim, ypaddim,paddim,ctfxpaddim;

    // Padding factor
    double pad;

    // Flag whether to use automated
    bool do_auto;

    // Flag whether to throw out anisotropic CTFs
    bool do_discard_anisotropy;

    /// Flag for calculating Wiener filters
    bool do_wiener;

    /// Wiener filter constant
    double wiener_constant;

    /// Replace ctf.param file sampling rate by this
    bool replaceSampling;

    /// New ctf sampling rate
    double samplingRate;

    /// Compute 1D CTF using avg(defocusU + defocusV)/2 as defocus.
    /// This approach speed up the computation and is recommended for very larg data sets
    bool do1Dctf;

    /// Available memory, set mmap to on if more memory is needed
    double memory;// (Gigabits)

    ///do not use ram but map data to a file
    bool mmapOn;

    /// Matrix with denominator term of Wiener filter
    MultidimArray<double> Mwien;

    ///auxiliary matrices to speed up process
    MultidimArray<double> diff;
    MultidimArray<int> dd;

    // Array with  CTF profiles for all micrographs and all groups
    //std::vector<MultidimArray<double> > mics_ctf2d, groups_ctf2d;
    MultidimArray<double> mics_ctf2d;//, groups_ctf2d;

    // Vector with average defocus value per micrograph and per group
    std::vector<double> mics_defocus;

    //Metadata with ctfs sorted by average defocus
    MetaData sortedCtfMD;

    // Vector with number of images per defocus group
    //std::vector<int> mics_count;

    // Vector of vector with image names for all micrographs
    //std::vector< std::vector <FileName> > mics_fnimgs;

    // Pointers
    std::vector< std::vector <int> > pointer_group2mic;

public:

    /** Read parameters from command line. */
    void readParams();

    /** Define parameters. */
    void defineParams();

    /** Show. */
    void show();

    /** Usage. */
    void usage();

    /** Check anisotropy of a single CTF */
    bool isIsotropic(CTFDescription &ctf);

    /** Produce side information. */
    void produceSideInfo();

    /** Produce the CTF groups automatically */
    void autoRun();

    /** Split based onto defocus values  given in a docfile */
    void manualRun();

    /** Write output */
    void writeOutputToDisc();

    /** Do the job */
    void run();

};
//@}
#endif
