/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2010)
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

#ifndef _TOMO_EXTRACT_SUBVOLUME_H
#define _TOMO_EXTRACT_SUBVOLUME_H

#include <data/xmipp_fftw.h>
#include <data/xmipp_fft.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/xmipp_image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/ctf.h>
#include <data/sampling.h>
#include <data/symmetries.h>
#include "symmetrize.h"
#include <pthread.h>
#include <vector>

/** tomo_extract_subvolume parameters. */
class ProgTomoExtractSubvolume: public XmippMetadataProgram
{
public:
    /** Filenames reference selfile/image, fraction metadata & output rootname */
    FileName fn_doc, fn_root, fn_sym, fn_missing, fn_mask, fn_aux;

    /** Size of the output subvolume */
    int size;

    /** Coordinates of the center of the subvolume in the reference */
    Matrix1D<double> center_ref;

    /** Coordinates of the unique symmetry related centers in the reference */
    std::vector<Matrix1D<double> > centers_subvolumes;
    std::vector<Matrix2D<double> > rotations_subvolumes;

    /** Number of input volumes */
    int nr_exp_images;

    /** Number of subvolumes per input volume */
    int nr_output_volumes;

    // Metadata with input volumes names and their alignment parameters
    //MetaData DF;

    // Docfile with output subvolumes names and their alignment parameters
    // This file can be fed directly to ml_tomo again...
    MetaData DFout;

    // Selfile with output subvolumes
    //MetaData SFout;

    /** Minimum distance between unique subvolumes */
    double mindist;

    // Symmetry setup
    int symmetry, sym_order;
    SymList SL;

    //Some local variables
    FileName fn_out;
    Image<double> vol, volout;
    Matrix1D<double> center, doccenter;
    Matrix2D<double> A, R, I;
    double rot, tilt, psi, rotp, tiltp, psip;
    int x0, xF;

public:

    /// Define the arguments accepted
    void defineParams();

    /// Read arguments from command line
    void readParams();

    /// Show
    void show();

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    void postProcess();
    void preProcess();


    //void processImage(int imgno_start, int imgno_end);

};
//@}
#endif
