/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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

#ifndef MAKESPECTRA_H
#define MAKESPECTRA_H

#include "rotational_spectrum.h"
#include "progs.h"
#include "image.h"

#include <vector>

/// @defgroup MakeSpectra Make spectra.

/** Make spectra parameters.
 * @ingroup MakeSpectra
 */
class Prog_make_spectra_prm: public Prog_parameters
{
public:
    /** Output filename.
     */
    FileName fn_out;

    /** Rotational spectrum.
     */
    Rotational_Spectrum rot_spt;

public:
    /** Set of harmonics.
     */
    std::vector< matrix1D< double > > Harmonics;

    /** Set of images.
     */
    std::vector< FileName > Img_name;

    /** Empty constructor.
     */
    Prog_make_spectra_prm();

    /** Read parameters from command line.
     */
    void read(int argc, char** argv);

    /** Produce side info.
     */
    void produce_side_info();

    /** Show parameters. This function calls show_specific.
     */
    void show();

    /** Show specific.
     */
    void show_specific();

    /** Usage. This function calls usage_specific.
     */
    void usage();

    /** Show specific parameters.
     */
    void usage_specific();

    /** Process image.
     */
    void process_img(ImageXmipp& img);

    /** Finish processing.
     */
    void finish_processing();
};

#endif
