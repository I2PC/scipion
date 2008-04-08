/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.uam.es)
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

#ifndef _WIENER_FILTER_HH
#define _WIENER_FILTER_HH

#include <data/args.h>
#include <data/selfile.h>
#include <data/filters.h>
#include "fourier_filter.h"
#include "ctf.h"

/**@defgroup WienerFilter ctf_wiener_filter 
   @ingroup ReconsLibraryPrograms */
//@{
/// Wiener Filter parameters
class WienerFilterParams
{
public:
    /// Filenames
    FileName fn_sel, fn_ctfdat, fn_root;

    /// Wiener filter constant
    double wiener_constant;

    /// Flag for phase-flipped data
    bool phase_flipped;

    // Size of the images
    int dim;

    // vector with all CTF filenames
    std::vector<FileName> all_fn_ctfs;

    // Vector with all ctf matrices
    std::vector<Matrix2D<double> > all_Mctfs;

    // Sum of all CTF^2 + wiener_constant
    Matrix2D<double> Msum;

    // Vector with number of images per defocus group
    std::vector<int> count_defocus;


public:

    /** Read parameters from command line. */
    void read(int argc, char **argv);

    /** Show. */
    void show();

    /** Usage. */
    void usage();

    /** Produce side information.
        The CTFdat is read. */
    void produceSideInfo();

    /** Produce the wiener filter images. */
    void run();

};
//@}
#endif
