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
#ifndef _PROG_SYMMETRIZE_HH
#  define _PROG_SYMMETRIZE_HH

#include <data/xmipp_funcs.h>
#include <data/xmipp_image.h>
#include <data/mask.h>
#include <data/symmetries.h>
#include <data/xmipp_program.h>

/**@defgroup SymmetrizeProgram symmetrize (Symmetrize a volume or image)
   @ingroup ReconsLibrary */
//@{
/* Test parameters --------------------------------------------------------- */
/// Symmetrize Parameters
class ProgSymmetrize : public XmippMetadataProgram
{
public:
    /// symmetry file
    FileName        fn_sym;
    /// Do not generate subgroup
    bool            do_not_generate_subgroup;
    /// wrap or don't wrap input file during symmetrisation
    bool            wrap;
    /// Sum or average the result
    bool            sum;
public:
    /** Read parameters from command line. */
    void readParams();

    /** Define Parameters */
    void defineParams();

    /** Show parameters */
    void show();

    /** Run */
    void preProcess();

    /// Process image or volume
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
public:
    // Symmetry description for volumes
    SymList SL;
    // Symmetry descriptio for images
    int symorder;
};

/** Symmetrize volume.*/
void symmetrizeVolume(const SymList &SL, const MultidimArray<double> &V_in,
                      MultidimArray<double> &V_out,
                      bool wrap=true, bool do_outside_avg=false, bool sum=false);

/** Symmetrize image.*/
void symmetrizeImage(int symorder, const MultidimArray<double> &I_in,
                      MultidimArray<double> &I_out,
                      bool wrap=true, bool do_outside_avg=false, bool sum=false);
//@}
#endif
