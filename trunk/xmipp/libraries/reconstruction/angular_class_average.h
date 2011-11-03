/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.csic.es (2008)
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

#ifndef _ANGULAR_CLASS_AVERAGE_H
#define _ANGULAR_CLASS_AVERAGE_H

#include <data/xmipp_fftw.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/polar.h>

#include <data/xmipp_program.h>

#define AVG_OUPUT_SIZE 10

/**@defgroup ClassAverage Create class averages from projection
   matching docfiles
   @ingroup ReconsLibrary */
//@{
/** angular_class_average parameters. */
class ProgAngularClassAverage : public XmippProgram
{

public:

    /** Input and library docfiles */
    MetaData         DF, DFlib;
    /** metadata with classes which have experimental images applied to them */
    MetaData         DFclassesExp;
    /** Output rootnames */
    FileName         fn_out, fn_out1, fn_out2, fn_wien;
    /** Column numbers */
    std::string      col_select;
    /** Upper and lower absolute and relative selection limits */
    double           limit0, limitF, limitR;
    /** Flags wether to use limit0, limitF and limitR selection */
    bool             do_limit0, do_limitF, do_limitR0, do_limitRF;
    /** Flag whether to apply mirror operations. By default set to True */
    bool             do_mirrors;
    /** Flag whether also to write out class averages of random halves of the data */
    bool             do_split;
    /** Image dimensions before and after padding (only for Wiener correction) */
    int              dim, paddim;
    /** Padding factor */
    double           pad;
    /** One empty image with correct dimensions */
    Image<double>    Iempty;
    /** Do NOT skip writing of selfiles */
    bool             write_selfiles;
    /** Number of 3d references */
    int              number_3dref;
    /** Delete auxiliary files from previous execution.
     * Alloc disk space for output stacks */
    bool             do_preprocess;
    /** Create block with average images filenames */
    bool             do_postprocess;
    /** Add output to existing files */
    bool             do_add;
    /** Wiener filter image */
    MultidimArray<double> Mwien;
    /** Selfiles containing all class averages */
    MetaData         SFclasses, SFclasses1, SFclasses2;

    /** Re-alignment of classes */

    /** Input file */
    FileName inFile, refGroup;
    /** Inner and outer radius for rotational alignment */
    int Ri, Ro;
    /** Number of iterations */
    int nr_iter;
    /** Convergence criterion */
    double eps;
    /** Search shift (shifts larger than this will be set to 0)*/
    double max_shift;
    /** Maximum allowed shift in last iteration (shifts larger than this will be set to 0)*/
    double max_shift_change, max_psi_change;
    /** transformers for all rings */
    Polar_fftw_plans global_plans;
    FourierTransformer global_transformer;
    MultidimArray<double> corr;

public:

    /// Read arguments from command line
    void readParams();

    /// Read arguments from command line
    void defineParams();

    /** Run. */
    void run();

    /** Make shiftmask and calculate nr_psi */
    void produceSideInfo();

    /** Delete auxiliary files from previous execution.
     * Alloc disk space for output stacks */
    void preprocess();

    /** Read all the metadata blocks with images assigned to different projection directions
     * and creates a single metadata block with the name of the projection direction averages. */
    void postprocess();

    /** Convert from cartesian to FT of polar coordinates */
    void getPolar(MultidimArray<double> &img, Polar<std::complex <double> > &fP,
                  bool conjugated=false, float xoff = 0., float yoff = 0.);

    /** Re-align all images in a class */
    void reAlignClass(Image<double> &avg1,
                      Image<double> &avg2,
                      MetaData      &SFclass1,
                      MetaData      &SFclass2,
                      std::vector< Image<double> > imgs,
                      std::vector<int> splits,
                      std::vector<int> numbers,
                      size_t dirno,
                      double * my_output);

    /** Apply Wiener filter */
    void applyWienerFilter(MultidimArray<double> &img);

    /** Process a single class */
    void processOneClass(size_t &dirno,
                         double * my_output);

    /** Write average and selfiles to disc */
    void writeToDisc(Image<double> avg,
                     size_t     dirno,
                     MetaData   SF,
                     FileName   fn,
                     bool       write_selfile,
                     FileName   oext="stk");

    /** Add an image to the selfile and docfile of all class averages */
    void addClassAverage(size_t dirno,
                         double w,
                         double w1,
                         double w2);

    /** Write selfile and docfile of all class averages to disc */
    void finalWriteToDisc();
};
//@}
#endif
