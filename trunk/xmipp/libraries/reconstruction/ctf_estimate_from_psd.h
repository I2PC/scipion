/***************************************************************************
 *
 * Authors:      Javier Angel Velazquez Muriel    javi@cnb.csic.es
 *               Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef _ADJUST_CTF_HH
#define _ADJUST_CTF_HH

#include "fourier_filter.h"

/**@defgroup AdjustParametricCTF adjust_ctf (Adjust CTF parameters to PSD)
   @ingroup ReconsLibrary */
//@{
/** Adjust CTF parameters. */
class ProgCTFEstimateFromPSD: public XmippProgram
{
public:
    /// CTF filename
    FileName             fn_psd;
    /// Downsample performed
    double				 downsampleFactor;
    /// CTF amplitude to model
    Image<double>        ctftomodel;
    /// CTF amplitude to model
    Image<double>        enhanced_ctftomodel;
    /// CTF amplitude to model
    Image<double>        enhanced_ctftomodel_fullsize;
    /// CTF model
    CTFDescription       initial_ctfmodel;
    /// Show convergence values
    bool                 show_optimization;
    /// X dimension of particle projections (-1=the same as the psd)
    int                  ctfmodelSize;
    /// Bootstrap estimation
    bool                 bootstrap;
    /// Fast defocus estimate
    bool                 fastDefocusEstimate;
    /// Regularization factor for the phase direction and unwrapping estimates (used in Zernike estimate)
    double               lambdaPhase;
    /// Size of the average window used during phase direction and unwrapping estimates (used in Zernike estimate)
    int                  sizeWindowPhase;
    /// Minimum frequency to adjust
    double               min_freq;
    /// Maximum frequency to adjust
    double               max_freq;
    /// Sampling rate
    double               Tm;
    /// Defocus range
    double               defocus_range;

    /// Enhancement filter low freq
    double               f1;
    /// Enhancement filter high freq
    double               f2;
    /// Weight of the enhanced image
    double               enhanced_weight;

    /// Set of parameters for the complete adjustment of the CTF
    Matrix1D<double>     adjust;
    
    /// Model simplification
    int                  modelSimplification;
public:
    /// Read parameters
    void readParams();

    /// Read parameters
    void readBasicParams(XmippProgram *program);

    /// Show parameters
    void show();

    /// Define basic parameters
    static void defineBasicParams(XmippProgram * program);

    /// Define Parameters
    void defineParams();

    /// Produce side information
    void produce_side_info();

    /** Generate half-plane model at a given size.
        It is assumed that ROUT_Adjust_CTF has been already run */
    void generate_model_halfplane(int Ydim, int Xdim, MultidimArray<double> &model);

    /** Generate quadrant model at a given size.
        It is assumed that ROUT_Adjust_CTF has been already run */
    void generate_model_quadrant(int Ydim, int Xdim, MultidimArray<double> &model);

    /** Run */
    void run();
};

/** Core of the Adjust CTF routine.
    This is the routine which does everything. It returns the fitting error
    committed in the best fit.*/
double ROUT_Adjust_CTF(ProgCTFEstimateFromPSD &prm, CTFDescription &output_ctfmodel, 
    bool standalone = true);
//@}
#endif
