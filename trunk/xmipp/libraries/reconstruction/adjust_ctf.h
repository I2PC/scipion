/***************************************************************************
 *
 * Authors:      Javier Angel Velazquez Muriel    javi@cnb.uam.es
 *               Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#ifndef _ADJUST_CTF_HH
#define _ADJUST_CTF_HH

#include "fourier_filter.h"

/**@name Adjust parametric CTF */
//@{
/** Adjust CTF parameters. */
class Adjust_CTF_Parameters
{
public:
    /// CTF filename
    FileName             fn_ctf;
    /// Model to which the current one must be similar
    FileName             fn_similar_model;
    /// CTF amplitude to model
    ImageXmipp           ctftomodel;
    /// CTF amplitude to model
    ImageXmipp           enhanced_ctftomodel;
    /// CTF amplitude to model
    ImageXmipp           enhanced_ctftomodel_fullsize;
    /// CTF model
    XmippCTF             initial_ctfmodel;
    /// Show convergence values
    bool                 show_optimization;
    /// Allow astigmatic noise
    bool                 astigmatic_noise;
    /// X dimension of particle projections (-1=the same as the psd)
    int   ctfmodelSize;

    /// Minimum frequency to adjust
    double               min_freq;
    /// Maximum frequency to adjust
    double               max_freq;
    /// Sampling rate
    double               Tm;
    /// Defocus range
    double               defocus_range;

    /// Initial Chromatic aberration
    double               initial_Ca;

    /// Enhancement filter low freq
    double               f1;
    /// Enhancement filter high freq
    double               f2;
    /// Weight of the enhanced image
    double               enhanced_weight;

    /// Set of parameters for the complete adjustment of the CTF
    matrix1D<double>     adjust;
public:
    /// Read parameters from file
    void read(const FileName &fn_param);

    /// Write to a file
    void write(const FileName &fn, bool rewrite = true);

    /// Show parameters
    void show();

    /// Usage
    void Usage();

    /// Produce side information
    void produce_side_info();

    /** Generate half-plane model at a given size.
        It is assumed that ROUT_Adjust_CTF has been already run */
    void generate_model_halfplane(int Ydim, int Xdim, matrix2D<double> &model);

    /** Generate quadrant model at a given size.
        It is assumed that ROUT_Adjust_CTF has been already run */
    void generate_model_quadrant(int Ydim, int Xdim, matrix2D<double> &model);
};

/** Core of the Adjust CTF routine.
    This is the routine which does everything. It returns the fitting error
    committed in the best fit.*/
double ROUT_Adjust_CTF(Adjust_CTF_Parameters &prm, bool standalone = true);
//@}
#endif
