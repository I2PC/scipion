/***************************************************************************
*
* Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
* Updated by:  J.M. de la Rosa Trevin  (jmdelarosa@cnb.csic.es)
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

#ifndef _MEAN_SHIFT_H
#define _MEAN_SHIFT_H

#include <string>
#include <data/filters.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_threads.h>
#include <data/histogram.h>

/// @defgroup Denoise Image denoising
/// @ingroup ReconsLibrary

/// Parameters for denoise program
/// @ingroup Denoise
class MeanShiftFilter: public XmippFilter
{
public:
    double sigma_r, sigma_s;
    int iters, numThreads;
    bool fast, save_iters;
    /** Verbosity used by denoiser.
     * By default will be 0, programs using denoiser should set this
     */
    int verbose;

    static void defineParams(XmippProgram *program);
    void readParams(XmippProgram *program);

    /** Apply the filter
     */
    void apply(MultidimArray< double >& img);
};

#endif
