/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
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

#include "inria.h"

#ifdef _HAVE_INRIA
#include <typedefs.h>
#include <extrema.h>
#include <recline.h>

#include <data/volume.h>

/* Derivative -------------------------------------------------------------- */
void compute_derivative(const Matrix3D<double> &in_vol,
                        Matrix3D<double> &out_vol, char *type,
                        double sigma)
{
    // Prepare some variables
    int borderLengths[3] = {1, 1, 1};
    recursiveFilterType filterType = ALPHA_DERICHE;
    float filterCoefs[3];
    filterCoefs[0] = filterCoefs[1] = filterCoefs[2] = sigma;
    int Dims[3];
    Dims[2] = ZSIZE(in_vol);
    Dims[1] = YSIZE(in_vol);
    Dims[0] = XSIZE(in_vol);
    out_vol.resize(in_vol);

    // Translate derivative type
    derivativeOrder derivative_type[3];
#define COUNT(c,N) \
    N=0; \
    if (type[0]==c) {N++;} \
    if (type[0]!='\0') {if (type[1]==c) N++;} \
    if (type[0]!='\0') {if (type[1]!='\0') {if (type[2]==c) N++;}}
#define SET_DERIVATIVE_TYPE_FOR(c,idx) {\
        int N; COUNT(c, N); \
        switch (N) { \
        case 0: derivative_type[idx]=SMOOTHING; break; \
        case 1: derivative_type[idx]=DERIVATIVE_1; break; \
        case 2: derivative_type[idx]=DERIVATIVE_2; break; \
        case 3: derivative_type[idx]=DERIVATIVE_3; break; \
        } \
    }
    SET_DERIVATIVE_TYPE_FOR('X', 0);
    SET_DERIVATIVE_TYPE_FOR('Y', 1);
    SET_DERIVATIVE_TYPE_FOR('Z', 2);

    // Compute derivative
    Matrix3D<float> aux;
    typeCast(in_vol, aux);
    // This auxiliar volume is used due to the fact that INRIA library
    // seems to have a problem with input double volumes
    if (partial_derivative_3D(VOL_ARRAY(aux), FLOAT,
                              VOL_ARRAY(out_vol), DOUBLE,
                              Dims, borderLengths, filterCoefs, filterType, derivative_type) == 0)
        REPORT_ERROR(1, "Processing of derivative failed.");
}

/* Gradient ---------------------------------------------------------------- */
void compute_gradient(const Matrix3D<double> &in_vol,
                      Vectorial_Matrix3D &out_vol, double sigma)
{
    out_vol.resize(in_vol);
    compute_derivative(in_vol, out_vol.X(), "X", sigma);
    compute_derivative(in_vol, out_vol.Y(), "Y", sigma);
    compute_derivative(in_vol, out_vol.Z(), "Z", sigma);
}
#endif
