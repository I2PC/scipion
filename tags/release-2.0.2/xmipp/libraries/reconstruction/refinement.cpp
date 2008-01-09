/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "refinement.h"

#include <data/funcs.h>
#include <data/geometry.h>
#include <data/fft.h>

//-------------------------------------------------------------------------
/* Correlate two projections and find the maximun of the correlation matrix -- */
/* If the maximun if moved further away than max_step returns 0 */
void calculate_and_find_correlation_max_proj(Projection const &proj1,
        Projection const &proj2,
        Projection &proj_temp,
        double &shift_X, double &shift_Y,
        double const max_step,
        int ref_trans_after, int imagen_no)
{
//       #define DEBUG_calculate_and_find_correlation_max_proj
#ifdef DEBUG_calculate_and_find_correlation_max_proj
    cout << "\n  (cal_find_corr_proj) imagen_no:  "   << imagen_no  << endl;
#endif


    proj_temp().resize(proj1());
    calculate_and_find_correlation_max_mat(IMGMATRIX(proj1),
                                           IMGMATRIX(proj2),
                                           IMGMATRIX(proj_temp),
                                           shift_X, shift_Y, max_step);

}//calculate_and_find_correlation_max end
#undef DEBUG_calculate_and_find_correlation_max_proj
