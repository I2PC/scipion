/***************************************************************************
 * Authors:     AUTHOR_NAME (jlvilas@cnb.csic.es)
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

#ifndef ASSIGNMENT_TILT_PAIR_H_
#define ASSIGNMENT_TILT_PAIR_H_
#define PI 3.14159265

#include <data/xmipp_program.h>
#include <math.h>
#include <alglib/src/ap.h>

class ProgassignmentTiltPair: public XmippProgram
{


public:
    /** Filenames */
    FileName fntilt, fnuntilt,  fndir, fnmic;

    /** Maximum shift*/
    int mshift;

    /** Particle size*/
    double particle_size;

    /** threshold*/
    double thr;

    /** Tilt estimation angle*/
    double tiltest;

    MetaData mdPartial;

public:

    void defineParams();

    void readParams();

    void search_affine_transform(float u1x, float u1y, float u2x, float u2y, float u3x, float u3y, float t1x,
    		float t1y, float t2x, float t2y, float t3x, float t3y,
    		Matrix1D<double> ux, Matrix1D<double> uy, size_t Xdim, size_t Ydim, struct Delaunay_T &delaunay_tilt, int &bestInliers,
    		Matrix2D<double> &A_coarse, Matrix1D<double> &T_coarse, bool contingency, int thrs);

    void findMaximumMinimum(const float u1, const float u2, const float u3, double &u_max, double &u_min);
    bool checkwindow(const float t1, const float t2, const float t3,
    				const double u_max, const double u_min);
    void run();


};
#endif
