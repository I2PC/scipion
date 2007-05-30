/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "blobs.h"

#include <data/funcs.h>
#include <data/geometry.h>

/* Value of a blob --------------------------------------------------------- */
double kaiser_value(double r, double a, double alpha, int m)
{
    double rda, rdas, arg, w;

    rda = r / a;
    rdas = rda * rda;
    if (rdas <= 1.0)
    {
        arg = alpha * sqrt(1.0 - rdas);
        if (m == 0)
        {
            w = bessi0(arg) / bessi0(alpha);
        }
        else if (m == 1)
        {
            if (alpha == 0.0)
                w = 1.0 - rdas;
            else
                w = sqrt(1.0 - rdas) * bessi1(arg) / bessi1(alpha);
        }
        else if (m == 2)
        {
            if (alpha == 0.0)
                w = (1.0 - rdas) * (1.0 - rdas);
            else
                w = (1.0 - rdas) * bessi2(arg) / bessi2(alpha);
        }
        else REPORT_ERROR(1, "m out of range in kaiser_value()");

    }
    else
        w = 0.0;

    return w;
}

/* Line integral through a blob -------------------------------------------- */
/* Value of line integral through Kaiser-Bessel radial function
   (n >=2 dimensions) at distance s from center of function.
   Parameter m = 0, 1, or 2. */
double kaiser_proj(double s, double a, double alpha, int m)
{
    double sda, sdas, w, arg, p;

    sda = s / a;
    sdas = sda * sda;
    w = 1.0 - sdas;
    if (w > 1.0e-10)
    {
        arg = alpha * sqrt(w);
        if (m == 0)
        {
            if (alpha == 0.0)
                p = 2.0 * a * sqrt(w);
            else
                p = (2.0 * a / alpha) * sinh(arg) / bessi0(alpha);

        }
        else if (m == 1)
        {
            if (alpha == 0.0)
                p = 2.0 * a * w * sqrt(w) * (2.0 / 3.0);
            else
                p = (2.0 * a / alpha) * sqrt(w) * (cosh(arg) - sinh(arg) / arg)
                    / bessi1(alpha);

        }
        else if (m == 2)
        {
            if (alpha == 0.0)
                p = 2.0 * a * w * w * sqrt(w) * (8.0 / 15.0);
            else
                p = (2.0 * a / alpha) * w *
                    ((3.0 / (arg * arg) + 1.0) * sinh(arg) - (3.0 / arg) * cosh(arg)) / bessi2(alpha);
        }
        else REPORT_ERROR(2, "m out of range in kaiser_proj()");

    }
    else
        p = 0.0;

    return p;
}

/* Fourier value of a blob ------------------------------------------------- */
double kaiser_Fourier_value(double w, double a, double alpha, int m)
{
    if (m != 2)
        REPORT_ERROR(2, "m out of range in kaiser_Fourier_value()");
    double sigma = sqrt(ABS(alpha * alpha - (2 * PI * a * w) * (2 * PI * a * w)));
    if (2*PI*a*w > alpha)
        return  pow(2*PI, 3 / 2)*pow(a, 3)*pow(alpha, 2)*bessj3_5(sigma)
                / (bessi3_5(alpha)*pow(sigma, 3.5));
    else
        return  pow(2*PI, 3 / 2)*pow(a, 3)*pow(alpha, 2)*bessi3_5(sigma)
                / (bessi3_5(alpha)*pow(sigma, 3.5));
}


/* Sum a blob on a simple grid --------------------------------------------- */
// Computes sum of the values of a unitary blob on grid points. The blob is
// supposed to be at the origin of the absolute coordinate system
double sum_blob_SimpleGrid(const struct blobtype &blob, const SimpleGrid &grid,
                           const Matrix2D<double> *D)
{
    SPEED_UP_temps;
    Matrix1D<double> gr(3), ur(3), corner1(3), corner2(3);
    double         actual_radius;
    int          i, j, k;
    double        sum = 0.0;

// Compute the limits of the blob in the grid coordinate system
    grid.universe2grid(vectorR3(-blob.radius, -blob.radius, -blob.radius), corner1);
    grid.universe2grid(vectorR3(blob.radius, blob.radius, blob.radius), corner2);
    if (D != NULL)
        box_enclosing(corner1, corner2, *D, corner1, corner2);

// Compute the sum in the points inside the grid
// The integer part of the vectors is taken for not picking points
// just in the border of the blob, which we know they are 0.
    for (i = (int)corner1.X(); i <= (int)corner2.X(); i++)
        for (j = (int)corner1.Y(); j <= (int)corner2.Y(); j++)
            for (k = (int)corner1.Z(); k <= (int)corner2.Z(); k++)
            {
                VECTOR_R3(gr, i, j, k);
                grid.grid2universe(gr, ur);
                if (D != NULL)
                    M3x3_BY_V3x1(ur, *D, ur);
                actual_radius = ur.module();
                if (actual_radius < blob.radius)
                    sum += kaiser_value(actual_radius,
                                        blob.radius, blob.alpha, blob.order);
            }
    return sum;
}

/* Volume integral of a blob ----------------------------------------------- */
double  basvolume(double a, double alpha, int m, int n)
{
    double  hn, tpi, v;
    hn = 0.5 * n;
    tpi = 2.0 * PI;

    if (alpha == 0.0)
    {
        if ((n / 2)*2 == n)           /* n even                               */
            v = pow(tpi, hn) * in_zeroarg(n / 2 + m) / in_zeroarg(m);
        else                        /* n odd                                */
            v = pow(tpi, hn) * inph_zeroarg(n / 2 + m) / in_zeroarg(m);

    }
    else
    {                        /* alpha > 0.0                          */
        if ((n / 2)*2 == n)           /* n even                               */
            v = pow(tpi / alpha, hn) * i_n(n / 2 + m, alpha) / i_n(m, alpha);
        else                        /* n odd                                */
            v = pow(tpi / alpha, hn) * i_nph(n / 2 + m, alpha) / i_n(m, alpha);
    }

    return v * pow(a, (double)n);
}

/* Bessel function I_n (x),  n = 0, 1, 2, ...
 Use ONLY for small values of n     */

double i_n(int n, double x)
{
    int i;
    double i_ns1, i_n, i_np1;

    if (n == 0)   return bessi0(x);
    if (n == 1)   return bessi1(x);
    if (x == 0.0) return 0.0;
    i_ns1 = bessi0(x);
    i_n   = bessi1(x);

    for (i = 1; i < n; i++)
    {
        i_np1 = i_ns1 - (2 * i) / x * i_n;
        i_ns1 = i_n;
        i_n   = i_np1;
    }
    return i_n;
}

/*.....Bessel function I_(n+1/2) (x),  n = 0, 1, 2, ..........................*/

double i_nph(int n, double x)
{
    int i;
    double r2dpix;
    double i_ns1, i_n, i_np1;

    if (x == 0.0) return 0.0;
    r2dpix = sqrt(2.0 / (PI * x));
    i_ns1 = r2dpix * cosh(x);
    i_n   = r2dpix * sinh(x);

    for (i = 1; i <= n; i++)
    {
        i_np1 = i_ns1 - (2 * i - 1) / x * i_n;
        i_ns1 = i_n;
        i_n   = i_np1;
    }
    return i_n;
}

/*....Limit (z->0) of (1/z)^n I_n(z)..........................................*/
double in_zeroarg(int n)
{
    int i;
    double fact;
    fact = 1.0;
    for (i = 1; i <= n; i++)
    {
        fact *= 0.5 / i;
    }
    return fact;
}

/*.......Limit (z->0) of (1/z)^(n+1/2) I_(n+1/2) (z)..........................*/

double inph_zeroarg(int n)
{
    int i;
    double fact;
    fact = 1.0;
    for (i = 1; i <= n; i++)
    {
        fact *= 1.0 / (2 * i + 1.0);
    }
    return fact*sqrt(2.0 / PI);
}

/* Sum a blob on a complex grid -------------------------------------------- */
double sum_blob_Grid(const struct blobtype &blob, const Grid &grid,
                     const Matrix2D<double> *D)
{
    double sum = 0;
    for (int i = 0; i < grid.GridsNo(); i++)
        sum += sum_blob_SimpleGrid(blob, grid(i), D);
    return sum;
}

/* Zero freq --------------------------------------------------------------- */
double blob_freq_zero(struct blobtype b)
{
    return sqrt(b.alpha*b.alpha + 6.9879*6.9879) / (2*PI*b.radius);
}

/* Attenuation ------------------------------------------------------------- */
double blob_att(double w, struct blobtype b)
{
    return blob_Fourier_val(w, b) / blob_Fourier_val(0, b);
}

/* Number of operations ---------------------------------------------------- */
double blob_ops(double w, struct blobtype b)
{
    return pow(b.alpha*b.alpha + 6.9879*6.9879, 1.5) / b.radius;
}

/* Optimal grid ------------------------------------------------------------ */
double optimal_CC_grid_relative_size(struct blobtype b)
{
    return 1 / blob_freq_zero(b);
}
double optimal_BCC_grid_relative_size(struct blobtype b)
{
    return sqrt(8.0) / (2*blob_freq_zero(b));
}
double optimal_FCC_grid_relative_size(struct blobtype b)
{
    return sqrt(27.0) / (4*blob_freq_zero(b));
}

/* Best blob in region ----------------------------------------------------- */
//#define DEBUG
blobtype best_blob(double alpha_0, double alpha_F, double inc_alpha,
                   double a_0, double a_F, double inc_a, double w, double *target_att,
                   int target_length)
{
    blobtype retval;
    retval.order = 2;
    int alpha_size = FLOOR((alpha_F - alpha_0) / inc_alpha + 1);
    int a_size = FLOOR((a_F - a_0) / inc_a + 1);
    ImageXmipp att(alpha_size, a_size);
    ImageXmipp ops(alpha_size, a_size);
    int i, j, best_i, best_j;
    double a, alpha, best_a = -1, best_alpha = -1;
    double best_ops = 1e10, best_att = 0;
    for (i = 0, alpha = alpha_0; i < alpha_size; alpha += inc_alpha, i++)
        for (j = 0, a = a_0; j < a_size; a += inc_a, j++)
        {
            retval.radius = a;
            retval.alpha = alpha;
            att(i, j) = blob_att(w, retval);
            ops(i, j) = blob_ops(w, retval);
            if (j > 0)
                for (int n = target_length - 1; n >= 0; n--)
                    if (att(i, j - 1) > target_att[n] && att(i, j) < target_att[n])
                    {
                        att(i, j) = 0;
                        if (ops(i, j - 1) < best_ops && att(i, j - 1) >= best_att)
                        {
                            best_i = i;
                            best_j = j - 1;
                            best_alpha = alpha;
                            best_a = a - inc_a;
                            best_ops = ops(i, j - 1);
                            best_att = target_att[n];
                        }
                    }
        }
#ifdef DEBUG
    att.write((string)"att" + floatToString(w) + ".xmp");
    ops.write((string)"ops" + floatToString(w) + ".xmp");
#endif
    retval.radius = best_a;
    retval.alpha = best_alpha;
    return retval;
}

/* Footprint of a blob ----------------------------------------------------- */
void footprint_blob(
    ImageOver &blobprint,         // blob foot_print table
    const struct blobtype &blob,  // blob description
    int   istep,            // number of foot-print samples per one sample
    // on projection plane in u,v directions
    int   normalise)        // set to 1 if you want normalise. Usually
// you may omit it and no normalisation is performed
{
// Resize output image and redefine the origin of it
    int footmax = CEIL(blob.radius);
    blobprint.init(-footmax, footmax, istep, -footmax, footmax, istep);

// Run for Imge class indexes
    for (int i = STARTINGY(blobprint()); i <= FINISHINGY(blobprint()); i++)
        for (int j = STARTINGX(blobprint()); j <= FINISHINGX(blobprint()); j++)
        {
            // Compute oversampled index and blob value
            double vi, ui;
            IMG2OVER(blobprint, i, j, vi, ui);
            double r = sqrt(vi * vi + ui * ui);
            IMGPIXEL(blobprint, i, j) = blob_proj(r, blob);
        }

// Adjust the footprint structure
    if (normalise) blobprint() /= blobprint().sum();
}

//#define DEBUG
//#define DEBUG_MORE
/* Blobs -> Voxels for a SimpleGrid ---------------------------------------- */
// This function will construct a table of blob values (something like the
// footprint)
#define DEFORM_BLOB_WHEN_IN_CRYSTAL
void blobs2voxels_SimpleGrid(const Matrix3D<double> &vol_blobs,
                             const SimpleGrid &grid, const struct blobtype &blob,
                             Matrix3D<double> *vol_voxels, const Matrix2D<double> *D = NULL, int istep = 50,
                             Matrix3D<double> *vol_corr = NULL, const Matrix3D<double> *vol_mask = NULL,
                             bool FORW = true, int eq_mode = VARTK)
{
    Matrix2D<double> Dinv;                   // Inverse of D
    Matrix1D<double> act_coord(3);           // Coord: Actual position inside
    // the voxel volume without deforming
    Matrix1D<double> real_position(3);       // Coord: actual position after
    // applying the V transformation
    Matrix1D<double> beginZ(3);              // Coord: Voxel coordinates of the
    // blob at the 3D point
    // (z0,YY(lowest),XX(lowest))
    Matrix1D<double> beginY(3);              // Coord: Voxel coordinates of the
    // blob at the 3D point
    // (z0,y0,XX(lowest))
    Matrix1D<double> corner2(3), corner1(3); // Coord: Corners of the
    // blob in the voxel volume
    Matrix1D<double> gcurrent(3);            // Position in g of current point
    Matrix1D<double> blob_table;             // Something like a blobprint
    // but with the values of the
    // blob in space
    double         d;                        // Distance between the center
    // of the blob and a voxel position
    int           id;                        // index inside the blob value
    // table for tha blob value at
    // a distance d
    double         intx, inty, intz;         // Nearest integer voxel
    int           i, j, k;                   // Index within the blob volume
    int           process;                   // True if this blob has to be
    // processed
    double         vol_correction;           // Correction to apply to the
    // volume when "projecting" back
    SPEED_UP_temps;

    // Some aliases
#define x0 STARTINGX(*vol_voxels)
#define xF FINISHINGX(*vol_voxels)
#define y0 STARTINGY(*vol_voxels)
#define yF FINISHINGY(*vol_voxels)
#define z0 STARTINGZ(*vol_voxels)
#define zF FINISHINGZ(*vol_voxels)

#ifdef DEBUG
    bool condition = !FORW;
    if (condition)
    {
        (*vol_voxels)().printShape();
        cout << endl;
        cout << "x0= " << x0 << " xF= " << xF << endl;
        cout << "y0= " << y0 << " yF= " << yF << endl;
        cout << "z0= " << z0 << " zF= " << zF << endl;
        cout << grid;
    }
#endif

    // Invert deformation matrix ............................................
    if (D != NULL) Dinv = D->inv();

    // Compute a blob value table ...........................................
    blob_table.resize((int)(blob.radius*istep + 1));
    for (i = 0; i < blob_table.xdim; i++)
    {
        VEC_ELEM(blob_table, i) = blob_val((double)i / istep, blob);
#ifdef DEBUG_MORE
        if (condition)
            cout << "Blob (" << i << ") r=" << (double)i / istep <<
            " val= " << VEC_ELEM(blob_table, i) << endl;
#endif
    }

    // Convert the whole grid ...............................................
    // Corner of the plane defined by Z. These coordinates are in the
    // universal coord. system
    grid.grid2universe(grid.lowest, beginZ);

    Matrix1D<double> grid_index(3);
    for (k = (int) ZZ(grid.lowest); k <= (int) ZZ(grid.highest); k++)
    {
        // Corner of the row defined by Y
        beginY = beginZ;
        for (i = (int) YY(grid.lowest); i <= (int) YY(grid.highest); i++)
        {
            // First point in the row
            act_coord = beginY;
            for (j = (int) XX(grid.lowest); j <= (int) XX(grid.highest); j++)
            {
                VECTOR_R3(grid_index, j, i, k);
#ifdef DEBUG
                if (condition)
                {
                    printf("Dealing blob at (%d,%d,%d) = %f\n", j, i, k, VOL_ELEM(vol_blobs, k, i, j));
                    cout << "Center of the blob      "
                    << act_coord.transpose() << endl;
                }
#endif

                // Place act_coord in its right place
                if (D != NULL)
                {
                    M3x3_BY_V3x1(real_position, *D, act_coord);
#ifdef DEBUG
                    if (condition)
                        cout << "Center of the blob moved to "
                        //ROB, the "moved" coordinates are in
                        // real_position not in act_coord
                        << act_coord.transpose() << endl;
                    << real_position.transpose() << endl;
#endif
                    // ROB This is OK if blob.radius is in Cartesian space as I
                    // think is the case
                }
                else real_position = act_coord;

                // These two corners are also real valued
                process = true;
//ROB
//This is OK if blob.radius is in Cartesian space as I think is the case
                V3_PLUS_CT(corner1, real_position, -blob.radius);
                V3_PLUS_CT(corner2, real_position, blob.radius);
#ifdef DEFORM_BLOB_WHEN_IN_CRYSTAL
                //ROB
                //we do not need this, it is already in Cartesian space
                //if (D!=NULL)
                //   box_enclosing(corner1,corner2, *D, corner1, corner2);
#endif

                if (XX(corner1) >= xF) process = false;
                if (YY(corner1) >= yF) process = false;
                if (ZZ(corner1) >= zF) process = false;
                if (XX(corner2) <= x0) process = false;
                if (YY(corner2) <= y0) process = false;
                if (ZZ(corner2) <= z0) process = false;
#ifdef DEBUG
                if (!process && condition) cout << "   It is outside output volume\n";
#endif
                if (!grid.is_interesting(real_position))
                {
#ifdef DEBUG
                    if (process && condition) cout << "   It is not interesting\n";
#endif
                    process = false;
                }

#ifdef DEBUG
                if (condition)
                {
                    cout << "Corner 1 for this point " << corner1.transpose() << endl;
                    cout << "Corner 2 for this point " << corner2.transpose() << endl;
                }
#endif

                if (process)
                {
                    // Clip the corners to the volume borders
                    XX(corner1) = ROUND(CLIP(XX(corner1), x0, xF));
                    YY(corner1) = ROUND(CLIP(YY(corner1), y0, yF));
                    ZZ(corner1) = ROUND(CLIP(ZZ(corner1), z0, zF));
                    XX(corner2) = ROUND(CLIP(XX(corner2), x0, xF));
                    YY(corner2) = ROUND(CLIP(YY(corner2), y0, yF));
                    ZZ(corner2) = ROUND(CLIP(ZZ(corner2), z0, zF));
#ifdef DEBUG
                    if (condition)
                    {
                        cout << "Clipped and rounded Corner 1 " << corner1.transpose() << endl;
                        cout << "Clipped and rounded Corner 2 " << corner2.transpose() << endl;
                    }
#endif

                    if (!FORW)
                        switch (eq_mode)
                        {
                        case VARTK:
                            vol_correction = 0;
                            break;
                        case VMAXARTK:
                            vol_correction = -1e38;
                            break;
                        }

                    // Effectively convert
                    long N_eq;
                    N_eq = 0;
                    for (intz = ZZ(corner1); intz <= ZZ(corner2); intz++)
                        for (inty = YY(corner1); inty <= YY(corner2); inty++)
                            for (intx = XX(corner1); intx <= XX(corner2); intx++)
                            {
                                int iz = (int)intz, iy = (int)inty, ix = (int)intx;
                                if (vol_mask != NULL)
                                    if (!VOL_ELEM(*vol_mask, iz, iy, ix)) continue;

                                // Compute distance to the center of the blob
                                VECTOR_R3(gcurrent, intx, inty, intz);
#ifdef DEFORM_BLOB_WHEN_IN_CRYSTAL
                                // ROB
                                //if (D!=NULL)
                                //   M3x3_BY_V3x1(gcurrent,Dinv,gcurrent);
#endif
                                V3_MINUS_V3(gcurrent, real_position, gcurrent);
                                d = sqrt(XX(gcurrent) * XX(gcurrent) +
                                         YY(gcurrent) * YY(gcurrent) +
                                         ZZ(gcurrent) * ZZ(gcurrent));
                                if (d > blob.radius) continue;
                                id = (int)(d * istep);
#ifdef DEBUG_MORE
                                if (condition)
                                {
                                    cout << "At (" << intx << ","
                                    << inty << "," << intz << ") distance=" << d;
                                    cout.flush();
                                }
#endif

                                // Add at that position the corresponding blob value

                                if (FORW)
                                {
                                    VOL_ELEM(*vol_voxels, iz, iy, ix) +=
                                        VOL_ELEM(vol_blobs, k, i, j) *
                                        VEC_ELEM(blob_table, id);
#ifdef DEBUG_MORE
                                    if (condition)
                                    {
                                        cout << " adding " << VOL_ELEM(vol_blobs, k, i, j)
                                        << " * " << VEC_ELEM(blob_table, id) << " = "
                                        << VOL_ELEM(vol_blobs, k, i, j)*
                                        VEC_ELEM(blob_table, id) << endl;
                                        cout.flush();
                                    }
#endif
                                    if (vol_corr != NULL)
                                        VOL_ELEM(*vol_corr, iz, iy, ix) +=
                                            VEC_ELEM(blob_table, id) * VEC_ELEM(blob_table, id);
                                }
                                else
                                {
                                    double contrib = VOL_ELEM(*vol_corr, iz, iy, ix) *
                                                     VEC_ELEM(blob_table, id);
                                    switch (eq_mode)
                                    {
                                    case VARTK:
                                        vol_correction += contrib;
                                        N_eq++;
                                        break;
                                    case VMAXARTK:
                                        if (contrib > vol_correction)
                                            vol_correction = contrib;
                                        break;

                                    }
#ifdef DEBUG_MORE
                                    if (condition)
                                    {
                                        cout << " adding " << VOL_ELEM(*vol_corr, iz, iy, ix)
                                        << " * " << VEC_ELEM(blob_table, id) << " = "
                                        << contrib << endl;
                                        cout.flush();
                                    }
#endif
                                }
                            }
                    if (N_eq == 0) N_eq = 1;
                    if (!FORW)
                    {
                        VOL_ELEM(vol_blobs, k, i, j) += vol_correction / N_eq;
#ifdef DEBUG_MORE
                        cout << " correction= " << vol_correction << endl
                        << " Number of eqs= " << N_eq << endl
                        << " Blob after correction= "
                        << VOL_ELEM(vol_blobs, k, i, j) << endl;
#endif
                    }
                }

                // Prepare for next iteration
                XX(act_coord) = XX(act_coord) + grid.relative_size * MAT_ELEM(grid.basis, 0, 0);
                YY(act_coord) = YY(act_coord) + grid.relative_size * MAT_ELEM(grid.basis, 1, 0);
                ZZ(act_coord) = ZZ(act_coord) + grid.relative_size * MAT_ELEM(grid.basis, 2, 0);
            }
            XX(beginY) = XX(beginY) + grid.relative_size * MAT_ELEM(grid.basis, 0, 1);
            YY(beginY) = YY(beginY) + grid.relative_size * MAT_ELEM(grid.basis, 1, 1);
            ZZ(beginY) = ZZ(beginY) + grid.relative_size * MAT_ELEM(grid.basis, 2, 1);
        }
        XX(beginZ) = XX(beginZ) + grid.relative_size * MAT_ELEM(grid.basis, 0, 2);
        YY(beginZ) = YY(beginZ) + grid.relative_size * MAT_ELEM(grid.basis, 1, 2);
        ZZ(beginZ) = ZZ(beginZ) + grid.relative_size * MAT_ELEM(grid.basis, 2, 2);
    }
}
#undef x0
#undef y0
#undef z0
#undef xF
#undef yF
#undef zF
#undef DEBUG
#undef DEBUG_MORE

/* Voxel volume shape ------------------------------------------------------ */
//#define DEBUG
void voxel_volume_shape(const GridVolume &vol_blobs,
                        const struct blobtype &blob, const Matrix2D<double> *D,
                        Matrix1D<int> &corner1, Matrix1D<int> &size)
{
    Matrix1D<double>  Gcorner1(3),  Gcorner2(3);     // lowest and highest coord.

    corner1.resize(3);
    size.resize(3);

    // Look for the lowest and highest volume coordinate
    vol_blobs.grid().voxel_corners(Gcorner1, Gcorner2, D);

    // Add blob radius in each direction, and find the furthest integer
    // samples => compute size
    Gcorner1 = CEILnD(Gcorner1 - blob.radius);
    Gcorner2 = FLOORnD(Gcorner2 + blob.radius);

    XX(size) = (int)XX(Gcorner2) - (int)XX(Gcorner1) + 1;
    YY(size) = (int)YY(Gcorner2) - (int)YY(Gcorner1) + 1;
    ZZ(size) = (int)ZZ(Gcorner2) - (int)ZZ(Gcorner1) + 1;
#ifdef DEBUG
    cout << "Gcorner1  " << Gcorner1.transpose() << endl;
    cout << "Gcorner2  " << Gcorner2.transpose() << endl;
    cout << "Size of voxel volume " << (int)ZZ(size) << " x "
    << (int)YY(size) << " x " << (int)XX(size) << endl;
#endif

#ifdef NEVER_DEFINED
    // In principle this limitation has been substituted by a direct
    // specification of the output volume size, and should no longer
    // be valid. However it is a beatiful piece of code to be removed
    // already
    if (limit != 0 && XX(size) != limit)
    {
        double diff = ((double)XX(size) - (double)limit) / 2;
        if (diff == (int)diff)
        {
            Gcorner1 += diff;
            Gcorner2 -= diff;
        }
        else
        {
            Gcorner1 += (diff - 0.5);
            Gcorner2 -= (diff + 0.5);
        }
        XX(size) = (int)XX(Gcorner2) - (int)XX(Gcorner1) + 1;
        YY(size) = (int)YY(Gcorner2) - (int)YY(Gcorner1) + 1;
        ZZ(size) = (int)ZZ(Gcorner2) - (int)ZZ(Gcorner1) + 1;
#ifdef DEBUG
        cout << "Limiting to " << limit << " diff = " << diff << endl;
        cout << "New Gcorner1  " << Gcorner1.transpose() << endl;
        cout << "New Gcorner2  " << Gcorner2.transpose() << endl;
#endif
    }
#endif

    type_cast(Gcorner1, corner1);

#ifdef DEBUG
    cout << "Final size of voxel volume " << (int)ZZ(size) << " x "
    << (int)YY(size) << " x " << (int)XX(size) << endl;
    cout << "Corner1= " << corner1.transpose() << endl;
#endif
}
#undef DEBUG

/* Blobs -> Voxels for a Grid ---------------------------------------------- */
//#define DEBUG
void blobs2voxels(const GridVolume &vol_blobs,
                  const struct blobtype &blob, Matrix3D<double> *vol_voxels,
                  const Matrix2D<double> *D, int Zdim, int Ydim, int Xdim)
{

    // Resize and set starting corner .......................................
    if (Zdim == 0 || Ydim == 0 || Xdim == 0)
    {
        Matrix1D<int> size, corner;
        voxel_volume_shape(vol_blobs, blob, D, corner, size);
        (*vol_voxels).initZeros(ZZ(size), YY(size), XX(size));
        (*vol_voxels).startingX() = XX(corner);
        (*vol_voxels).startingY() = YY(corner);
        (*vol_voxels).startingZ() = ZZ(corner);
    }
    else
    {
        (*vol_voxels).initZeros(Zdim, Ydim, Xdim);
        (*vol_voxels).setXmippOrigin();
    }

    // Convert each subvolume ...............................................
    for (int i = 0; i < vol_blobs.VolumesNo(); i++)
    {
        blobs2voxels_SimpleGrid(vol_blobs(i)(), vol_blobs.grid(i),
                                blob, vol_voxels, D);
#ifdef DEBUG
        cout << "Blob grid no " << i << " stats: ";
        vol_blobs(i)().print_stats();
        cout << endl;
        cout << "So far vol stats: ";
        (*vol_voxels).print_stats();
        cout << endl;
        VolumeXmipp save;
        save() = *vol_voxels;
        save.write((string)"PPPvoxels" + ItoA(i));
#endif
    }

    // Now normalise the resulting volume ..................................
    double inorm = 1.0 / sum_blob_Grid(blob, vol_blobs.grid(), D); // Aqui tambien hay que multiplicar ****!!!!
    FOR_ALL_ELEMENTS_IN_MATRIX3D(*vol_voxels)
    VOL_ELEM(*vol_voxels, k, i, j) *= inorm;

    // Set voxels outside interest region to minimum value .................
    double R = vol_blobs.grid(0).get_interest_radius();
    if (R != -1)
    {
        double R2 = (R - 6) * (R - 6);

        // Compute minimum value within sphere
        double min_val = VOL_ELEM(*vol_voxels, 0, 0, 0);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(*vol_voxels)
        if (j*j + i*i + k*k <= R2 - 4)
            min_val = MIN(min_val, VOL_ELEM(*vol_voxels, k, i, j));

        // Substitute minimum value
        R2 = (R - 2) * (R - 2);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(*vol_voxels)
        if (j*j + i*i + k*k >= R2)
            VOL_ELEM(*vol_voxels, k, i, j) = min_val;
    }
}
#undef DEBUG

/* Blobs -> Coefs ---------------------------------------------------------- */
//#define DEBUG
void blobs2space_coefficients(const GridVolume &vol_blobs,
                              const struct blobtype &blob, Matrix3D<double> *vol_coefs)
{

    // Compute vol_coefs shape
    Matrix1D<int> corner1, size;
    voxel_volume_shape(vol_blobs, blob, NULL, corner1, size);
    double g_2 = vol_blobs.grid(0).relative_size / 2;
    XX(corner1) = (int)(FLOOR(XX(corner1)) / g_2);
    XX(size) = (int)(CEIL(XX(size)) / g_2);
    YY(corner1) = (int)(FLOOR(YY(corner1)) / g_2);
    YY(size) = (int)(CEIL(YY(size)) / g_2);
    ZZ(corner1) = (int)(FLOOR(ZZ(corner1)) / g_2);
    ZZ(size) = (int)(CEIL(ZZ(size)) / g_2);
    (*vol_coefs).initZeros(ZZ(size), YY(size), XX(size));
    (*vol_coefs).startingX() = XX(corner1);
    (*vol_coefs).startingY() = YY(corner1);
    (*vol_coefs).startingZ() = ZZ(corner1);

    // Set all blob coefficients at the right position
    for (int n = 0; n < vol_blobs.VolumesNo(); n++)
    {
        int ZZ_lowest = (int)ZZ(vol_blobs.grid(n).lowest);
        int YY_lowest = (int)YY(vol_blobs.grid(n).lowest);
        int XX_lowest = (int)XX(vol_blobs.grid(n).lowest);
        int ZZ_highest = (int)ZZ(vol_blobs.grid(n).highest);
        int YY_highest = (int)YY(vol_blobs.grid(n).highest);
        int XX_highest = (int)XX(vol_blobs.grid(n).highest);
        for (int k = ZZ_lowest; k <= ZZ_highest; k++)
            for (int i = YY_lowest; i <= YY_highest; i++)
                for (int j = XX_lowest; j <= XX_highest; j++)
                {
                    Matrix1D<double> grid_index(3), univ_position(3),
                    coef_position(3);
                    VECTOR_R3(grid_index, j, i, k);
                    vol_blobs.grid(n).grid2universe(grid_index, univ_position);
                    V3_BY_CT(coef_position, univ_position, 1.0 / g_2);
                    (*vol_coefs)((int)ZZ(coef_position), (int)YY(coef_position),
                                 (int)XX(coef_position)) = vol_blobs(n)(k, i, j);
#ifdef DEBUG
                    cout << "Blob value at (" << j << "," << i << ","
                    << k << ") (" << XX(univ_position)
                    << "," << YY(univ_position) << ","
                    << ZZ(univ_position) << ") ("
                    << XX(coef_position) << "," << YY(coef_position) << ","
                    << ZZ(coef_position) << ") --> "
                    << vol_blobs(n)((int)k, (int)i, (int)j) << endl;
#endif
                }
    }
}
#undef DEBUG

/* Voxels -> blobs --------------------------------------------------------- */
/* This function is very similar to the ART reconstruction process from
   projections, look there for another point of view of the ART process. */
#define FORWARD  true
#define BACKWARD false
//#define DEBUG
void ART_voxels2blobs_single_step(
    GridVolume &vol_in,                 // Input solution
    GridVolume *vol_out,                // Output solution, might be the same
    // as the input one
    const struct blobtype &blob,        // blob
    const Matrix2D<double> *D,          // deformation matrix
    double lambda,                      // ART lambda
    Matrix3D<double> *theo_vol,         // Theoretical volume
    const Matrix3D<double> *read_vol,   // Volume we want to translate to blobs
    Matrix3D<double> *corr_vol,         // Normalizing volume
    const Matrix3D<double> *mask_vol,   // Mask volume, 1 if that voxel must
    // be counted as a true equation
    double &mean_error,                 // Output mean error
    double &max_error,                  // Output maximum error in a voxel
    int eq_mode                         // Equation mode
)
{

    // Resize output volumes ................................................
    if (read_vol != NULL)
    {
        (*theo_vol).initZeros(*read_vol);
    }
    else if (mask_vol != NULL)
    {
        (*theo_vol).initZeros(*mask_vol);
    }
    else
    {
        REPORT_ERROR(1,
                     "ART_voxels2blobs_single_step: Mask and voxel volumes are empty");
    }
    (*corr_vol).initZeros(*theo_vol);

    // Translate actual blob volume to voxels ...............................
    for (int i = 0; i < vol_in.VolumesNo(); i++)
    {
        blobs2voxels_SimpleGrid(vol_in(i)(), vol_in.grid(i), blob, theo_vol, D,
                                50, corr_vol, mask_vol, FORWARD, eq_mode);
#ifdef DEBUG
        cout << "Blob grid no " << i << " stats: ";
        vol_in(i)().print_stats();
        cout << endl;
        cout << "So far vol stats: ";
        (*theo_vol)().print_stats();
        cout << endl;
#endif
    }

    // Now normalise the resulting volume ..................................
    double norm = sum_blob_Grid(blob, vol_in.grid(), D); // Aqui tambien hay que multiplicar ****!!!!
    FOR_ALL_ELEMENTS_IN_MATRIX3D(*theo_vol)
    VOL_ELEM(*theo_vol, k, i, j) /= norm;

#ifdef DEBUG
    VolumeXmipp save, save2;
    save() = *theo_vol;
    save.write("PPPtheovol.vol");
    cout << "Theo stats:";
    save().print_stats();
    cout << endl;
    save() = *corr_vol;
    save.write("PPPcorr2vol.vol");
    save2().resize(save());
#endif

    // Compute differences ..................................................
    mean_error = 0;
    double read_val;
    if (read_vol != NULL) read_val = VOL_ELEM(*read_vol, 0, 0, 0);
    else                read_val = 0;
    max_error = ABS(read_val - VOL_ELEM(*theo_vol, 0, 0, 0));

    double diff;
    int N = 0;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(*theo_vol)
    {
        // Compute difference volume and error
        if (read_vol != NULL) read_val = VOL_ELEM(*read_vol, k, i, j);
        else                read_val = 0;

        if (mask_vol == NULL)
        {
            diff = read_val - VOL_ELEM(*theo_vol, k, i, j);
            N++;
        }
        else
            if (VOL_ELEM(*mask_vol, k, i, j) == 1)
            {
                diff = read_val - VOL_ELEM(*theo_vol, k, i, j);
                N++;
            }
            else
                diff = 0;

        max_error = MAX(max_error, ABS(diff));
        mean_error += diff * diff;
#ifdef DEBUG
        save(k, i, j) = diff;
        save2(k, i, j) = read_val;
#endif


        // Compute the correction volume
        if (ABS(VOL_ELEM(*corr_vol, k, i, j)) < 1e-2)
            VOL_ELEM(*corr_vol, k, i, j) = SGN(VOL_ELEM(*corr_vol, k, i, j));
        VOL_ELEM(*corr_vol, k, i, j) =
            lambda * diff / VOL_ELEM(*corr_vol, k, i, j);
    }
#ifdef DEBUG
    save.write("PPPdiffvol.vol");
    cout << "Diff stats:";
    save().print_stats();
    cout << endl;
    save2.write("PPPreadvol.vol");
    cout << "Read stats:";
    save2().print_stats();
    cout << endl;
    save() = *corr_vol;
    save.write("PPPcorrvol.vol");
    cout << "Corr stats:";
    save().print_stats();
    cout << endl;
#endif

    mean_error /= MAX(N, 1); // At worst, divided by 1

    // Backprojection of correction volume ..................................
    for (int i = 0; i < vol_in.VolumesNo(); i++)
    {
        blobs2voxels_SimpleGrid((*vol_out)(i)(), (*vol_out).grid(i), blob,
                                theo_vol, D, 50, corr_vol, mask_vol, BACKWARD, eq_mode);
#ifdef DEBUG
        cout << "Blob grid no " << i << " stats: ";
        vol_in(i)().print_stats();
        cout << endl;
#endif
    }
#ifdef DEBUG
    char c;
    cout << "Press any key to continue\n";
    cin >> c;
#endif
}
#undef DEBUG

//#define DEBUG
void voxels2blobs(const Matrix3D<double> *vol_voxels,
                  const struct blobtype &blob,
                  GridVolume &vol_blobs, int grid_type, double grid_relative_size,
                  double lambda, const Matrix3D<double> *vol_mask,
                  const Matrix2D<double> *D, double final_error_change,
                  int tell, double R)
{
    VolumeXmipp theo_vol, corr_vol;
    double mean_error, mean_error_1, max_error;
    int it = 1;

    tell = SHOW_CONVERSION;

    // Resize output volume .................................................
    Grid grid_blobs;
    if (R == -1)
    {
        Matrix1D<double> corner1(3), corner2(3);
        XX(corner1) = STARTINGX(*vol_voxels);
        YY(corner1) = STARTINGY(*vol_voxels);
        ZZ(corner1) = STARTINGZ(*vol_voxels);
        XX(corner2) = FINISHINGX(*vol_voxels);
        YY(corner2) = FINISHINGY(*vol_voxels);
        ZZ(corner2) = FINISHINGZ(*vol_voxels);

        switch (grid_type)
        {
        case (CC):
                        grid_blobs = Create_CC_grid(grid_relative_size, corner1, corner2);
            break;
        case (FCC):
                        grid_blobs = Create_FCC_grid(grid_relative_size, corner1, corner2);
            break;
        case (BCC):
                        grid_blobs = Create_BCC_grid(grid_relative_size, corner1, corner2);
            break;
        }
    }
    else
{
        switch (grid_type)
        {
        case (CC):
                        grid_blobs = Create_CC_grid(grid_relative_size, R);
            break;
        case (FCC):
                        grid_blobs = Create_FCC_grid(grid_relative_size, R);
            break;
        case (BCC):
                        grid_blobs = Create_BCC_grid(grid_relative_size, R);
            break;
        }
    }
    vol_blobs.adapt_to_grid(grid_blobs);

    // Convert ..............................................................
    cout << "Converting voxel volume to blobs\n";
    ART_voxels2blobs_single_step(vol_blobs, &vol_blobs, blob, D, lambda,
                                 &(theo_vol()), vol_voxels,
                                 &(corr_vol()), vol_mask, mean_error_1, max_error, VARTK);
    if (tell && SHOW_CONVERSION)
        cout << "   Finished iteration: " << it++
        << " Mean Error= " << mean_error_1
        << " Max_Error= " << max_error
        << endl;
    else
        cout << "0%";
    cout.flush();
#ifdef DEBUG
    theo_vol.write("PPPtheo.vol");
    corr_vol.write("PPPcorr.vol");
    cout << "Press any key\n";
    char c;
    cin >> c;
#endif

    double change;
    bool end_condition;
    do
{
        ART_voxels2blobs_single_step(vol_blobs, &vol_blobs, blob, D, lambda,
                                     &(theo_vol()), vol_voxels,
                                     &(corr_vol()), vol_mask, mean_error, max_error, VARTK);
#ifdef DEBUG
        theo_vol.write("PPPtheo.vol");
        corr_vol.write("PPPcorr.vol");
        cout << "Press any key\n";
        char c;
        cin >> c;
#endif
        change = ABS(mean_error - mean_error_1) / mean_error_1;
        mean_error_1 = mean_error;
        if (tell && SHOW_CONVERSION)
            cout << "   Finished iteration: " << it++
            << " Mean Error= " << mean_error
            << " Max_Error= " << max_error
            << endl;
        else
        {
            printf("\r");
            cout << 100 - change << "%";
        }
        cout.flush();
        end_condition = change <= final_error_change;
    }
    while (!end_condition);
    cout << endl;
}
#undef DEBUG
