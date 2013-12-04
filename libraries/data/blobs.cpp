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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "blobs.h"

#include "xmipp_funcs.h"
#include "geometry.h"

pthread_mutex_t blobs_conv_mutex = PTHREAD_MUTEX_INITIALIZER;

int * slices_status;
int slices_processed;

/* Value of a blob --------------------------------------------------------- */
double kaiser_value(double r, double a, double alpha, int m)
{
    double rda, rdas, arg, w;

    rda = r / a;
    if (rda <= 1.0)
    {
        rdas = rda * rda;
        arg = alpha * sqrt(1.0 - rdas);
        if (m == 0)
        {
            w = bessi0(arg) / bessi0(alpha);
        }
        else if (m == 1)
        {
            w = sqrt (1.0 - rdas);
            if (alpha != 0.0)
                w *= bessi1(arg) / bessi1(alpha);
        }
        else if (m == 2)
        {
            w = sqrt (1.0 - rdas);
            w = w * w;
            if (alpha != 0.0)
                w *= bessi2(arg) / bessi2(alpha);
        }
        else if (m == 3)
        {
            w = sqrt (1.0 - rdas);
            w = w * w * w;
            if (alpha != 0.0)
                w *= bessi3(arg) / bessi3(alpha);
        }
        else if (m == 4)
        {
            w = sqrt (1.0 - rdas);
            w = w * w * w *w;
            if (alpha != 0.0)
                w *= bessi4(arg) / bessi4(alpha);
        }
        else
            REPORT_ERROR(ERR_VALUE_INCORRECT, "m out of range in kaiser_value()");

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
        else
            REPORT_ERROR(ERR_VALUE_INCORRECT, "m out of range in kaiser_proj()");

    }
    else
        p = 0.0;

    return p;
}

/* Fourier value of a blob ------------------------------------------------- */
double kaiser_Fourier_value(double w, double a, double alpha, int m)
{
    if (m != 2 && m !=0)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "m out of range in kaiser_Fourier_value()");
    double sigma = sqrt(ABS(alpha * alpha - (2. * PI * a * w) * (2. * PI * a * w)));
    if (m == 2)
    {
        if (2.*PI*a*w > alpha)
            return  pow(2.*PI, 3. / 2.)*pow(a, 3.)*pow(alpha, 2.)*bessj3_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 3.5));
        else
            return  pow(2.*PI, 3. / 2.)*pow(a, 3.)*pow(alpha, 2.)*bessi3_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 3.5));
    }
    else if (m == 0)
    {
        if (2*PI*a*w > alpha)
            return  pow(2.*PI, 3. / 2.)*pow(a, 3)*bessj1_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 1.5));
        else
            return  pow(2.*PI, 3. / 2.)*pow(a, 3)*bessi1_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 1.5));
    }
    else
    	REPORT_ERROR(ERR_ARG_INCORRECT,"Invalid blob order");
}


/* Sum a blob on a simple grid --------------------------------------------- */
// Computes sum of the values of a unitary blob on grid points. The blob is
// supposed to be at the origin of the absolute coordinate system
double sum_blob_SimpleGrid(const struct blobtype &blob, const SimpleGrid &grid,
                           const Matrix2D<double> *D)
{
    SPEED_UP_temps012;
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
    for (i = (int)XX(corner1); i <= (int)XX(corner2); i++)
        for (j = (int)YY(corner1); j <= (int)YY(corner2); j++)
            for (k = (int)ZZ(corner1); k <= (int)ZZ(corner2); k++)
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

    if (n == 0)
        return bessi0(x);
    if (n == 1)
        return bessi1(x);
    if (x == 0.0)
        return 0.0;
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

    if (x == 0.0)
        return 0.0;
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
    for (size_t i = 0; i < grid.GridsNo(); i++)
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
    Image<double> att, ops;
    att().resize(alpha_size, a_size);
    ops().resize(alpha_size, a_size);
    int i, j;
    double a, alpha, best_a = -1, best_alpha = -1;
    double best_ops = 1e10, best_att = 0;
    for (i = 0, alpha = alpha_0; i < alpha_size; alpha += inc_alpha, i++)
        for (j = 0, a = a_0; j < a_size; a += inc_a, j++)
        {
            retval.radius = a;
            retval.alpha = alpha;
            A2D_ELEM(att(), i, j) = blob_att(w, retval);
            A2D_ELEM(ops(), i, j) = blob_ops(w, retval);
            if (j > 0)
                for (int n = target_length - 1; n >= 0; n--)
                    if (A2D_ELEM(att(), i, j - 1) > target_att[n] &&
                        A2D_ELEM(att(), i, j) < target_att[n])
                    {
                        A2D_ELEM(att(), i, j) = 0;
                        if (A2D_ELEM(ops(), i, j - 1) < best_ops &&
                            A2D_ELEM(att(), i, j - 1) >= best_att)
                        {
                            best_alpha = alpha;
                            best_a = a - inc_a;
                            best_ops = A2D_ELEM(ops(), i, j - 1);
                            best_att = target_att[n];
                        }
                    }
        }
#ifdef DEBUG
    att.write((std::string)"att" + floatToString(w) + ".xmp");
    ops.write((std::string)"ops" + floatToString(w) + ".xmp");
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
    if (normalise)
        blobprint() /= blobprint().sum();
}

//#define DEBUG
//#define DEBUG_MORE
/* Blobs -> Voxels for a SimpleGrid ---------------------------------------- */
// This function will construct a table of blob values (something like the
// footprint)
#define DEFORM_BLOB_WHEN_IN_CRYSTAL
void * blobs2voxels_SimpleGrid( void * data )
{
    ThreadBlobsToVoxels * thread_data = (ThreadBlobsToVoxels *) data;

    const MultidimArray<double> *vol_blobs = thread_data->vol_blobs;
    const SimpleGrid *grid = thread_data->grid;
    const struct blobtype *blob = thread_data->blob;
    MultidimArray<double> *vol_voxels = thread_data->vol_voxels;
    const Matrix2D<double> *D = thread_data->D;
    int istep = thread_data->istep;
    MultidimArray<double> *vol_corr = thread_data->vol_corr;
    const MultidimArray<double> *vol_mask = thread_data->vol_mask;
    ;
    bool FORW = thread_data->FORW;
    int eq_mode = thread_data->eq_mode;

    int min_separation = thread_data->min_separation;

    int z_planes = (int)(ZZ(grid->highest) - ZZ(grid->lowest) + 1);

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
    MultidimArray<double> blob_table;             // Something like a blobprint
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
    double         vol_correction=0;         // Correction to apply to the
    // volume when "projecting" back
    SPEED_UP_temps012;

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
        std::cout << std::endl;
        std::cout << "x0= " << x0 << " xF= " << xF << std::endl;
        std::cout << "y0= " << y0 << " yF= " << yF << std::endl;
        std::cout << "z0= " << z0 << " zF= " << zF << std::endl;
        std::cout << grid;
    }
#endif

    // Invert deformation matrix ............................................
    if (D != NULL)
        Dinv = D->inv();

    // Compute a blob value table ...........................................
    blob_table.resize((int)(blob->radius*istep + 1));
    for (size_t i = 0; i < blob_table.xdim; i++)
    {
        A1D_ELEM(blob_table, i) = kaiser_value((double)i/istep, blob->radius, blob->alpha, blob->order);

#ifdef DEBUG_MORE

        if (condition)
            std::cout << "Blob (" << i << ") r=" << (double)i / istep <<
            " val= " << A1D_ELEM(blob_table, i) << std::endl;
#endif

    }

    int assigned_slice;

    do
    {
        assigned_slice = -1;
        do
        {
            pthread_mutex_lock(&blobs_conv_mutex);
            if( slices_processed == z_planes )
            {
                pthread_mutex_unlock(&blobs_conv_mutex);
                return (void*)NULL;
            }

            for(int w = 0 ; w < z_planes ; w++ )
            {
                if( slices_status[w]==0 )
                {
                    slices_status[w] = -1;
                    assigned_slice = w;
                    slices_processed++;

                    for( int in = (w-min_separation+1) ; in <= (w+min_separation-1 ) ; in ++ )
                    {
                        if( in != w )
                        {
                            if( ( in >= 0 ) && ( in < z_planes ))
                            {
                                if( slices_status[in] != -1 )
                                    slices_status[in]++;
                            }
                        }
                    }
                    break;
                }
            }

            pthread_mutex_unlock(&blobs_conv_mutex);
        }
        while( assigned_slice == -1);

        // Convert the whole grid ...............................................
        // Corner of the plane defined by Z. These coordinates are in the
        // universal coord. system
        Matrix1D<double> aux( grid->lowest );
        k = (int)(assigned_slice + ZZ( grid->lowest ));
        ZZ(aux) = k;
        grid->grid2universe(aux, beginZ);

        Matrix1D<double> grid_index(3);

        // Corner of the row defined by Y
        beginY = beginZ;
        for (i = (int) YY(grid->lowest); i <= (int) YY(grid->highest); i++)
        {
            // First point in the row
            act_coord = beginY;
            for (j = (int) XX(grid->lowest); j <= (int) XX(grid->highest); j++)
            {
                VECTOR_R3(grid_index, j, i, k);
#ifdef DEBUG

                if (condition)
                {
                    printf("Dealing blob at (%d,%d,%d) = %f\n", j, i, k, A3D_ELEM(*vol_blobs, k, i, j));
                    std::cout << "Center of the blob      "
                    << act_coord.transpose() << std::endl;
                }
#endif

                // Place act_coord in its right place
                if (D != NULL)
                {
                    M3x3_BY_V3x1(real_position, *D, act_coord);
#ifdef DEBUG

                    if (condition)
                        std::cout << "Center of the blob moved to "
                        //ROB, the "moved" coordinates are in
                        // real_position not in act_coord
                        << act_coord.transpose() << std::endl;
                    << real_position.transpose() << std::endl;
#endif
                    // ROB This is OK if blob.radius is in Cartesian space as I
                    // think is the case
                }
                else
                    real_position = act_coord;

                // These two corners are also real valued
                process = true;
                //ROB
                //This is OK if blob.radius is in Cartesian space as I think is the case
                V3_PLUS_CT(corner1, real_position, -blob->radius);
                V3_PLUS_CT(corner2, real_position, blob->radius);
#ifdef DEFORM_BLOB_WHEN_IN_CRYSTAL
                //ROB
                //we do not need this, it is already in Cartesian space
                //if (D!=NULL)
                //   box_enclosing(corner1,corner2, *D, corner1, corner2);
#endif

                if (XX(corner1) >= xF)
                    process = false;
                if (YY(corner1) >= yF)
                    process = false;
                if (ZZ(corner1) >= zF)
                    process = false;
                if (XX(corner2) <= x0)
                    process = false;
                if (YY(corner2) <= y0)
                    process = false;
                if (ZZ(corner2) <= z0)
                    process = false;
#ifdef DEBUG

                if (!process && condition)
                    std::cout << "   It is outside output volume\n";
#endif

                if (!grid->is_interesting(real_position))
                {
#ifdef DEBUG
                    if (process && condition)
                        std::cout << "   It is not interesting\n";
#endif

                    process = false;
                }

#ifdef DEBUG
                if (condition)
                {
                    std::cout << "Corner 1 for this point " << corner1.transpose() << std::endl;
                    std::cout << "Corner 2 for this point " << corner2.transpose() << std::endl;
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
                        std::cout << "Clipped and rounded Corner 1 " << corner1.transpose() << std::endl;
                        std::cout << "Clipped and rounded Corner 2 " << corner2.transpose() << std::endl;
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
                                    if (!A3D_ELEM(*vol_mask, iz, iy, ix))
                                        continue;

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
                                if (d > blob->radius)
                                    continue;
                                id = (int)(d * istep);
#ifdef DEBUG_MORE

                                if (condition)
                                {
                                    std::cout << "At (" << intx << ","
                                    << inty << "," << intz << ") distance=" << d;
                                    std::cout.flush();
                                }
#endif

                                // Add at that position the corresponding blob value

                                if (FORW)
                                {
                                    A3D_ELEM(*vol_voxels, iz, iy, ix) +=
                                        A3D_ELEM(*vol_blobs, k, i, j) *
                                        A1D_ELEM(blob_table, id);
#ifdef DEBUG_MORE

                                    if (condition)
                                    {
                                        std::cout << " adding " << A3D_ELEM(*vol_blobs, k, i, j)
                                        << " * " << A1D_ELEM(blob_table, id) << " = "
                                        << A3D_ELEM(*vol_blobs, k, i, j)*
                                        A1D_ELEM(blob_table, id) << std::endl;
                                        std::cout.flush();
                                    }
#endif
                                    if (vol_corr != NULL)
                                        A3D_ELEM(*vol_corr, iz, iy, ix) +=
                                            A1D_ELEM(blob_table, id) * A1D_ELEM(blob_table, id);
                                }
                                else
                                {
                                    double contrib = A3D_ELEM(*vol_corr, iz, iy, ix) *
                                                     A1D_ELEM(blob_table, id);
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
                                        std::cout << " adding " << A3D_ELEM(*vol_corr, iz, iy, ix)
                                        << " * " << A1D_ELEM(blob_table, id) << " = "
                                        << contrib << std::endl;
                                        std::cout.flush();
                                    }
#endif

                                }
                            }
                    if (N_eq == 0)
                        N_eq = 1;
                    if (!FORW)
                    {
                        A3D_ELEM(*vol_blobs, k, i, j) += vol_correction / N_eq;
#ifdef DEBUG_MORE

                        std::cout << " correction= " << vol_correction << std::endl
                        << " Number of eqs= " << N_eq << std::endl
                        << " Blob after correction= "
                        << A3D_ELEM(*vol_blobs, k, i, j) << std::endl;
#endif

                    }
                }

                // Prepare for next iteration
                XX(act_coord) = XX(act_coord) + grid->relative_size * (grid->basis)( 0, 0);
                YY(act_coord) = YY(act_coord) + grid->relative_size * (grid->basis)( 1, 0);
                ZZ(act_coord) = ZZ(act_coord) + grid->relative_size * (grid->basis)( 2, 0);
            }
            XX(beginY) = XX(beginY) + grid->relative_size * (grid->basis)( 0, 1);
            YY(beginY) = YY(beginY) + grid->relative_size * (grid->basis)( 1, 1);
            ZZ(beginY) = ZZ(beginY) + grid->relative_size * (grid->basis)( 2, 1);
        }

        pthread_mutex_lock(&blobs_conv_mutex);

        for( int in = (assigned_slice-min_separation+1) ; in <= (assigned_slice+min_separation-1); in ++ )
        {
            if( in != assigned_slice )
            {
                if( ( in >= 0 ) && ( in < z_planes))
                {
                    if( slices_status[in] != -1 )
                        slices_status[in]--;
                }
            }
        }

        pthread_mutex_unlock(&blobs_conv_mutex);
    }
    while(1);

    return (void*)NULL;
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
    Gcorner1 -= blob.radius;
    Gcorner1.selfCEIL();
    Gcorner2 += blob.radius;
    Gcorner2.selfFLOOR();

    XX(size) = (int)XX(Gcorner2) - (int)XX(Gcorner1) + 1;
    YY(size) = (int)YY(Gcorner2) - (int)YY(Gcorner1) + 1;
    ZZ(size) = (int)ZZ(Gcorner2) - (int)ZZ(Gcorner1) + 1;
#ifdef DEBUG

    std::cout << "Gcorner1  " << Gcorner1.transpose() << std::endl;
    std::cout << "Gcorner2  " << Gcorner2.transpose() << std::endl;
    std::cout << "Size of voxel volume " << (int)ZZ(size) << " x "
    << (int)YY(size) << " x " << (int)XX(size) << std::endl;
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

        std::cout << "Limiting to " << limit << " diff = " << diff << std::endl;
        std::cout << "New Gcorner1  " << Gcorner1.transpose() << std::endl;
        std::cout << "New Gcorner2  " << Gcorner2.transpose() << std::endl;
#endif

    }
#endif

    typeCast(Gcorner1, corner1);

#ifdef DEBUG

    std::cout << "Final size of voxel volume " << (int)ZZ(size) << " x "
    << (int)YY(size) << " x " << (int)XX(size) << std::endl;
    std::cout << "Corner1= " << corner1.transpose() << std::endl;
#endif
}
#undef DEBUG

/* Blobs -> Voxels for a Grid ---------------------------------------------- */
//#define DEBUG
void blobs2voxels(const GridVolume &vol_blobs,
                  const struct blobtype &blob, MultidimArray<double> *vol_voxels,
                  const Matrix2D<double> *D, int threads, int Zdim, int Ydim, int Xdim)
{

    // Resize and set starting corner .......................................
    if (Zdim == 0 || Ydim == 0 || Xdim == 0)
    {
        Matrix1D<int> size, corner;
        voxel_volume_shape(vol_blobs, blob, D, corner, size);
        (*vol_voxels).initZeros(ZZ(size), YY(size), XX(size));
        STARTINGX(*vol_voxels) = XX(corner);
        STARTINGY(*vol_voxels) = YY(corner);
        STARTINGZ(*vol_voxels) = ZZ(corner);
    }
    else
    {
        (*vol_voxels).initZeros(Zdim, Ydim, Xdim);
        (*vol_voxels).setXmippOrigin();
    }

    pthread_t * th_ids = new pthread_t [threads];
    ThreadBlobsToVoxels * threads_d = new ThreadBlobsToVoxels [threads];

    // Convert each subvolume ...............................................
    for (size_t i = 0; i < vol_blobs.VolumesNo(); i++)
    {
        int min_distance = (int)ceil((2*(vol_blobs.grid(i)).relative_size ) / blob.radius ) + 1;

        slices_status = (int *)malloc(sizeof(int)*(int)((ZZ((&(vol_blobs.grid(i)))->highest)-ZZ((&(vol_blobs.grid(i)))->lowest)+1)));
        memset(slices_status,0,sizeof(int)*(int)((ZZ((&(vol_blobs.grid(i)))->highest)-ZZ((&(vol_blobs.grid(i)))->lowest)+1)));
        slices_processed = 0;

        for( int c = 0 ; c < threads ; c++ )
        {
            threads_d[c].vol_blobs = &(vol_blobs(i)());
            threads_d[c].grid = &(vol_blobs.grid(i));
            threads_d[c].blob = &blob;
            threads_d[c].vol_voxels = vol_voxels;
            threads_d[c].D = D;
            threads_d[c].istep = 50;
            threads_d[c].vol_corr = NULL;
            threads_d[c].vol_mask = NULL;
            threads_d[c].FORW = true;
            threads_d[c].eq_mode = VARTK;
            threads_d[c].thread_id = c;
            threads_d[c].threads_num = threads;
            threads_d[c].min_separation = min_distance;

            pthread_create( (th_ids+c), NULL, blobs2voxels_SimpleGrid, (void *)(threads_d+c) );
        }

        // Wait for threads to finish
        for( int c = 0 ; c < threads ; c++ )
        {
            pthread_join(*(th_ids+c),NULL);
        }

#ifdef DEBUG
        std::cout << "Blob grid no " << i << " stats: ";
        vol_blobs(i)().printStats();
        std::cout << std::endl;
        std::cout << "So far vol stats: ";
        (*vol_voxels).printStats();
        std::cout << std::endl;
        Image<double> save;
        save() = *vol_voxels;
        save.write((std::string)"PPPvoxels" + integerToString(i));
#endif

        free( slices_status );
    }

    // Now normalise the resulting volume ..................................
    double inorm = 1.0 / sum_blob_Grid(blob, vol_blobs.grid(), D); // Aqui tambien hay que multiplicar ****!!!!
    FOR_ALL_ELEMENTS_IN_ARRAY3D(*vol_voxels)
    A3D_ELEM(*vol_voxels, k, i, j) *= inorm;

    // Set voxels outside interest region to minimum value .................
    double R = vol_blobs.grid(0).get_interest_radius();
    if (R != -1)
    {
        double R2 = (R - 6) * (R - 6);

        // Compute minimum value within sphere
        double min_val = A3D_ELEM(*vol_voxels, 0, 0, 0);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(*vol_voxels)
        if (j*j + i*i + k*k <= R2 - 4)
            min_val = XMIPP_MIN(min_val, A3D_ELEM(*vol_voxels, k, i, j));

        // Substitute minimum value
        R2 = (R - 2) * (R - 2);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(*vol_voxels)
        if (j*j + i*i + k*k >= R2)
            A3D_ELEM(*vol_voxels, k, i, j) = min_val;
    }

    delete[] th_ids;
    delete[] threads_d;
}
#undef DEBUG

/* Blobs -> Coefs ---------------------------------------------------------- */
//#define DEBUG
void blobs2space_coefficients(const GridVolume &vol_blobs,
                              const struct blobtype &blob, MultidimArray<double> *vol_coefs)
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
    STARTINGX(*vol_coefs) = XX(corner1);
    STARTINGY(*vol_coefs) = YY(corner1);
    STARTINGZ(*vol_coefs) = ZZ(corner1);

    // Set all blob coefficients at the right position
    for (size_t n = 0; n < vol_blobs.VolumesNo(); n++)
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
                    A3D_ELEM((*vol_coefs),
                             (int)ZZ(coef_position),
                             (int)YY(coef_position),
                             (int)XX(coef_position) ) = A3D_ELEM(vol_blobs(n)(), k, i, j);
#ifdef DEBUG

                    std::cout << "Blob value at (" << j << "," << i << ","
                    << k << ") (" << XX(univ_position)
                    << "," << YY(univ_position) << ","
                    << ZZ(univ_position) << ") ("
                    << XX(coef_position) << "," << YY(coef_position) << ","
                    << ZZ(coef_position) << ") --> "
                    << A3D_ELEM(vol_blobs(n)(), (int)k, (int)i, (int)j) << std::endl;
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
    MultidimArray<double> *theo_vol,         // Theoretical volume
    const MultidimArray<double> *read_vol,   // Volume we want to translate to blobs
    MultidimArray<double> *corr_vol,         // Normalizing volume
    const MultidimArray<double> *mask_vol,   // Mask volume, 1 if that voxel must
    // be counted as a true equation
    double &mean_error,                 // Output mean error
    double &max_error,                  // Output maximum error in a voxel
    int eq_mode,                         // Equation mode
    int threads
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
        REPORT_ERROR(ERR_MATRIX_EMPTY,
                     "ART_voxels2blobs_single_step: Mask and voxel volumes are empty");
    }
    (*corr_vol).initZeros(*theo_vol);

    pthread_t * th_ids = (pthread_t *)malloc( threads * sizeof( pthread_t));
    ThreadBlobsToVoxels * threads_d = (ThreadBlobsToVoxels *) malloc ( threads * sizeof( ThreadBlobsToVoxels ) );

    // Translate actual blob volume to voxels ...............................
    for (size_t i = 0; i < vol_in.VolumesNo(); i++)
    {
        int min_distance = (int)ceil((2*(vol_in.grid(i)).relative_size ) / blob.radius ) + 1;

        slices_status = (int *)malloc(sizeof(int)*(int)((ZZ((&(vol_in.grid(i)))->highest)-ZZ((&(vol_in.grid(i)))->lowest)+1)));
        memset(slices_status,0,sizeof(int)*(int)((ZZ((&(vol_in.grid(i)))->highest)-ZZ((&(vol_in.grid(i)))->lowest)+1)));
        slices_processed = 0;

        for( int c = 0 ; c < threads ; c++ )
        {
            threads_d[c].vol_blobs = &(vol_in(i)());
            threads_d[c].grid = &(vol_in.grid(i));
            threads_d[c].blob = &blob;
            threads_d[c].vol_voxels = theo_vol;
            threads_d[c].D = D;
            threads_d[c].istep = 50;
            threads_d[c].vol_corr = corr_vol;
            threads_d[c].vol_mask = mask_vol;
            threads_d[c].FORW = true;
            threads_d[c].eq_mode = eq_mode;
            threads_d[c].thread_id = c;
            threads_d[c].threads_num = threads;
            threads_d[c].min_separation = min_distance;

            pthread_create( (th_ids+c), NULL, blobs2voxels_SimpleGrid, (void *)(threads_d+c) );
        }

        // Wait for threads to finish
        for( int c = 0 ; c < threads ; c++ )
        {
            pthread_join(*(th_ids+c),NULL);
        }

        free( slices_status );
        //        blobs2voxels_SimpleGrid(vol_in(i)(), vol_in.grid(i), blob, theo_vol, D,
        //                                50, corr_vol, mask_vol, FORWARD, eq_mode);
#ifdef DEBUG

        std::cout << "Blob grid no " << i << " stats: ";
        vol_in(i)().printStats();
        std::cout << std::endl;
        std::cout << "So far vol stats: ";
        (*theo_vol)().printStats();
        std::cout << std::endl;
#endif

    }

    // Now normalise the resulting volume ..................................
    double norm = sum_blob_Grid(blob, vol_in.grid(), D); // Aqui tambien hay que multiplicar ****!!!!
    FOR_ALL_ELEMENTS_IN_ARRAY3D(*theo_vol)
    A3D_ELEM(*theo_vol, k, i, j) /= norm;

#ifdef DEBUG

    Image<double> save, save2;
    save() = *theo_vol;
    save.write("PPPtheovol.vol");
    std::cout << "Theo stats:";
    save().printStats();
    std::cout << std::endl;
    save() = *corr_vol;
    save.write("PPPcorr2vol.vol");
    save2().resize(save());
#endif

    // Compute differences ..................................................
    mean_error = 0;
    double read_val;
    if (read_vol != NULL)
        read_val = A3D_ELEM(*read_vol, 0, 0, 0);
    else
        read_val = 0;
    max_error = ABS(read_val - A3D_ELEM(*theo_vol, 0, 0, 0));

    double diff;
    int N = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(*theo_vol)
    {
        // Compute difference volume and error
        if (read_vol != NULL)
            read_val = A3D_ELEM(*read_vol, k, i, j);
        else
            read_val = 0;

        if (mask_vol == NULL)
        {
            diff = read_val - A3D_ELEM(*theo_vol, k, i, j);
            N++;
        }
        else
            if (A3D_ELEM(*mask_vol, k, i, j) == 1)
            {
                diff = read_val - A3D_ELEM(*theo_vol, k, i, j);
                N++;
            }
            else
                diff = 0;

        max_error = XMIPP_MAX(max_error, ABS(diff));
        mean_error += diff * diff;
#ifdef DEBUG

        save(k, i, j) = diff;
        save2(k, i, j) = read_val;
#endif


        // Compute the correction volume
        if (ABS(A3D_ELEM(*corr_vol, k, i, j)) < 1e-2)
            A3D_ELEM(*corr_vol, k, i, j) = SGN(A3D_ELEM(*corr_vol, k, i, j));
        A3D_ELEM(*corr_vol, k, i, j) =
            lambda * diff / A3D_ELEM(*corr_vol, k, i, j);
    }
#ifdef DEBUG
    save.write("PPPdiffvol.vol");
    std::cout << "Diff stats:";
    save().printStats();
    std::cout << std::endl;
    save2.write("PPPreadvol.vol");
    std::cout << "Read stats:";
    save2().printStats();
    std::cout << std::endl;
    save() = *corr_vol;
    save.write("PPPcorrvol.vol");
    std::cout << "Corr stats:";
    save().printStats();
    std::cout << std::endl;
#endif

    mean_error /= XMIPP_MAX(N, 1); // At worst, divided by 1

    // Backprojection of correction volume ..................................
    for (size_t i = 0; i < vol_in.VolumesNo(); i++)
    {
        slices_status = (int *)malloc(sizeof(int)*(int)((ZZ((&(vol_out->grid(i)))->highest)-ZZ((&(vol_out->grid(i)))->lowest)+1)));
        memset(slices_status,0,sizeof(int)*(int)((ZZ((&(vol_out->grid(i)))->highest)-ZZ((&(vol_out->grid(i)))->lowest)+1)));
        slices_processed = 0;

        for( int c = 0 ; c < threads ; c++ )
        {
            threads_d[c].vol_blobs = &((*vol_out)(i)());
            threads_d[c].grid = &(vol_out->grid(i));
            threads_d[c].blob = &blob;
            threads_d[c].vol_voxels = theo_vol;
            threads_d[c].D = D;
            threads_d[c].istep = 50;
            threads_d[c].vol_corr = corr_vol;
            threads_d[c].vol_mask = mask_vol;
            threads_d[c].FORW = false;
            threads_d[c].eq_mode = eq_mode;
            threads_d[c].thread_id = c;
            threads_d[c].threads_num = threads;
            threads_d[c].min_separation = 1;

            pthread_create( (th_ids+c), NULL, blobs2voxels_SimpleGrid, (void *)(threads_d+c) );
        }

        // Wait for threads to finish
        for( int c = 0 ; c < threads ; c++ )
        {
            pthread_join(*(th_ids+c),NULL);
        }

        free( slices_status );
        //        blobs2voxels_SimpleGrid((*vol_out)(i)(), (*vol_out).grid(i), blob,
        //                                theo_vol, D, 50, corr_vol, mask_vol, BACKWARD, eq_mode);
#ifdef DEBUG

        std::cout << "Blob grid no " << i << " stats: ";
        vol_in(i)().printStats();
        std::cout << std::endl;
#endif

    }
#ifdef DEBUG
    char c;
    std::cout << "Press any key to continue\n";
    std::cin >> c;
#endif
}
#undef DEBUG

//#define DEBUG
void voxels2blobs(const MultidimArray<double> *vol_voxels,
                  const struct blobtype &blob,
                  GridVolume &vol_blobs, int grid_type, double grid_relative_size,
                  double lambda, const MultidimArray<double> *vol_mask,
                  const Matrix2D<double> *D, double final_error_change,
                  int tell, double R, int threads)
{
    Image<double> theo_vol, corr_vol;
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
    std::cout << "Converting voxel volume to blobs\n";
    ART_voxels2blobs_single_step(vol_blobs, &vol_blobs, blob, D, lambda,
                                 &(theo_vol()), vol_voxels,
                                 &(corr_vol()), vol_mask, mean_error_1, max_error, VARTK, threads);
    if (tell && SHOW_CONVERSION)
        std::cout << "   Finished iteration: " << it++
        << " Mean Error= " << mean_error_1
        << " Max_Error= " << max_error
        << std::endl;
    else
        std::cout << "0%";
    std::cout.flush();
#ifdef DEBUG

    theo_vol.write("PPPtheo.vol");
    corr_vol.write("PPPcorr.vol");
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
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
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
#endif

        change = ABS(mean_error - mean_error_1) / mean_error_1;
        mean_error_1 = mean_error;
        if (tell && SHOW_CONVERSION)
            std::cout << "   Finished iteration: " << it++
            << " Mean Error= " << mean_error
            << " Max_Error= " << max_error
            << std::endl;
        else
        {
            printf("\r");
            std::cout << 100 - change << "%";
        }
        std::cout.flush();
        end_condition = change <= final_error_change;
    }
    while (!end_condition);
    std::cout << std::endl;
}
#undef DEBUG
