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

#ifndef _BASIS_HH
#define _BASIS_HH

#include "blobs.h"
#include "splines.h"

#include <data/image.h>

const int BLOB_SUBSAMPLING = 10;
const int PIXEL_SUBSAMPLING = 1;
const int SPLINE_SUBSAMPLING = 1;

/**@defgroup BasisFunction Basis function
   @ingroup ReconsLibrary
    This class defines the basis function to use for the reconstruction.
    Currently, valid basis functions are blobs and voxels.
*/
//@{
/** Basis class. */
class Basis
{
public:
    /// Type of basis function
    typedef enum {blobs, voxels, splines} t_basis_function;

    /// Basis function to use
    t_basis_function type;

    /// Sampling rate
    double Tm;

    /// Blob parameters
    struct blobtype blob;

    /** Volume deformation matrix.
        See the documentation of Basic_ART_Parameters for further explanation. */
    Matrix2D<double> *D;
public:
    /// Empty constructor. By default, blobs
    Basis()
    {
        set_default();
    }

    /// Default values
    void set_default();

    /** Read parameters from a command line.
        This function reads the parameters from a command line
        defined by argc and argv. An exception might be thrown by any
        of the internal conversions, this would mean that there is
        an error in the command line and you might show a usage message. */
    void read(int argc, char **argv);

    /** Read parameters from a file.
        An exception is thrown if the file cannot be open */
    void read(const FileName &fn);

    /** Produce side information.
        You must provide the grid in which this basis function will live */
    void produce_side_info(const Grid &grid);

    /// Usage
    void usage() const;

    /// Show
    friend std::ostream & operator << (std::ostream &out, const Basis &basis);

    /** Set sampling rate. */
    void set_sampling_rate(double _Tm);

    /** Set D.
        D is the deformation matrix used for crystals. */
    void set_D(Matrix2D<double> *_D)
    {
        D = _D;
    }

    /** Max length of the basis.
        This is the maximum distance between the center of the basis and
        its further point. */
    double max_length() const
    {
        switch (type)
        {
        case blobs:
            return blob.radius;
            break;
        case voxels:
            return sqrt(3.0) * 0.5;
            break;
        case splines:
            return sqrt(3.0) * 2.0;
            break;
        }
        return 0.0;
    }

    /** Change basis to voxels.
        This function takes a grid volume in the basis indicated in this object
        and translates it into voxels of the given size. If the volume is
        already in voxels a padding is done so that the output is of the given
        size and the basis volume is in the center. */
    void changeToVoxels(GridVolume &vol_basis, Matrix3D<double> *vol_voxels,
                        int Zdim, int Ydim, int Xdim, int threads = 1 ) const;

    /** Change basis from voxels.
        A voxel volume is provided, then the output vol_basis will be shaped
        to represent this voxel volume. If the volume is already in voxels, then
        the mask and the radius mask are applied. */
    void changeFromVoxels(const Matrix3D<double> &vol_voxels,
                          GridVolume &vol_basis, int grid_type, double grid_relative_size,
                          const Matrix3D<double> *vol_mask,
                          const Matrix2D<double> *D, double R,int threads=1) const;

    /** Basis value at a given point. */
    double value_at(const Matrix1D<double> &r) const
    {
        double module_r;
        switch (type)
        {
        case (blobs):
            module_r = sqrt(XX(r) * XX(r) + YY(r) * YY(r) + ZZ(r) * ZZ(r));
            return blob_val(module_r, blob);
            break;
        case (voxels):
            if (-0.5 <= XX(r) && XX(r) < 0.5 &&
                -0.5 <= YY(r) && YY(r) < 0.5 &&
                -0.5 <= ZZ(r) && ZZ(r) < 0.5) return 1.0;
            else return 0.0;
            break;
        case (splines):
            if (-2 <= XX(r) && XX(r) < 2 &&
                -2 <= YY(r) && YY(r) < 2 &&
                -2 <= ZZ(r) && ZZ(r) < 2)
                return spatial_Bspline03LUT(r);
            else return 0.0;
            break;
        }
        return 0.0;
    }

    /** Projection at a given direction (u) with a given point (r). */
    double projection_at(const Matrix1D<double> &u, const Matrix1D<double> &r)
    const
{
        const double p0 = 1.0 / (2 * PIXEL_SUBSAMPLING) - 0.5;
        const double pStep = 1.0 / PIXEL_SUBSAMPLING;
        const double pAvg = 1.0 / (PIXEL_SUBSAMPLING * PIXEL_SUBSAMPLING);
        double module_r, px, py;
        Matrix1D<double> aux(3);
        int i, j;
        switch (type)
        {
            case (blobs):
                module_r = sqrt(XX(r) * XX(r) + YY(r) * YY(r) + ZZ(r) * ZZ(r));
                return blob_proj(module_r, blob);
                break;
            case (voxels):
            {
                double retval = 0;
                ZZ(aux) = ZZ(r);
                for (i = 0, px = p0; i < PIXEL_SUBSAMPLING; i++, px += pStep)
                {
                    XX(aux) = XX(r) + px;
                    for (j = 0, py = p0; j < PIXEL_SUBSAMPLING; j++, py += pStep)
                    {
                        YY(aux) = YY(r) + py;
                        retval += intersection_unit_cube(u, aux);
                    }
                }
                return retval*pAvg;
                break;
            }
            case (splines):
                return spatial_Bspline03_proj(r, u);
                break;
        }
        return 0.0;
    }
public:
    /// Blob footprint
    ImageOver       blobprint;

    /// Square of the footprint
    ImageOver       blobprint2;

    /// Sum of the basis on the grid points
    double          sum_on_grid;
};
//@}
#endif
