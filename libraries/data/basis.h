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

#ifndef _BASIS_HH
#define _BASIS_HH

#include "blobs.h"
#include "splines.h"
#include "xmipp_image_over.h"
#include "xmipp_program.h"

const int BLOB_SUBSAMPLING = 10;
const int PIXEL_SUBSAMPLING = 1;
const int SPLINE_SUBSAMPLING = 1;

/**@defgroup BasisFunction Basis function
   @ingroup DataLibrary
    This class defines the basis function to use for the reconstruction.
    Currently, valid basis functions are blobs and voxels.
*/
//@{
/** Basis class. */
class Basis
{
public:
    /// Type of basis function
    typedef enum {blobs, voxels, splines} tBasisFunction;

    /// Basis function to use
    tBasisFunction type;

    /// Footprint is convolved with a volume PSF // At this moment only used with blobs
    MultidimArray<double> *VolPSF; // If NULL then standard blob is used

    /// Sampling rate
    double Tm;

    /// Blob parameters
    struct blobtype blob;

    /// Relative size for the grid
    double grid_relative_size;

    /** Volume deformation matrix.
        See the documentation of BasicARTParameters for further explanation. */
    Matrix2D<double> *D;

    /// Blob footprint
    ImageOver       blobprint;

    /// Square of the footprint
    ImageOver       blobprint2;

    /// Sum of the basis on the grid points
    double          sum_on_grid;

    /// Auxiliary vector for projections
    Matrix1D<double> aux;

public:
    /// Empty constructor. By default, blobs
    Basis();

    /// Default values
    void setDefault();

    /// Basis name
    String basisName() const;

    /** Read parameters from a command line.
        This function reads the parameters from a command line
        defined by argc and argv. An exception might be thrown by any
        of the internal conversions, this would mean that there is
        an error in the command line and you might show a usage message. */
    void read(int argc, char **argv);

    /** Read parameters from a file.
        An exception is thrown if the file cannot be open */
    void read(const FileName &fn);

    /**  Definition of paramaters
     */
    static void defineParams(XmippProgram * program, const char* prefix=NULL, const char* comment=NULL);

    /** Read the parameters from the command line
     */
    void readParams(XmippProgram * program);

    /** Produce side information.
        You must provide the grid in which this basis function will live */
    void produceSideInfo(const Grid &grid);

    /// Show
    friend std::ostream & operator << (std::ostream &out, const Basis &basis);

    /** Set sampling rate. */
    void setSamplingRate(double _Tm);

    /** Set D.
        D is the deformation matrix used for crystals. */
    void setD(Matrix2D<double> *_D)
    {
        D = _D;
    }

    /** Max length of the basis.
        This is the maximum distance between the center of the basis and
        its further point. */
    double maxLength() const;

    /** Change basis to voxels.
        This function takes a grid volume in the basis indicated in this object
        and translates it into voxels of the given size. If the volume is
        already in voxels a padding is done so that the output is of the given
        size and the basis volume is in the center. */
    void changeToVoxels(GridVolume &vol_basis, MultidimArray<double> *vol_voxels,
                        int Zdim, int Ydim, int Xdim, int threads = 1 ) const;

    /** Change basis from voxels.
        A voxel volume is provided, then the output vol_basis will be shaped
        to represent this voxel volume. If the volume is already in voxels, then
        the mask and the radius mask are applied. */
    void changeFromVoxels(const MultidimArray<double> &vol_voxels,
                          GridVolume &vol_basis, int grid_type, double grid_relative_size,
                          const MultidimArray<double> *vol_mask,
                          const Matrix2D<double> *D, double R,int threads=1) const;

    /** Basis value at a given point. */
    double valueAt(const Matrix1D<double> &r) const;

    /** Projection at a given direction (u) with a given point (r). */
    double projectionAt(const Matrix1D<double> &u, const Matrix1D<double> &r)
    const;
};

//@}
#endif
