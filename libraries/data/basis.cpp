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

#include "basis.h"
#include "xmipp_fftw.h"

Basis::Basis()
{
    setDefault();
}

// Set default -------------------------------------------------------------
void Basis::setDefault()
{
    type        = blobs;
    blob.radius = 2;
    blob.order  = 2;
    blob.alpha  = 10.4;
    VolPSF   = NULL;
    D           = NULL;
    blobprint.clear();
    blobprint2.clear();
    aux.resizeNoCopy(3);
}

// Basis name --------------------------------------------------------------
String Basis::basisName() const
{
    switch (type)
    {
    case Basis::blobs:
        return "blobs";
        break;
    case Basis::voxels:
        return "voxels";
        break;
    case Basis::splines:
        return "splines";
        break;
    default:
    	return "Unknown";
    	break;
    }
}

void Basis::defineParams(XmippProgram * program, const char* prefix, const char* comment)
{
    char tempLine[256];

    if(prefix == NULL)
        sprintf(tempLine, "  [--basis <basis_type=blobs>] ");
    else
        sprintf(tempLine,"%s --basis <basis_type=blobs> ", prefix);
    if (comment != NULL)
        sprintf(tempLine, "%s : %s", tempLine, comment);
    else
        sprintf(tempLine, "%s : Basis function to use for the reconstruction", tempLine);

    program->addParamsLine(tempLine);
    program->addParamsLine("    where <basis_type> ");
    program->addParamsLine("      blobs <radius=2> <Bessel_order=2> <alpha_param=10.4> : Default blob parameters and grid relative size adjusted to use small blobs");
    program->addParamsLine("     voxels");
    program->addParamsLine("     splines");
    program->addParamsLine(" or   --big_blobs    : blob parameters and grid relative size adjusted to use big blobs");
    program->addParamsLine(" or   --visual_blobs : blobs optimal for direct visualization");
}

void Basis::readParams(XmippProgram * program)
{
    String basisType = program->getParam("--basis");

    if (basisType == "blobs")
    {
        blob.radius = program->getDoubleParam("--basis", "blobs", 0);
        blob.order  = program->getIntParam("--basis", "blobs", 1);
        blob.alpha  = program->getDoubleParam("--basis", "blobs", 2);
    }
    if (!program->checkParam("--basis")) // Default is for small blobs
        grid_relative_size = 1.41;

    if (program->checkParam("--big_blobs"))
    {
        blob.radius = 2;
        blob.order  = 2;
        blob.alpha  = 3.6;

        grid_relative_size = 2.26;
    }
    else if (program->checkParam("--visual_blobs"))
    {
        blob.radius = 2.4;
        blob.order  = 2;
        blob.alpha  = 13.3633;

        grid_relative_size = 1.41;
    }
    else
    {
        if (basisType == "voxels")
        {
            type = voxels;
            grid_relative_size = 1;
        }
        else if (basisType == "splines")
        {
            type = splines;
            grid_relative_size = 1;
        }
    }
}

// Set sampling rate -------------------------------------------------------
void Basis::setSamplingRate(double _Tm)
{
    if (_Tm == 0)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "Basis::set_sampling_rate: Sampling rate cannot be 0");
    Tm = _Tm;
    blob.radius /= Tm;
}

// Produce side information ------------------------------------------------
void Basis::produceSideInfo(const Grid &grid)
{
    switch (type)
    {
    case (blobs):
        {
            int subsampling = (VolPSF == NULL )? BLOB_SUBSAMPLING : PIXEL_SUBSAMPLING;

            footprint_blob(blobprint, blob, subsampling);
            sum_on_grid = sum_blob_Grid(blob, grid, D);
            blobprint()  /= sum_on_grid;

            if (VolPSF != NULL)
            {
                // let adjust to the same resolution and size both blobprint and VolPSF
                //                selfScaleToSize(LINEAR, *VolPSF, XSIZE(*VolPSF)*BLOB_SUBSAMPLING,
                //                                YSIZE(*VolPSF)*BLOB_SUBSAMPLING, ZSIZE(*VolPSF));

                blobprint().setXmippOrigin();

                if (XSIZE(blobprint()) < XSIZE(*VolPSF))
                    blobprint().selfWindow(STARTINGY(*VolPSF),
                                           STARTINGX(*VolPSF),
                                           STARTINGY(*VolPSF)+YSIZE(*VolPSF)-1,
                                           STARTINGX(*VolPSF)+XSIZE(*VolPSF)-1);
                else if (XSIZE(blobprint()) > XSIZE(*VolPSF))
                    VolPSF->selfWindow(STARTINGY(blobprint()),
                                       STARTINGX(blobprint()),
                                       STARTINGY(blobprint())+YSIZE(blobprint())-1,
                                       STARTINGX(blobprint())+XSIZE(blobprint())-1);

                // creation of the 3D blobprint from convolution of blobprint and PSF
                ImageOver footprintT;
                footprintT.init(*VolPSF, subsampling);

                convolutionFFT(*VolPSF, blobprint(), footprintT());

                blobprint = footprintT;
            }

            blobprint2()  = blobprint();
            blobprint2() *= blobprint();
            break;
        }
    case (voxels):  sum_on_grid = 1;
        break;
    case (splines): sum_on_grid = sum_spatial_Bspline03_Grid(grid);
        break;
    }

#ifdef DEBUG
    std::cout << "Sum of a basis on the grid=" << sum_on_grid << std::endl;
    std::cout << "D\n" << D << std::endl;
    ImageXmipp save;
    save() = blobprint();
    save.write("footprint.xmp");
#endif
}
#undef DEBUG

// Show --------------------------------------------------------------------
std::ostream & operator << (std::ostream & out, const Basis &basis)
{
    switch (basis.type)
    {
    case Basis::blobs:
        out << " Blobs:         radius=" << basis.blob.radius << " pixels"
        << " alpha=" << basis.blob.alpha
        << " order=" << basis.blob.order << std::endl;
        break;
    case Basis::voxels:
        out << "Voxels\n";
        break;
    case Basis::splines:
        out << "Splines\n";
        break;
    }
    return out;
}

double Basis::maxLength() const
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


// Change to voxels --------------------------------------------------------
void Basis::changeToVoxels(GridVolume &vol_basis, MultidimArray<double> *vol_voxels,
                           int Zdim, int Ydim, int Xdim, int threads ) const
{
    int xdiff, ydiff, zdiff;
    switch (type)
    {
    case blobs:
        blobs2voxels(vol_basis, blob, vol_voxels, D, threads , Zdim, Ydim, Xdim );
        break;
    case voxels:
        *vol_voxels = vol_basis(0)();
        xdiff = (Xdim - XSIZE(*vol_voxels)) / 2;
        ydiff = (Ydim - YSIZE(*vol_voxels)) / 2;
        zdiff = (Zdim - ZSIZE(*vol_voxels)) / 2;
        vol_voxels->selfWindow(
            STARTINGZ(*vol_voxels) - zdiff,
            STARTINGY(*vol_voxels) - ydiff,
            STARTINGX(*vol_voxels) - xdiff,
            (int)(FINISHINGZ(*vol_voxels) + Zdim - ZSIZE(*vol_voxels) - zdiff),
            (int)(FINISHINGY(*vol_voxels) + Ydim - YSIZE(*vol_voxels) - ydiff),
            (int)(FINISHINGX(*vol_voxels) + Xdim - XSIZE(*vol_voxels) - xdiff));
        break;
    case splines:
        spatial_Bspline032voxels(vol_basis, vol_voxels, Zdim, Ydim, Xdim);
        break;
    }
}

// Change from voxels ------------------------------------------------------
void Basis::changeFromVoxels(const MultidimArray<double> &vol_voxels,
                             GridVolume &vol_basis, int grid_type, double grid_relative_size,
                             const MultidimArray<double> *vol_mask,
                             const Matrix2D<double> *D, double R, int threads) const
{
    Grid grid;
    Matrix1D<double> corner1(3), corner2(3);
    double R2 = R * R;
    switch (type)
    {
    case blobs:
        voxels2blobs(&vol_voxels, blob, vol_basis,
                     grid_type, grid_relative_size, 0.05, vol_mask,
                     D, 0.01, 0, R, threads);
        break;
    case voxels:
        VECTOR_R3(corner1, STARTINGX(vol_voxels),
                  STARTINGY(vol_voxels), STARTINGZ(vol_voxels));
        VECTOR_R3(corner2, FINISHINGX(vol_voxels),
                  FINISHINGY(vol_voxels), FINISHINGZ(vol_voxels));
        grid = Create_CC_grid(1, corner1, corner2);
        vol_basis.adapt_to_grid(grid);
        vol_basis(0)() = vol_voxels;
        if (R != -1)
            FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_voxels)
            if (k*k + i*i + j*j > R2)
                vol_basis(0)()(k, i, j) = 0;
        if (vol_mask != NULL)
            FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_voxels)
            if ((*vol_mask)(k, i, j) == 0)
                vol_basis(0)()(k, i, j) = 0;
        break;
    case splines:
        REPORT_ERROR(ERR_NOT_IMPLEMENTED,"");
        /* TODO */
        break;
    }
}

double Basis::valueAt(const Matrix1D<double> & r) const
{
    double module_r;
    switch (type)
    {
    case (blobs):
        {
            module_r = sqrt(XX(r) * XX(r) + YY(r) * YY(r) + ZZ(r) * ZZ(r));
            return blob_val(module_r, blob);
            break;
        }
    case (voxels):
        {
            if (-0.5 <= XX(r) && XX(r) < 0.5 &&
                -0.5 <= YY(r) && YY(r) < 0.5 &&
                -0.5 <= ZZ(r) && ZZ(r) < 0.5)
                return 1.0;
            else
                return 0.0;
            break;
        }
    case (splines):
        {
            if (-2 <= XX(r) && XX(r) < 2 &&
                -2 <= YY(r) && YY(r) < 2 &&
                -2 <= ZZ(r) && ZZ(r) < 2)
                return spatial_Bspline03LUT(r);
            else
                return 0.0;
            break;
        }
    }
    return 0.0;
}

double Basis::projectionAt(const Matrix1D<double> & u, const Matrix1D<double> & r) const
{
    const double p0 = 1.0 / (2 * PIXEL_SUBSAMPLING) - 0.5;
    const double pStep = 1.0 / PIXEL_SUBSAMPLING;
    const double pAvg = 1.0 / (PIXEL_SUBSAMPLING * PIXEL_SUBSAMPLING);
    double module_r, px, py;
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
