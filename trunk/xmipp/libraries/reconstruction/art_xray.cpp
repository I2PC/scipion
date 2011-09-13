/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#include "art_xray.h"


void XrayARTRecons::defineParams(XmippProgram * program, const char* prefix, const char* comment)
{
    char tempLine[256];

    if(prefix == NULL)
        sprintf(tempLine, "  [--xray <psf_param_file>]   : X-ray mode activation");
    else
        sprintf(tempLine,"%s   [--xray <psf_param_file>]   : X-ray mode activation", prefix);
    if (comment != NULL)
        sprintf(tempLine, "%s : %s", tempLine, comment);

    program->addParamsLine(tempLine);


    //    XRayPSF::defineParams(program);
}

void XrayARTRecons::readParams(XmippProgram * program)
{
    ARTReconsBase::readParams(program);
    fnPSF = program->getParam("--xray");
    psf.read(fnPSF);
    //    psf.readParams(program);
}

void XrayARTRecons::produceSideInfo(GridVolume &vol_basis0)
{
    ARTReconsBase::produceSideInfo(vol_basis0);
    psf.calculateParams(artPrm.sampling*1.e-10); // sampling is read in angstrom


}

void XrayARTRecons::print(std::ostream &o) const
{
    o << "X-rays information ------------------------------------------" << std::endl;
    o << "Microscope parameters file: " << fnPSF.c_str() << std::endl;
}

void XrayARTRecons::singleStep(
    GridVolume              &vol_in,          // Input Reconstructed volume
    GridVolume              *vol_out,         // Output Reconstructed volume
    Projection              &theo_proj,       // Projection of the reconstruction
    // It is outside to make it visible
    // just if it needed for any
    // other purpose
    Projection              &read_proj,       // Real projection
    int sym_no,                               // Symmetry matrix index
    Projection              &diff_proj,       // Difference between read and
    // theoretical projection
    Projection              &corr_proj,       // Correcting projection
    Projection              &alig_proj,       // Translation alignement aux proj
    double                  &mean_error,      // Mean error over the pixels
    int                      numIMG,          // number of images in the set
    // in SIRT the correction must
    // be divided by this number
    double                   lambda,          // Lambda to be used
    int                      imagen_no,       // Projection number
    const FileName           &fn_ctf,         // CTF to apply
    const MultidimArray<int> *maskPtr,         // Mask to apply
    bool refine)
{



	 // Only works for blob volumes .............................................
	    if (artPrm.basis.type != Basis::voxels)
	        REPORT_ERROR(ERR_VALUE_INCORRECT,
	                     "This function only works with volxel volumes");

	    // Project structure .......................................................
	    // The correction image is reused in this call to store the normalizing
	    // projection, ie, the projection of an all-1 volume

//	    project_Crystal_Volume(vol_in, artPrm.basis, theo_proj,
//	                           corr_proj, YSIZE(read_proj()), XSIZE(read_proj()),
//	                           read_proj.rot(), read_proj.tilt(), read_proj.psi(), shift,
//	                           aint, bint, *(artPrm.D), *(artPrm.Dinv), this->unit_cell_mask,
//	                           FORWARD, artPrm.eq_mode);
	    double shift_X, shift_Y;

	    //   #define DEBUG_SHIFT
	#ifdef DEBUG_SHIFT

	    Matrix2D<double> A(3, 3);
	    A.initIdentity();
	    dMij(A, 0, 2) =  8;
	    dMij(A, 1, 2) =  -5;
	    std::cout << "A" << A;
	    //move read_proj
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(theo_proj())
	    dMij(theo_proj(), i, j) = dMij(read_proj(), i, j);
	    applyGeometry(IMGMATRIX(read_proj), A, IMGMATRIX(theo_proj), IS_NOT_INV, WRAP);
	#endif
	#undef DEBUG_SHIFT

	    if (artPrm.ref_trans_after != -1    &&
	        imagen_no > artPrm.ref_trans_after && imagen_no != 0)
	    {
//	        calculate_and_find_correlation_max_proj(read_proj, theo_proj,
//	                                                alig_proj,
//	                                                shift_X, shift_Y,
//	                                                artPrm.ref_trans_step,
//	                                                artPrm.ref_trans_after,
//	                                                imagen_no);

	        // Apply correction
	        Matrix2D<double> Correction(3, 3);
	        alig_proj().resize(read_proj());
	        Correction.initIdentity();

	        dMij(Correction, 0, 2) =  - shift_X;
	        dMij(Correction, 1, 2) =  - shift_Y;
	        //copy theo_proj to a temporal matrix
	        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(theo_proj())
	        dAij(alig_proj(), i, j) = dAij(read_proj(), i, j);
	        applyGeometry(LINEAR, IMGMATRIX(alig_proj), IMGMATRIX(read_proj), Correction,
	                      IS_NOT_INV, WRAP);
	    }
	    // Now compute differences .................................................
	    double applied_lambda = lambda / numIMG; // In ART mode, numIMG=1
	    mean_error = 0;
	    diff_proj().resize(read_proj());

	    FOR_ALL_ELEMENTS_IN_ARRAY2D(IMGMATRIX(read_proj))
	    {
	        // Compute difference image and error
	        IMGPIXEL(diff_proj, i, j) = IMGPIXEL(read_proj, i, j) - IMGPIXEL(theo_proj, i, j);
	        mean_error += IMGPIXEL(diff_proj, i, j) * IMGPIXEL(diff_proj, i, j);

	        // Compute the correction image
	        if (ABS(IMGPIXEL(corr_proj, i, j)) < 1)
	            IMGPIXEL(corr_proj, i, j) = SGN(IMGPIXEL(corr_proj, i, j));
	        IMGPIXEL(corr_proj, i, j) =
	            applied_lambda * IMGPIXEL(diff_proj, i, j) / IMGPIXEL(corr_proj, i, j);
	    }
	    mean_error /= XSIZE(diff_proj()) * YSIZE(diff_proj());

	    // Backprojection of correction plane ......................................
//	    project_Crystal_Volume(*vol_out, artPrm.basis, theo_proj,
//	                           corr_proj, YSIZE(read_proj()), XSIZE(read_proj()),
//	                           read_proj.rot(), read_proj.tilt(), read_proj.psi(), shift,
//	                           aint, bint, *(artPrm.D), *(artPrm.Dinv), this->unit_cell_mask,
//	                           BACKWARD, artPrm.eq_mode);
}

