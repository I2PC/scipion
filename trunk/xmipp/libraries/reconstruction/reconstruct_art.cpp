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

#include "reconstruct_art.h"
#include "denoise.h"
#include "fourier_filter.h"

#include <data/wavelet.h>

/* ------------------------------------------------------------------------- */
/* Plain ART Parameters                                                      */
/* ------------------------------------------------------------------------- */
/* Produce Side information ------------------------------------------------ */
void Plain_ART_Parameters::produce_Side_Info(const Basic_ART_Parameters &prm,
        GridVolume &vol_basis0)
{}

/* std::cout -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &o, const Plain_ART_Parameters &eprm)
{
    return o;
}

/* ------------------------------------------------------------------------- */
/* Update residual vector for WLS                                            */
/* ------------------------------------------------------------------------- */
void update_residual_vector(Basic_ART_Parameters &prm, GridVolume &vol_basis,
                            double &kappa, double &pow_residual_vol, double &pow_residual_imgs)
{
    GridVolume       residual_vol;
    Projection       read_proj, dummy_proj, new_proj;
    FileName         fn_resi, fn_tmp;
    double           sqrtweight, dim2, norma, normb, apply_kappa;
    ImageOver        *footprint = (ImageOver *) & prm.basis.blobprint;
    ImageOver        *footprint2 = (ImageOver *) & prm.basis.blobprint2;
    Matrix2D<double> *A = NULL;
    std::vector<Matrix2D<double> > newres_imgs;
    Matrix2D<int>    mask;

    residual_vol.resize(vol_basis);
    residual_vol.initZeros();

    // Calculate volume from all backprojected residual images
    std::cerr << "Backprojection of residual images " << std::endl;
    if (!(prm.tell&TELL_SHOW_ERROR)) init_progress_bar(prm.numIMG);

    for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++)
    {
	// Check whether to kill job
	exit_if_not_exists(prm.fn_control);
	
        // backprojection of the weighted residual image
        sqrtweight = sqrt(prm.residual_imgs[iact_proj].weight() / prm.sum_weight);

        read_proj = prm.residual_imgs[iact_proj];
        read_proj() *= sqrtweight;
        dummy_proj().resize(read_proj());
        dummy_proj.set_angles(prm.IMG_Inf[iact_proj].rot,
                              prm.IMG_Inf[iact_proj].tilt,
                              prm.IMG_Inf[iact_proj].psi);

        project_Volume(residual_vol, prm.basis, dummy_proj,
                       read_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                       prm.IMG_Inf[iact_proj].rot,
                       prm.IMG_Inf[iact_proj].tilt,
                       prm.IMG_Inf[iact_proj].psi, BACKWARD, prm.eq_mode,
                       prm.GVNeq, NULL, NULL, prm.ray_length, prm.threads);

        if (!(prm.tell&TELL_SHOW_ERROR))
            if (iact_proj % XMIPP_MAX(1, prm.numIMG / 60) == 0) progress_bar(iact_proj);
    }
    if (!(prm.tell&TELL_SHOW_ERROR)) progress_bar(prm.numIMG);

    // Convert to voxels: solely for output of power of residual volume
    VolumeXmipp      residual_vox;
    int Xoutput_volume_size = (prm.Xoutput_volume_size == 0) ?
                              prm.projXdim : prm.Xoutput_volume_size;
    int Youtput_volume_size = (prm.Youtput_volume_size == 0) ?
                              prm.projYdim : prm.Youtput_volume_size;
    int Zoutput_volume_size = (prm.Zoutput_volume_size == 0) ?
                              prm.projXdim : prm.Zoutput_volume_size;
    prm.basis.changeToVoxels(residual_vol, &(residual_vox()),
                             Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
    pow_residual_vol = residual_vox().sum2() / MULTIDIM_SIZE(residual_vox());
    residual_vox.clear();

    std::cerr << "Projection of residual volume; kappa = " << kappa << std::endl;
    if (!(prm.tell&TELL_SHOW_ERROR)) init_progress_bar(prm.numIMG);

    // Now that we have the residual volume: project in all directions
    pow_residual_imgs = 0.;
    new_proj().resize(read_proj());
    mask.resize(read_proj());
    BinaryCircularMask(mask, YSIZE(read_proj()) / 2, INNER_MASK);

    dim2 = (double)YSIZE(read_proj()) * XSIZE(read_proj());
    for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++)
    {
	// Check whether to kill job
	exit_if_not_exists(prm.fn_control);

        project_Volume(residual_vol, prm.basis, new_proj,
                       dummy_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                       prm.IMG_Inf[iact_proj].rot,
                       prm.IMG_Inf[iact_proj].tilt,
                       prm.IMG_Inf[iact_proj].psi, FORWARD, prm.eq_mode,
                       prm.GVNeq, A, NULL, prm.ray_length, prm.threads);

        sqrtweight = sqrt(prm.residual_imgs[iact_proj].weight() / prm.sum_weight);

        // Next lines like normalization in [EHL] (2.18)?
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(new_proj())
        {
            dMij(dummy_proj(), i, j) = XMIPP_MAX(1., dMij(dummy_proj(), i, j)); // to avoid division by zero
            dMij(new_proj(), i, j) /= dMij(dummy_proj(), i, j);
        }
        new_proj() *= sqrtweight * kappa;

        /*
        fn_tmp="residual_"+integerToString(iact_proj);
        dummy_proj()=1000*prm.residual_imgs[iact_proj]();
        dummy_proj.write(fn_tmp+".old");
        */

        prm.residual_imgs[iact_proj]() -= new_proj();
        pow_residual_imgs += prm.residual_imgs[iact_proj]().sum2();

        // Mask out edges of the images
        apply_binary_mask(mask, prm.residual_imgs[iact_proj](), prm.residual_imgs[iact_proj](), 0.);

        /*
        dummy_proj()=1000*new_proj();
        dummy_proj.write(fn_tmp+".change");
        dummy_proj()=1000*prm.residual_imgs[iact_proj]();
        dummy_proj.write(fn_tmp+".new");
        */

        if (!(prm.tell&TELL_SHOW_ERROR))
            if (iact_proj % XMIPP_MAX(1, prm.numIMG / 60) == 0) progress_bar(iact_proj);
    }

    pow_residual_imgs /= dim2;
    newres_imgs.clear();

    if (!(prm.tell&TELL_SHOW_ERROR)) progress_bar(prm.numIMG);

}

/* ------------------------------------------------------------------------- */
/* ART Single step                                                           */
/* ------------------------------------------------------------------------- */
void ART_single_step(
    GridVolume              &vol_in,         // Input Reconstructed volume
    GridVolume              *vol_out,        // Output Reconstructed volume
    Basic_ART_Parameters    &prm,            // blob, lambda
    Plain_ART_Parameters    &eprm,           // In this case, nothing
    Projection              &theo_proj,      // Projection of the reconstruction
                                             // It is outside to make it visible
                                             // just if it needed for any
                                             // other purpose
    Projection             &read_proj,       // Real projection
    int sym_no,                              // Symmetry matrix index
    Projection             &diff_proj,       // Difference between read and
                                             // theoretical projection
    Projection             &corr_proj,       // Correcting projection
    Projection             &alig_proj,       // Translation alignement aux proj
    double                 &mean_error,      // Mean error over the pixels
    int                     numIMG,          // number of images in the set
                                             // in SIRT the correction must
                                             // be divided by this number
    double                  lambda,          // Lambda to be used
    int                     act_proj,        // Projection number
    const FileName         &fn_ctf,          // CTF to apply
    const Matrix2D<int>    *maskPtr)         // Mask to apply
{
// Prepare to work with CTF ................................................
    FourierMask ctf;
    ImageOver *footprint = (ImageOver *) & prm.basis.blobprint;
    ImageOver *footprint2 = (ImageOver *) & prm.basis.blobprint2;
    bool remove_footprints = false;
    double weight, sqrtweight;

    if (fn_ctf != "" && !prm.unmatched)
    {
        if (prm.basis.type != Basis::blobs)
            REPORT_ERROR(1, "ART_single_step: This way of correcting for the CTF "
                         "only works with blobs");
        // It is a description of the CTF
        ctf.FilterShape = ctf.FilterBand = CTF;
        ctf.ctf.enable_CTFnoise = false;
        ctf.ctf.read(fn_ctf);
        ctf.ctf.Tm /= BLOB_SUBSAMPLING;
        ctf.ctf.Produce_Side_Info();

        // Create new footprints
        footprint = new ImageOver;
        footprint2 = new ImageOver;
        remove_footprints = true;

        // Enlarge footprint, bigger than necessary to avoid
        // aliasing
        *footprint = prm.basis.blobprint;
        (*footprint)().setXmippOrigin();
        double blob_radius = prm.basis.blob.radius;
        int finalsize = 2 * CEIL(30 + blob_radius) + 1;
        footprint->window(
            FIRST_XMIPP_INDEX(finalsize), FIRST_XMIPP_INDEX(finalsize),
            LAST_XMIPP_INDEX(finalsize), LAST_XMIPP_INDEX(finalsize));

        // Generate mask to the size of the footprint, correct phase
        // and apply CTF
        ctf.generate_mask((*footprint)());
        ctf.correct_phase();
        ctf.apply_mask_Space((*footprint)());

        // Remove unnecessary regions
        finalsize = 2 * CEIL(15 + blob_radius) + 1;
        footprint->window(
            FIRST_XMIPP_INDEX(finalsize), FIRST_XMIPP_INDEX(finalsize),
            LAST_XMIPP_INDEX(finalsize), LAST_XMIPP_INDEX(finalsize));
#ifdef DEBUG
        ImageXmipp save;
        save() = (*footprint)();
        save.write("PPPfootprint.xmp");
#endif

        // Create footprint2
        *footprint2 = *footprint;
        (*footprint2)() *= (*footprint2)();
    }

// Project structure .......................................................
    // The correction image is reused in this call to store the normalising
    // projection, ie, the projection of an all-1 volume
    Matrix2D<double> *A = NULL;
    if (prm.print_system_matrix) A = new Matrix2D<double>;
    corr_proj().initZeros();
    project_Volume(vol_in, prm.basis, theo_proj,
                   corr_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                   read_proj.rot(), read_proj.tilt(), read_proj.psi(), FORWARD, prm.eq_mode,
                   prm.GVNeq, A, maskPtr, prm.ray_length, prm.threads);

    if (fn_ctf != "" && prm.unmatched)
    {
        ctf.generate_mask(theo_proj());
        ctf.apply_mask_Space(theo_proj());
    }

    // Check whether to kill job
    exit_if_not_exists(prm.fn_control);

    // Print system matrix
    if (prm.print_system_matrix)
    {
        std::cout << "Equation system (Ax=b) ----------------------\n";
        std::cout << "Size: ";
        A->printShape();
        std::cout << std::endl;
        for (int i = 0; i < YSIZE(*A); i++)
        {
            bool null_row = true;
            for (int j = 0; j < YSIZE(*A); j++)
                if (DIRECT_MAT_ELEM(*A, i, j) != 0)
                {
                    null_row = false;
                    break;
                }
            if (!null_row)
            {
                std::cout << "pixel=" << integerToString(i, 3) << " --> "
                    << DIRECT_MULTIDIM_ELEM(read_proj(), i) << " = ";
                for (int j = 0; j < XSIZE(*A); j++)
                    std::cout << DIRECT_MAT_ELEM(*A, i, j) << " ";
                std::cout << std::endl;
            }
        }
        std::cout << "---------------------------------------------\n";
        delete A;
    }

    // Now compute differences .................................................
    double applied_lambda = lambda / numIMG; // In ART mode, numIMG=1

    mean_error = 0;
    diff_proj().resize(read_proj());

    // Weighted least-squares ART for Maximum-Likelihood refinement
    if (prm.WLS)
    {
        weight = read_proj.weight() / prm.sum_weight;
        sqrtweight = sqrt(weight);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(read_proj))
        {
            // Compute difference image and error
            IMGPIXEL(diff_proj, i, j) = IMGPIXEL(read_proj, i, j) - IMGPIXEL(theo_proj, i, j);
            mean_error += IMGPIXEL(diff_proj, i, j) * IMGPIXEL(diff_proj, i, j);

            // Subtract the residual image (stored in alig_proj!)
            IMGPIXEL(diff_proj, i, j) = sqrtweight * IMGPIXEL(diff_proj, i, j) - IMGPIXEL(alig_proj, i, j);

            // Calculate the correction and the updated residual images
            IMGPIXEL(corr_proj, i, j) =
                applied_lambda * IMGPIXEL(diff_proj, i, j) / (weight * IMGPIXEL(corr_proj, i, j) + 1.);
            IMGPIXEL(alig_proj, i, j) += IMGPIXEL(corr_proj, i, j);
            IMGPIXEL(corr_proj, i, j) *= sqrtweight;

        }
        mean_error /= XSIZE(diff_proj()) * YSIZE(diff_proj());
        mean_error *= weight;

    }
    else
    {
        long int Nmean=0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(read_proj))
        {
            if (maskPtr!=NULL)
                if ((*maskPtr)(i,j)<0.5) continue;
            // Compute difference image and error
            IMGPIXEL(diff_proj, i, j) = IMGPIXEL(read_proj, i, j) - IMGPIXEL(theo_proj, i, j);
            mean_error += IMGPIXEL(diff_proj, i, j) * IMGPIXEL(diff_proj, i, j);
            Nmean++;

            // Compute the correction image
            IMGPIXEL(corr_proj, i, j) = XMIPP_MAX(IMGPIXEL(corr_proj, i, j), 1);
            IMGPIXEL(corr_proj, i, j) =
                applied_lambda * IMGPIXEL(diff_proj, i, j) / IMGPIXEL(corr_proj, i, j);
        }
        mean_error /= Nmean;
    }

    // Backprojection of correction plane ......................................
    project_Volume(*vol_out, prm.basis, theo_proj,
                   corr_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                   read_proj.rot(), read_proj.tilt(), read_proj.psi(), BACKWARD, prm.eq_mode,
                   prm.GVNeq, NULL, maskPtr, prm.ray_length, prm.threads);

    // Remove footprints if necessary
    if (remove_footprints)
    {
        delete footprint;
        delete footprint2;
    }
}
#undef blob

/* Finish iterations ------------------------------------------------------- */
void finish_ART_iterations(const Basic_ART_Parameters &prm,
                           const Plain_ART_Parameters &eprm, GridVolume &vol_blobs)
{ }

/* Apply_symmetry ------------------------------------------------------- */
void apply_symmetry(GridVolume &vol_in, GridVolume *vol_out,
                    const Plain_ART_Parameters &eprm, int grid_type)
{
    REPORT_ERROR(1, "apply_symmetry: Function not implemented for single particles");
}
