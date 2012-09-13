/***************************************************************************
 *
 * Authors:    Sjors Scheres                 (scheres@cnb.uam.es)
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
#include "align2d.h"
#include "fourier_filter.h"

#include <data/funcs.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/args.h>
#include <data/funcs.h>
#ifdef isnan
#include <cmath>
#endif

// Read arguments ==========================================================
void Prog_align2d_prm::read(int argc, char **argv)
{

    fn_sel = getParameter(argc, argv, "-i");
    SF.read(fn_sel);
    // Filename for output reference
    fn_ave = fn_sel.without_extension() + ".xmp";
    fn_ave = fn_ave.insert_before_extension(".ref");
    // Reference image
    fn_ref = getParameter(argc, argv, "-ref", "");
    // Extension for output images (overwrite input if none)
    oext = getParameter(argc, argv, "-oext", "");
    // Write out document file?
    fn_doc = getParameter(argc, argv, "-doc", "");
    // Maximum shift (discard images that shift more in last iteration)
    max_shift = textToFloat(getParameter(argc, argv, "-max_shift", "0"));
    // Maximum rotational change (discard images that rotate more in last iteration)
    max_rot = textToFloat(getParameter(argc, argv, "-max_rot", "0"));
    // Inner and outer radii for rotational correlation
    Ri = textToInteger(getParameter(argc, argv, "-Ri", "0"));
    Ro = textToInteger(getParameter(argc, argv, "-Ro", "0"));
    int dim;
    SF.ImgSize(dim, dim);
    if (Ro == 0) Ro = (int)dim / 2;
    if (Ro <= Ri) REPORT_ERROR(1, "Align2D: Rout should be larger than Rin");
    // Expected resolution and sampling rate (for filtering)
    if (do_filter = checkParameter(argc, argv, "-filter"))
    {
        resol = textToFloat(getParameter(argc, argv, "-filter"));
        sam = textToFloat(getParameter(argc, argv, "-sampling"));
    }
    // Number of iterations
    Niter = textToInteger(getParameter(argc, argv, "-iter", "1"));
    // Only translational/rotational
    do_rot = !checkParameter(argc, argv, "-only_trans");
    do_trans = !checkParameter(argc, argv, "-only_rot");
    do_complete = checkParameter(argc, argv, "-complete");
    psi_interval = textToFloat(getParameter(argc, argv, "-psi_step", "10"));
}

// Show ====================================================================
void Prog_align2d_prm::show()
{
    std::cerr << " Input selfile         : " <<  SF.name() << std::endl;
    if (oext != "")
        std::cerr << " Output extension      : " << oext << std::endl;
    std::cerr << " Number of iterations  : " << Niter << std::endl;
    if (fn_ref != "")
        std::cerr << " Alignment Reference   : " <<  fn_ref << std::endl;
    else
        std::cerr << " Alignment Reference   : piramidal combination of images" << std::endl;

    if (do_filter)
    {
        std::cerr << " Low pass filter [Ang] : " <<  resol << std::endl;
        std::cerr << " Sampling rate   [Ang] : " <<  sam << std::endl;
    }
    if (Ri != 0 || Ro != 0)
        std::cerr << " Inner radius          : " << Ri  << std::endl;
    std::cerr << " Outer radius          : " << Ro  << std::endl;
    if (max_shift != 0)
        std::cerr << " Max. shift last iter. : " << max_shift << std::endl;
    if (max_rot != 0)
        std::cerr << " Max. rotat. last iter.: " << max_rot << std::endl;
    if (fn_doc != "")
        std::cerr << " Output document file  : " << fn_doc << std::endl;
    if (!do_rot)
        std::cerr << "Skip rotational alignment " << std::endl;
    if (!do_trans)
        std::cerr << "Skip translational alignment " << std::endl;
    if (do_complete)
    {
        std::cerr << "Use complete-search alignment with:" << std::endl;
        std::cerr << " Psi interval          : " << psi_interval  << std::endl;
    }
}

// usage ===================================================================
void Prog_align2d_prm::usage()
{
    std::cerr << "Usage:  " << std::endl;
    std::cerr << "  align2d [options]" << std::endl;
    std::cerr << "   -i <selfile>             : Selfile containing images to be aligned \n"
    << " [ -ref <image> ]           : reference image; if none: piramidal combination of subset of images \n"
    << " [ -oext <extension> ]      : For output images & selfile; if none: input will be overwritten \n"
    << " [ -iter <int=1> ]          : Number of iterations to perform \n"
    << " [ -Ri <inner radius=0> ]   : Region between radii Ri and Ro will be considered \n"
    << " [ -Ro <outer radius=dim/2>]     for rotational correlation\n"
    << " [ -filter <resol.> ]       : Fourier-filter images to expected resolution [Ang] \n"
    << " [ -sampling <pix. size> ]  : Sampling rate, i.e. pixel size [Ang]; required for filtering. \n"
    << " [ -max_shift <float> ]     : Discard images that shift more in the last iteration [pix]\n"
    << " [ -max_rot <float> ]       : Discard images that rotate more in the last iteration [deg]\n"
    << " [ -doc <docfile> ]         : write output document file with rotations & translations \n"
    << " [ -only_trans ]            : Skip rotational alignment \n"
    << " [ -only_rot ]              : Skip translational alignment \n"
    << " [ -complete ]              : Use complete-search alignment \n"
    << " [ -psi_step <float=10>]    : Sampling interval to search rotation [deg] \n"

    << std::endl;
}

// Rotational alignment ========================================================
bool Prog_align2d_prm::align_rot(ImageXmippT<double> &img, const Matrix2D<double> &Mref,
                                 const float &max_rot, const float &Rin, const float &Rout, const double &outside)
{

    Matrix2D<double> Mimg, Maux, A;
    Matrix1D<double> corr;
    Matrix2D<int>    mask;
    int nstep;
    int i, i_maxcorr, avewidth;
    double psi_actual, psi_max_coarse, psi_max_fine, sumcorr, psi_coarse_step = 15.;
    float psi;

    mask.resize(img().rowNumber(), img().colNumber());
    mask.setXmippOrigin();
    if (Rout <= Rin)  REPORT_ERROR(1, "Align2d_rot: Rout <= Rin");
    BinaryCrownMask(mask, Rin, Rout, INNER_MASK);

    Mimg.resize(img());
    Mimg.setXmippOrigin();
    A = img.get_transformation_matrix();
    applyGeometry(Mimg, A, img(), IS_INV, DONT_WRAP, outside);
    Maux.resize(Mimg);
    Maux.setXmippOrigin();

    // Optimize correlation in coarse steps
    nstep = (int)(360 / psi_coarse_step);
    corr.resize(nstep);
    for (i = 0; i < nstep; i++)
    {
        psi_actual = (double)i * psi_coarse_step;
        Mimg.rotate(psi_actual, Maux, DONT_WRAP);
        corr(i) = correlation(Mref, Maux, &mask);
    }
    corr.maxIndex(i_maxcorr);
    psi_max_coarse = (double)i_maxcorr * psi_coarse_step;

    // Optimize correlation in fine steps
    nstep = (int)(2 * psi_coarse_step - 1);
    corr.resize(nstep);
    for (i = 0; i < nstep; i++)
    {
        psi_actual = psi_max_coarse - psi_coarse_step + 1 + (double)i;
        Mimg.rotate(psi_actual, Maux, DONT_WRAP);
        corr(i) = correlation(Mref, Maux, &mask);
    }
    corr.maxIndex(i_maxcorr);
    psi_max_fine = psi_max_coarse - psi_coarse_step + (double)i_maxcorr;

    psi = 0.;
    sumcorr = 0.;
    // Weighted average over neighbours
    avewidth = XMIPP_MIN(i_maxcorr, nstep - i_maxcorr);
    if (avewidth > 0)
    {
        for (i = i_maxcorr - avewidth + 1; i < i_maxcorr + avewidth; i++)
        {
            psi_actual = psi_max_coarse - psi_coarse_step + 1 + (double)i;
            psi += corr(i) * psi_actual;
            sumcorr += corr(i);
        }
        psi /= sumcorr;
    }
    else psi = psi_max_fine;

    Mimg.clear();
    Maux.clear();
    mask.clear();

    psi = realWRAP(psi, -180., 180.);
    if ((max_rot < XMIPP_EQUAL_ACCURACY) || (ABS(psi) < max_rot))
    {
        // Store new rotation in the header of the image
        // beware: for untilted images (phi+psi) is rotated, and phi can be non-zero!
        // Add new rotation to psi only!
        psi += img.Psi();
        img.Psi() = realWRAP(psi, 0., 360.);
        return true;
    }
    else return false;
}

// translational alignment =====================================================
bool Prog_align2d_prm::align_trans(ImageXmippT<double> &img, const Matrix2D<double> &Mref, const float &max_shift,
                                   const double &outside)
{

    Matrix2D<double> Maux, Mcorr, A;
    int              dim, imax, jmax, i_actual, j_actual, dim2;
    double           max, xmax, ymax, sumcorr;
    float            xshift, yshift, shift;

    xshift = 0.;
    yshift = 0.;
    dim = img().rowNumber();
    Maux.resize(img());
    Maux.setXmippOrigin();
    Mcorr.resize(img());
    Mcorr.setXmippOrigin();

    // Apply transformation already present in its header
    A = img.get_transformation_matrix();
    applyGeometry(Maux, A, img(), IS_INV, DONT_WRAP, outside);

    // Calculate cross-correlation
    correlation_matrix(Maux, Mref, Mcorr);
    Mcorr.statisticsAdjust(0., 1.);
    Mcorr.maxIndex(imax, jmax);
    max = MAT_ELEM(Mcorr, imax, jmax);

    int              n_max = -1;
    bool             neighbourhood = true;
    while (neighbourhood)
    {
        n_max ++;
        for (int i = -n_max; i <= n_max; i++)
            for (int j = -n_max; j <= n_max; j++)
            {
                i_actual = i + imax;
                j_actual = j + jmax;
                if (i_actual < Mcorr.startingY()  || j_actual < Mcorr.startingX() ||
                    i_actual > Mcorr.finishingY() || j_actual > Mcorr.finishingX())
                    neighbourhood = false;
                else if (max / 1.414 > MAT_ELEM(Mcorr, i_actual, j_actual))
                    neighbourhood = false;
            }
    }
    // We have the neighbourhood => looking for the gravity centre

    xmax = ymax = sumcorr = 0.;
    for (int i = -n_max; i <= n_max; i++)
        for (int j = -n_max; j <= n_max; j++)
        {
            i_actual = i + imax;
            j_actual = j + jmax;
            if (i_actual >= Mcorr.startingY()  && j_actual >= Mcorr.startingX() &&
                i_actual <= Mcorr.finishingY() && j_actual <= Mcorr.finishingX())
            {
                ymax += i_actual * MAT_ELEM(Mcorr, i_actual, j_actual);
                xmax += j_actual * MAT_ELEM(Mcorr, i_actual, j_actual);
                sumcorr += MAT_ELEM(Mcorr, i_actual, j_actual);
            }
        }
    
    //if  (sumcorr  !=  sumcorr)
    xmax /= sumcorr;
    ymax /= sumcorr;
    xshift = (float) - xmax;
    yshift = (float) - ymax;
    //Protect shifts from wrong calculation of correlation maxima
    //I do not know why this happends, but when many-may particles
    //are selected and many iterations performed sooner or later
    //a NAN appears.
    if(sumcorr==0)
    {
        yshift=xshift=0.f;
    }
    //The following code is not portable but 
    // I do not know a better why to detect nans
#ifdef isnan
    if (!isnormal(xmax) || !isnormal(ymax))
    {
        xshift = 0.f;
        yshift = 0.f;
    }
#else
    if((xshift>(float)dim && xshift<(float)dim) ||
       (yshift>(float)dim && yshift<(float)dim))
    {
        yshift=xshift=0.f;
    }   
#endif    
    Maux.clear();
    Mcorr.clear();

    shift = sqrt(xshift * xshift + yshift * yshift);
    if ((max_shift < XMIPP_EQUAL_ACCURACY) || (shift < max_shift))
    {
        // Store shift in the header of the image
        img.Xoff() += xshift * DIRECT_MAT_ELEM(A, 0, 0) + yshift * DIRECT_MAT_ELEM(A, 0, 1);
        img.Yoff() += xshift * DIRECT_MAT_ELEM(A, 1, 0) + yshift * DIRECT_MAT_ELEM(A, 1, 1);
        img.Xoff() = realWRAP(img.Xoff(), (float) - dim / 2., (float)dim / 2.);
        img.Yoff() = realWRAP(img.Yoff(), (float) - dim / 2., (float)dim / 2.);
        return true;
    }
    else return false;

}

// Complete search alignment ========================================================
bool Prog_align2d_prm::align_complete_search(ImageXmippT<double> &img, const Matrix2D<double> &Mref,
        const float &max_shift, const float &max_rot, const float &psi_interval,
        const float &Rin, const float &Rout, const double &outside)
{

    Matrix2D<double> Mimg, Maux, Mcorr, Mref2, A;
    Matrix2D<int>    mask;
    int dim, nstep, imax, jmax;
    double psi_actual, corr, maxcorr, xshift, yshift, shift, psi_max;
    bool OK;

    dim = img().rowNumber();
    mask.resize(dim, dim);
    Mref2 = Mref;
    Maux.resize(Mimg);
    Mimg = img();
    Mcorr.resize(img());

    mask.setXmippOrigin();
    Mref2.setXmippOrigin();
    Maux.setXmippOrigin();
    Mimg.setXmippOrigin();
    Mcorr.setXmippOrigin();

    A = img.get_transformation_matrix();
    applyGeometry(Mimg, A, img(), IS_INV, DONT_WRAP, outside);

    if (Rout <= Rin)  REPORT_ERROR(1, "Align2d: Rout <= Rin");
    BinaryCrownMask(mask, Rin, Rout, INNER_MASK);
    apply_binary_mask(mask, Mref, Mref2);

    // Optimize correlation in coarse steps
    nstep = (int)(360 / psi_interval);
    maxcorr = 0.;
    xshift = 0.;
    yshift = 0.;
    psi_max = 0.;

    for (int i = 0; i < nstep; i++)
    {
        psi_actual = (double)i * psi_interval;
        Mref2.rotate(-psi_actual, Maux, DONT_WRAP);
        correlation_matrix(Mimg, Maux, Mcorr);
        Mcorr.maxIndex(imax, jmax);
        corr = MAT_ELEM(Mcorr, imax, jmax);
        if (corr > maxcorr)
        {
            maxcorr = corr;
            psi_max = psi_actual;
            xshift = -jmax;
            yshift = -imax;
        }
    }

    OK = true;
    psi_max = realWRAP(psi_max, -180., 180.);
    shift = sqrt(xshift * xshift + yshift * yshift);
    if ((max_shift > XMIPP_EQUAL_ACCURACY) && (shift > max_shift)) OK = false;
    if ((max_rot > XMIPP_EQUAL_ACCURACY) && (ABS(psi_max) > max_rot)) OK = false;

    if (OK)
    {
        // Store rotation & translation in the header of the image
        img.Psi() += psi_max;
        img.Psi() = realWRAP(img.Psi(), 0., 360.);
        img.Xoff() += xshift * DIRECT_MAT_ELEM(A, 0, 0) + yshift * DIRECT_MAT_ELEM(A, 0, 1);
        img.Yoff() += xshift * DIRECT_MAT_ELEM(A, 1, 0) + yshift * DIRECT_MAT_ELEM(A, 1, 1);
        img.Xoff() = realWRAP(img.Xoff(), (float) - dim / 2., (float)dim / 2.);
        img.Yoff() = realWRAP(img.Yoff(), (float) - dim / 2., (float)dim / 2.);
        return true;
    }
    else return false;

}

// PsPc piramidal combination of images ========================================
void Prog_align2d_prm::do_pspc()
{

    int               barf, imgno, nlev, n_piram, nlevimgs;
    float             xshift, yshift, psi, zero = 0.;
    Matrix2D<double>  Mref, Maux;
    Matrix2D<int>     mask;

    // Set-up matrices, etc.
    Mref.resize(images[0]());
    Maux.resize(Mref);
    Mref.setXmippOrigin();
    Maux.setXmippOrigin();

    // Calculate average image of non-aligned images to center pspc-reference
    ImageXmipp med, sig;
    double min, max;
    SF.get_statistics(med, sig, min, max, true);
    med().setXmippOrigin();

    // Use piramidal combination of images to construct an initial reference
    nlev = SF.ImgNo();
    int i = 0;
    while (nlev > 1)
    {
        nlev = nlev / 2;
        i++;
    }
    nlev = i;
    n_piram = (int)pow(2., (double)nlev);

    // Copy n_piram to a new temporary array of ImageXmipp
    std::vector<ImageXmipp>  imgpspc;
    for (imgno = 0;imgno < n_piram;imgno++) imgpspc.push_back(images[imgno]);

    std::cerr << "  Piramidal combination of " << n_piram << " images" << std::endl;
    init_progress_bar(n_piram);
    barf = XMIPP_MAX(1, (int)(1 + (n_piram / 60)));

    imgno = 0;
    for (int lev = nlev; lev > 0; lev--)
    {

        nlevimgs = (int)pow(2., (double)lev);
        for (int j = 0; j < nlevimgs / 2; j++)
        {

            applyGeometry(Mref, imgpspc[2*j].get_transformation_matrix(), imgpspc[2*j](), IS_INV, DONT_WRAP);
            Mref.setXmippOrigin();

            if (do_complete)
            {
                align_complete_search(imgpspc[2*j+1], Mref, zero, zero, psi_interval, Ri, Ro);
            }
            else
            {
                /** FIRST **/
                if (do_trans) align_trans(imgpspc[2*j+1], Mref, zero);
                if (do_rot)   align_rot(imgpspc[2*j+1], Mref, zero, Ri, Ro);
                /** SECOND **/
                if (do_trans) align_trans(imgpspc[2*j+1], Mref, zero);
                if (do_rot)   align_rot(imgpspc[2*j+1], Mref, zero, Ri, Ro);
            }

            // Re-calculate average image
            applyGeometry(Maux, imgpspc[2*j+1].get_transformation_matrix(), imgpspc[2*j+1](), IS_INV, DONT_WRAP);
            Maux.setXmippOrigin();
            Mref = (Mref + Maux) / 2;
            imgpspc[j]() = Mref;
            imgpspc[j].Psi() = 0.;
            imgpspc[j].Xoff() = 0.;
            imgpspc[j].Yoff() = 0.;

            imgno++;
            if (imgno % barf == 0) progress_bar(imgno);

        }//loop over images

    }//loop over levels
    progress_bar(n_piram);
    imgpspc.clear();

    // Center pspc reference wrt to average image
    Matrix1D<double> center(2);
    Matrix2D<double> Mcorr;
    Mcorr.resize(Mref);
    int imax, jmax;
    correlation_matrix(Mref, med(), Mcorr);
    Mcorr.maxIndex(imax, jmax);
    XX(center) = -jmax;
    YY(center) = -imax;
    Mref.selfTranslate(center);

    // Write out inter-mediate reference
    Iref() = Mref;
    FileName fn_tmp;
    fn_tmp = fn_sel.without_extension() + ".xmp";
    fn_tmp = fn_tmp.insert_before_extension(".pspc");
    Iref.write(fn_tmp);
    Iref.write(fn_ave);
}

// Alignment of all images by iterative refinement  ========================================
void Prog_align2d_prm::refinement()
{

    int               dim, n_refined, barf;
    float             curr_max_shift, curr_max_rot, xshift, yshift, psi, zero = 0.;
    double            outside, dummy;
    Matrix2D<double>  Mref, Maux, Msum;
    Matrix2D<int>     mask;

    // Set-up matrices, etc.
    dim = Iref().colNumber();
    mask.resize(dim, dim);
    BinaryCircularMask(mask, (int)dim / 2, OUTSIDE_MASK);
    Mref = Iref();
    Mref.setXmippOrigin();
    Maux.resize(Mref);
    Msum.resize(Mref);
    Maux.setXmippOrigin();
    Msum.setXmippOrigin();
    Msum.initZeros();
    barf = XMIPP_MAX(1, (int)(1 + (n_images / 60)));

    n_refined = 0;
    std::cerr << "  Alignment:  iteration " << 1 << " of " << Niter << " (with " << n_images << " images)" << std::endl;
    for (int iter = 0; iter < Niter; iter++)
    {

        if (iter > 0) std::cerr << "  Refinement: iteration " << iter + 1 << " of " << Niter << std::endl;
        if (iter == (Niter - 1))
        {
            curr_max_rot = max_rot;
            curr_max_shift = max_shift;
        }
        else
        {
            curr_max_rot = zero;
            curr_max_shift = zero;
        }

        init_progress_bar(n_images);
        for (int imgno = 0; imgno < n_images; imgno++)
        {

            // Following to compute best value for "outside", assuming no large shifts...
            computeStats_within_binary_mask(mask, images[imgno](), dummy, dummy, outside, dummy);

            if (iter != 0)
            {
                // Subtract current image from the reference
                applyGeometry(Maux, images[imgno].get_transformation_matrix(), images[imgno](), IS_INV, DONT_WRAP, outside);
                Maux.setXmippOrigin();
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mref)
                {
                    dMij(Mref, i, j) *= n_refined;
                    dMij(Mref, i, j) -= dMij(Maux, i, j);
                    dMij(Mref, i, j) /= (n_refined - 1);
                }
            }

            // Align translationally and rotationally
            if (do_complete)
            {
                success[imgno] = align_complete_search(images[imgno], Mref, curr_max_shift,
                                                       curr_max_rot, psi_interval, Ri, Ro, outside);
            }
            else
            {
                if (do_trans) success[imgno] = align_trans(images[imgno], Mref, curr_max_shift);
                if (do_rot && success[imgno]) success[imgno] = align_rot(images[imgno], Mref, curr_max_rot, Ri, Ro, outside);
            }

            if (!success[imgno])
            {
                n_refined--;
            }
            else
            {
                if (iter == 0)
                {
                    // Add refined images to form a new reference
                    n_refined++;
                    applyGeometry(Maux, images[imgno].get_transformation_matrix(), images[imgno](), IS_INV, DONT_WRAP, outside);
                    Maux.setXmippOrigin();
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Msum)
                    {
                        dMij(Msum, i, j) += dMij(Maux, i, j);
                    }
                }
                else
                {
                    // Add refined image to reference again
                    applyGeometry(Maux, images[imgno].get_transformation_matrix(), images[imgno](), IS_INV, DONT_WRAP, outside);
                    Maux.setXmippOrigin();
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mref)
                    {
                        dMij(Mref, i, j) *= (n_refined - 1);
                        dMij(Mref, i, j) += dMij(Maux, i, j);
                        dMij(Mref, i, j) /= n_refined;
                    }
                }
            }
            if (imgno % barf == 0) progress_bar(imgno);

        } // loop over all images
        progress_bar(n_images);

        if (n_images > n_refined) std::cerr << "  Discarded " << n_images - n_refined << " images." << std::endl;
        if (iter == 0) FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Msum)
        {
            dMij(Mref, i, j) = dMij(Msum, i, j) / (n_refined);
        }

        // Write out inter-mediate reference
        Mref.setXmippOrigin();
        Iref() = Mref;
        Iref.write(fn_ave);

    } // loop over iterations
}

// Write out results  ========================================================
void Prog_align2d_prm::calc_correlation(const Matrix2D<double> &Mref, const float &Rin, const float &Rout)
{

    Matrix2D<double> Maux, Mimg;
    Matrix2D<int>    mask;
    int dim = Mref.rowNumber();
    int nstep = (int)(360 / psi_interval);
    Matrix1D<double> ccf;
    double psi_actual;
    FileName fn_tmp;

    ccf.resize(nstep);
    Maux.resize(Mref);
    Mimg.resize(Mref);
    mask.resize(dim, dim);
    Mimg.setXmippOrigin();
    Maux.setXmippOrigin();
    mask.setXmippOrigin();
    if (Rout <= Rin)  REPORT_ERROR(1, "Align2d: Rout <= Rin");
    BinaryCrownMask(mask, Rin, Rout, INNER_MASK);

    for (int imgno = 0; imgno < n_images; imgno++)
    {
        applyGeometry(Mimg, images[imgno].get_transformation_matrix(), images[imgno](), IS_INV, DONT_WRAP);
        Mimg.setXmippOrigin();
        corr[imgno] = correlation_index(Mref, Mimg, &mask);

        if (do_rot)
        {
            for (int i = 0; i < nstep; i++)
            {
                psi_actual = (double)i * psi_interval;
                Mimg.rotate(psi_actual, Maux, DONT_WRAP);
                ccf(i) += correlation_index(Mref, Maux, &mask);
            }
        }
    }
    if (do_rot)
    {
        ccf /= n_images;
        fn_tmp = fn_sel.without_extension() + ".corr";
        if (oext != "") fn_tmp = fn_tmp.insert_before_extension("_" + oext);
        //Output rotation correlation file
        std::ofstream out(fn_tmp.c_str(), std::ios::out);
        out << "# Angle [deg]   Corr.Coeff." << std::endl;
        for (int i = 0; i < nstep; i++)
        {
            out << i*psi_interval << "  " << ccf(i) << std::endl;
        }
        out.close();
    }
}

// Write out results  ========================================================
void Prog_align2d_prm::align2d()
{

    // Read in all images
    double zero = 0.;
    FileName fn_img;
    ImageXmipp Itmp;
    n_images = 0;
    SF.go_beginning();
    while ((!SF.eof()))
    {
        fn_img = SF.NextImg();
        if (fn_img == "") continue;
        Itmp.read(fn_img);
        Itmp().setXmippOrigin();
        images.push_back(Itmp);
        corr.push_back(zero);
        success.push_back(true);
        n_images++;
    }

    // Filter if necessary
    if (do_filter)
    {
        FourierMask fmask;
        fmask.w1 = sam / resol;
        fmask.raised_w = 0.1;
        fmask.FilterShape = RAISED_COSINE;
        fmask.FilterBand = LOWPASS;
        Matrix2D< std::complex<double> > fft;
        FourierTransform(images[0](), fft);
        fmask.generate_mask(fft);
        for (int imgno = 0;imgno < n_images;imgno++)
        {
            FourierTransform(images[imgno](), fft);
            fmask.apply_mask_Fourier(fft);
            InverseFourierTransform(fft, images[imgno]());
        }
    }

    // Get Reference (either from file or from piramidal combination of images)
    if (fn_ref != "")
    {
        Iref.read(fn_ref, false, false, true);
    }
    else do_pspc();

    // Circular mask around reference image
    Matrix2D<int> mask;
    mask.resize(Iref().rowNumber(), Iref().colNumber());
    mask.setXmippOrigin();
    BinaryCrownMask(mask, 0, (int)(Iref().colNumber() / 2), INNER_MASK);
    Iref().setXmippOrigin();
    apply_binary_mask(mask, Iref(), Iref());

    // Do actual alignment & iteratively refine the reference
    refinement();

    // Write out images & selfile
    FileName          fn_out;
    SelFile           SFo;
    SFo.reserve(n_images);
    for (int imgno = 0; imgno < n_images; imgno++)
    {
        fn_img = images[imgno].name();
        if (oext != "") fn_out = fn_img.without_extension() + "." + oext;
        else fn_out = fn_img;
        if (do_filter)
        {
            Itmp.read(fn_img);
            Itmp.Xoff() = images[imgno].Xoff();
            Itmp.Yoff() = images[imgno].Yoff();
            Itmp.Psi() = images[imgno].Psi();
            Itmp.write(fn_out);
        }
        else images[imgno].write(fn_out);
        if (success[imgno]) SFo.insert(fn_out, SelLine::ACTIVE);
        else SFo.insert(fn_out, SelLine::DISCARDED);
    }
    fn_out = fn_sel;
    if (oext != "") fn_out = fn_out.insert_before_extension("_" + oext);
    SFo.write(fn_out);


    // Calculate average and stddev image and write out
//  if(SFo.ImgNo(SelLine::ACTIVE)){
    std::cerr << "Calculating average, correlations and writing out results ..." << std::endl;
    ImageXmipp med, sig;
    double min, max;
    if (SFo.ImgNo(SelLine::ACTIVE))
    {
        SFo.get_statistics(med, sig, min, max, true);
    }
    else
    {
        med().resize(Iref().rowNumber(), Iref().colNumber());
        med().initZeros();
        sig().resize(Iref().rowNumber(), Iref().colNumber());
        sig().initZeros();
    }
    fn_img = fn_sel.without_extension() + ".xmp";
    if (oext != "") fn_img = fn_img.insert_before_extension("_" + oext);
    fn_img = fn_img.insert_before_extension(".med");
    med.weight() = SFo.ImgNo(SelLine::ACTIVE);
    med.write(fn_img);
    fn_img = fn_sel.without_extension() + ".xmp";
    if (oext != "") fn_img = fn_img.insert_before_extension("_" + oext);
    fn_img = fn_img.insert_before_extension(".sig");
    sig.write(fn_img);

    // Calculate correlation wrt average image for document file
    med().setXmippOrigin();
    calc_correlation(med(), Ri, Ro);
    // Write out docfile
    if (fn_doc != "")
    {
        DocFile           DFo;
        DFo.reserve(n_images);
        DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Corr (6)");
        Matrix1D<double>  dataline(6);
        for (int imgno = 0; imgno < n_images; imgno++)
        {
            dataline(0) = images[imgno].Phi();
            dataline(1) = images[imgno].Theta();
            dataline(2) = images[imgno].Psi();
            dataline(3) = images[imgno].Xoff();
            dataline(4) = images[imgno].Yoff();
            dataline(5) = corr[imgno];
            DFo.append_comment(images[imgno].name());
            DFo.append_data_line(dataline);
        }
        DFo.write(fn_doc);
//      }
    }

    images.clear();
    corr.clear();
    success.clear();
}
