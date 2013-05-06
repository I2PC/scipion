/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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

#include "reconstruct_wbp.h"
#include <data/metadata_extension.h>

ProgRecWbp::ProgRecWbp()
{
    iter = NULL;
}

ProgRecWbp::~ProgRecWbp()
{
    delete iter;
}

// Read arguments ==========================================================
void ProgRecWbp::readParams()
{
    fn_sel = getParam("-i");
    fn_out = getParam("-o");
    fn_sym = getParam("--sym");
    threshold = getDoubleParam("--threshold");
    diameter = 2 * getIntParam("--radius");

    sampling = getDoubleParam("--filsam");
    do_all_matrices = checkParam("--use_each_image");
    do_weights = checkParam("--weight");
}

// Show ====================================================================
void ProgRecWbp::show()
{
    if (verbose > 0)
    {
        // To screen
        std::cerr
        << " ================================================================="
        << std::endl;
        std::cerr << " Weighted-back projection (arbitrary geometry) "
        << std::endl;
        std::cerr
        << " ================================================================="
        << std::endl;
        std::cerr << " Input selfile             : " << fn_sel << std::endl;
        std::cerr << " Output volume             : " << fn_out << std::endl;
        if (diameter > 0)
            std::cerr << " Reconstruction radius     : " << diameter / 2
            << std::endl;
        std::cerr << " Relative filter threshold : " << threshold << std::endl;
        if (fn_sym != "")
            std::cerr << " Symmetry file:            : " << fn_sym << std::endl;
        if (do_all_matrices)
            std::cerr
            << " --> Use all projection directions in arbitrary geometry filter"
            << std::endl;
        else
            std::cerr << " --> Use sampled directions for filter, sampling = "
            << sampling << std::endl;
        if (do_weights)
            std::cerr << " --> Use weights stored in the image headers"
            << std::endl;
        std::cerr
        << " -----------------------------------------------------------------"
        << std::endl;
    }
}

// Usage ====================================================================
void ProgRecWbp::defineParams()
{
    // To screen
    addUsageLine(
        "Generate 3D reconstruction from projections using the Weighted BackProjection algorithm.");
    addUsageLine(
        "+This program allows you to generate 3D reconstructions from projections ");
    addUsageLine(
        "+using the Weighted BackProjection algorithm with arbitrary geometry as ");
    addUsageLine(
        "+described by Radermacher, M. (1992) Weighted back-projection methods. ");
    addUsageLine("+Electron Tomography, ed. by J. Frank, Plenum Press.");
    addSeeAlsoLine(
        "angular_projection_matching, angular_discrete_assign, angular_continuous_assign");
    addParamsLine(
        "   -i <input_selfile>          : selection file with input images and Euler angles");
    addParamsLine(
        " [ -o <name=\"wbp.vol\"> ]     : filename for output volume ");
    addParamsLine(
        " [ --doc <docfile=\"\"> ]      : Ignore headers and get angles from this docfile ");
    addParamsLine(
        " [ --radius <int=-1> ]         : Reconstruction radius. int=-1 means radius=dim/2 ");
    addParamsLine(
        "                               : The volume will be zero outside this radius");
    addParamsLine(" [ --sym <sym=\"\"> ]          : Enforce symmetry ");
    addParamsLine(
        "                               :+A symmetry file or point-group description.");
    addParamsLine(
        "                               :+Valid point-group descriptions are: C1, Ci, Cs, Cn ");
    addParamsLine(
        "                               :+(from here on n must be an integer number with no ");
    addParamsLine(
        "                               :+more than 2 digits) Cnv, Cnh, Sn, Dn, Dnv, Dnh, T, ");
    addParamsLine(
        "                               :+Td, Th, O, Oh I, I1, I2, I3, I4, I5, Ih. For a full ");
    addParamsLine(
        "                               :+description of symmetries look at http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry");
    addParamsLine(
        " [ --threshold+ <float=0.005> ] : Lower (relative) threshold for filter values ");
    addParamsLine(
        "                               :+This is to avoid divisions by zero and consequent ");
    addParamsLine(
        "                               :+enhancement of noise. The threshold is given as a relative ");
    addParamsLine(
        "                               :+value with respect to the total number of images. The absolute");
    addParamsLine(
        "                               :+threshold will be calculated internally as the relative threshold");
    addParamsLine(
        "                               :+multiplied by the total number of (symmetry generated) projections.");
    addParamsLine(
        " [ --filsam+ <float=5>]        : Angular sampling rate for geometry filter ");
    addParamsLine(
        "                               :+Instead of summing over all experimental images to calculate ");
    addParamsLine(
        "                               :+the arbitrary geometry filter, we bin all images in representative ");
    addParamsLine(
        "                               :+projection directions, which are sampled every filsam degrees.");
    addParamsLine(
        " [ --use_each_image+]          : Use each image instead of sampled representatives for filter ");
    addParamsLine(
        "                               :+If this option is given, all experimental images will be used ");
    addParamsLine(
        "                               :+in the summation to calculate the arbitrary geometry filter. For ");
    addParamsLine(
        "                               :+large datasets this may be considerably slower than the default ");
    addParamsLine(
        "                               :+option of using representative projection directions");
    addParamsLine(
        " [ --weight]                   : Use weights stored in image headers or the input metadata");
    addExampleLine("xmipp_reconstruct_wbp -i images.sel -o reconstruction.vol");
}

void ProgRecWbp::run()
{
    show();
    produceSideInfo();
    apply2DFilterArbitraryGeometry();
    finishProcessing();
}

void ProgRecWbp::finishProcessing()
{
    free(mat_g);
    free(mat_f);
    if (verbose > 0)
        std::cerr << "Fourier pixels for which the threshold was not reached: "
        << (float) (count_thr * 100.) / (SF.size() * dim * dim) << " %"
        << std::endl;
    reconstructedVolume.write(fn_out);
}

void ProgRecWbp::setIO(const FileName &fnIn, const FileName &fnOut)
{
    fn_sel = fnIn;
    fn_out = fnOut;
}

void ProgRecWbp::produceSideInfo()
{
    // Read-in stuff
    SF.read(fn_sel);
    if (do_weights)
    {
        SF.removeObjects(MDValueEQ(MDL_WEIGHT, 0.0));
        if (SF.size() == 0)
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "there is no input file with weight!=0");
    }

    size_t zdim, ndim;
    getImageSize(SF, dim, dim, zdim, ndim);
    if (fn_sym != "")
        SL.readSymmetryFile(fn_sym);
    if (diameter <= 0)
        diameter = dim;

    // Fill arrays of transformation matrices
    if (do_all_matrices)
        getAllMatrices(SF);
    else
        getSampledMatrices(SF);

    time_bar_size = SF.size();
    time_bar_step = CEIL((double)time_bar_size / 60.0);
    time_bar_done = 0;
}

void ProgRecWbp::getAnglesForImage(size_t id, double &rot, double &tilt,
                                   double &psi, double &xoff, double &yoff, double &flip, double &weight)
{
    if (id != BAD_OBJID)
    {
        SF.getValue(MDL_ANGLE_ROT, rot, id);
        SF.getValue(MDL_ANGLE_TILT, tilt, id);
        SF.getValue(MDL_ANGLE_PSI, psi, id);
        SF.getValue(MDL_SHIFT_X, xoff, id);
        SF.getValue(MDL_SHIFT_Y, yoff, id);
        flip = 0;
        SF.getValue(MDL_FLIP, flip, id);
        weight = 1;
        SF.getValue(MDL_WEIGHT, weight, id);
    }
}

void ProgRecWbp::getSampledMatrices(MetaData &SF)
{
    FileName fn_tmp;
    Matrix2D<double> A(3, 3);
    Matrix2D<double> L(4, 4), R(4, 4);
    double newrot, newtilt, newpsi, rot, tilt, dum, weight, totimgs = 0.;
    std::vector<double> count_imgs;

    if (verbose > 0)
        std::cerr << "--> Sampling the filter ..." << std::endl;

    // Create an (asymmetric part of an) even projection direction distribution
    std::vector<double> rotList, tiltList;
    make_even_distribution(rotList, tiltList, sampling, SL, true);
    size_t NN = rotList.size();
    count_imgs.resize(NN);
    // Each experimental image contributes to the nearest of these directions
    FileName fn_img;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        SF.getValue(MDL_IMAGE, fn_img, __iter.objId);
        getAnglesForImage(__iter.objId, rot, tilt, dum, dum, dum, dum, weight);
        int idx = find_nearest_direction(rot, tilt, rotList, tiltList, SL, L,
                                         R);
        if (do_weights)
            count_imgs[idx] += weight;
        else
            count_imgs[idx] += 1.;
    }

    // Now calculate transformation matrices for all representative directions
    no_mats = 0;
    int SL_SymsNo = SL.symsNo();
    int SL_SymsNo_1 = SL_SymsNo + 1;
    for (size_t i = 0; i < NN; i++)
        if (count_imgs[i] > 0.)
            no_mats += SL_SymsNo_1;
    mat_g = (WBPInfo*) malloc(no_mats * sizeof(WBPInfo));

    no_mats = 0;
    for (size_t i = 0; i < NN; i++)
    {
        int count_i = (int)count_imgs[i];
        if (count_i > 0)
        {
            Euler_angles2matrix(rotList[i], -tiltList[i], 0.0, A);
            mat_g[no_mats].x = MAT_ELEM(A, 2, 0);
            mat_g[no_mats].y = MAT_ELEM(A, 2, 1);
            mat_g[no_mats].z = MAT_ELEM(A, 2, 2);
            mat_g[no_mats].count = count_i;
            totimgs += count_i;
            no_mats++;

            // Expand symmetric directions
            for (int j = 0; j < SL_SymsNo; j++)
            {
                SL.getMatrices(j, L, R, false);
                Euler_apply_transf(L, R, rot, tilt, 0., newrot, newtilt, newpsi);
                Euler_angles2matrix(newrot, -newtilt, newpsi, A);
                mat_g[no_mats].x = MAT_ELEM(A, 2, 0);
                mat_g[no_mats].y = MAT_ELEM(A, 2, 1);
                mat_g[no_mats].z = MAT_ELEM(A, 2, 2);
                mat_g[no_mats].count = count_i;
                totimgs += count_i;
                no_mats++;
            }
        }
    }

    // Adjust relative threshold
    threshold *= totimgs;
}

// Fill array with transformation matrices needed for arbitrary geometry filter
void ProgRecWbp::getAllMatrices(MetaData &SF)
{
    Matrix2D<double> A(3, 3);
    Matrix2D<double> L(4, 4), R(4, 4);
    double rot, tilt, psi, weight, dum, newrot, newtilt, newpsi, totimgs = 0.;
    int NN;

    no_mats = 0;

    NN = SF.size();
    NN *= (SL.symsNo() + 1);
    mat_g = (WBPInfo*) malloc(NN * sizeof(WBPInfo));
    FileName fn_img;

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        SF.getValue(MDL_IMAGE, fn_img, __iter.objId);
        getAnglesForImage(__iter.objId, rot, tilt, psi, dum, dum, dum, weight);
        Euler_angles2matrix(rot, -tilt, psi, A);
        mat_g[no_mats].x = MAT_ELEM(A, 2, 0);
        mat_g[no_mats].y = MAT_ELEM(A, 2, 1);
        mat_g[no_mats].z = MAT_ELEM(A, 2, 2);
        if (do_weights)
            mat_g[no_mats].count = weight;
        else
            mat_g[no_mats].count = 1.;
        totimgs += mat_g[no_mats].count;
        no_mats++;
        // Also add symmetry-related projection directions
        for (int i = 0; i < SL.symsNo(); i++)
        {
            SL.getMatrices(i, L, R);
            L.resize(3, 3);
            R.resize(3, 3);
            Euler_apply_transf(L, R, rot, -tilt, psi, newrot, newtilt, newpsi);
            Euler_angles2matrix(newrot, newtilt, newpsi, A);
            mat_g[no_mats].x = MAT_ELEM(A, 2, 0);
            mat_g[no_mats].y = MAT_ELEM(A, 2, 1);
            mat_g[no_mats].z = MAT_ELEM(A, 2, 2);
            if (do_weights)
                mat_g[no_mats].count = weight;
            else
                mat_g[no_mats].count = 1.;
            totimgs += mat_g[no_mats].count;
            no_mats++;
        }
    }

    // Adjust relative threshold
    threshold *= totimgs;
}

// Simple backprojection of a single image
void ProgRecWbp::simpleBackprojection(Projection &img,
                                      MultidimArray<double> &vol, int diameter)
{
	//this should be int not size_t ROB
    int i, j, k, l, m;
    Matrix2D<double> A(3, 3);
    double dim2, x, y, z, xp, yp;
    double value1, value2, scalex, scaley, value;
    double radius2, x2, y2, z2, z2_plus_y2;

    // Use minus-tilt, because code copied from OldXmipp
    Euler_angles2matrix(img.rot(), -img.tilt(), img.psi(), A);
    A = A.inv();

    radius2 = diameter / 2.;
    radius2 = radius2 * radius2;
    dim2 = dim / 2;
    double a00 = MAT_ELEM(A,0,0);
    double a01 = MAT_ELEM(A,0,1);
    double a10 = MAT_ELEM(A,1,0);
    double a11 = MAT_ELEM(A,1,1);
    double a20 = MAT_ELEM(A,2,0);
    double a21 = MAT_ELEM(A,2,1);

    double dim1 = dim - 1;
    const MultidimArray<double> mImg = img();
    int idim;
    idim = dim;//cast to int from size_t
    for (i = 0; i < idim; i++)
    {
        z = -i + dim2; /*** Z points upwards ***/
        z2 = z * z;
        if (z2 > radius2)
            continue;
        double xpz = z * a20 + dim2;
        double ypz = z * a21 + dim2;
        for (j = 0; j < idim; j++)
        {
            y = j - dim2;
            y2 = y * y;
            z2_plus_y2 = z2 + y2;
            if (z2_plus_y2 > radius2)
                continue;
            x = 0 - dim2; /***** X for k == 0 *****/
            xp = x * a00 + y * a10 + xpz;
            yp = x * a01 + y * a11 + ypz;
            if (yp >= dim1 || yp < 0.0)
                continue;
            l = (int) yp;
            scaley = yp - l;
            double scale1y = 1. - scaley;
            for (k = 0; k < idim; k++, xp += a00, yp += a01, x++)
            {
                x2 = x * x;
                if (x2 + z2_plus_y2 > radius2)
                    continue;
                if (xp >= dim1 || xp < 0.0)
                    continue;

                /**** interpolation ****/
                m = (int) xp;
                scalex = xp - m;
                double scale1x = 1. - scalex;
                value1 = scalex * dAij(mImg, l, m + 1)
                         + scale1x * dAij(mImg, l, m);
                value2 = scalex * dAij(mImg, l + 1, m + 1)
                         + scale1x * dAij(mImg, l + 1, m);
                value = scaley * value2 + scale1y * value1;
                dAkij(vol, i, j, k) += value;
            }
        }
    }
}

// Calculate the filter in 2D and apply ======================================
void ProgRecWbp::filterOneImage(Projection &proj, Tabsinc &TSINC)
{
    MultidimArray<std::complex<double> > IMG;
    Matrix2D<double> A(3, 3);
    double factor, argum, weight, x, y;

    factor = diameter;

    //Euler_angles2matrix(proj.rot(), -proj.tilt(), proj.psi(), A);
    //A = A.inv();
    Euler_angles2matrix(-proj.rot(), proj.tilt(), -proj.psi(), A);
    FourierTransform(proj(), IMG);
    CenterFFT(IMG, true);

    // loop over all transformation matrices
    double a00 = MAT_ELEM(A,0,0);
    double a01 = MAT_ELEM(A,0,1);
    double a10 = MAT_ELEM(A,1,0);
    double a11 = MAT_ELEM(A,1,1);
    double a20 = MAT_ELEM(A,2,0);
    double a21 = MAT_ELEM(A,2,1);
    for (int k = 0; k < no_mats; k++)
    {
        mat_f[k].x = a00 * mat_g[k].x + a10 * mat_g[k].y + a20 * mat_g[k].z;
        mat_f[k].y = a01 * mat_g[k].x + a11 * mat_g[k].y + a21 * mat_g[k].z;
    }

    double K = ((double) diameter) / dim;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(IMG)
    {
        y = K * i;
        x = K * j;
        weight = 0.;
        for (int k = 0; k < no_mats; k++)
        {
            argum = x * mat_f[k].x + y * mat_f[k].y;
            double daux;
            TSINCVALUE(TSINC, argum, daux);
            weight += mat_g[k].count * daux;
        }

        if (fabs(weight) < threshold)
        {
            count_thr++;
            A2D_ELEM(IMG, i, j) /= SGN(weight) * (threshold * factor);
        }
        else
        {
            A2D_ELEM(IMG, i, j) /= (weight * factor);
        }
    }

    // Calculate back-projection with the filtered projection
    CenterFFT(IMG, false);
    InverseFourierTransform(IMG, proj());
}

bool ProgRecWbp::getImageToProcess(size_t &objId, size_t &objIndex)
{
    if (time_bar_done == 0)
        iter = new MDIterator(SF);
    else
        iter->moveNext();

    ++time_bar_done;
    objIndex = iter->objIndex;
    return ((objId = iter->objId) != BAD_OBJID);
}

void ProgRecWbp::showProgress()
{
    if (verbose > 0 && time_bar_done % time_bar_step == 0)
        progress_bar(time_bar_done);
}

// Calculate the filter in 2D and apply ======================================
void ProgRecWbp::apply2DFilterArbitraryGeometry()
{
    double rot, tilt, psi, xoff, yoff, flip, weight;
    Projection proj;
    Matrix2D<double> L(4, 4), R(4, 4), A;
    FileName fn_img;

    MultidimArray<double> &mReconstructedVolume = reconstructedVolume();
    mReconstructedVolume.initZeros(dim, dim, dim);
    mReconstructedVolume.setXmippOrigin();
    count_thr = 0;

    // Initialize time bar
    if (verbose > 0)
    {
        std::cerr << "--> Back-projecting ..." << std::endl;
        init_progress_bar(time_bar_size);
    }

    mat_f = (WBPInfo*) malloc(no_mats * sizeof(WBPInfo));
    Tabsinc TSINC(0.0001, dim);

    size_t objId, objIndex;
    while (getImageToProcess(objId, objIndex))
    {
        SF.getValue(MDL_IMAGE, fn_img, objId);
        proj.read(fn_img, false);
        getAnglesForImage(objId, rot, tilt, psi, xoff, yoff, flip, weight);
        proj.setRot(rot);
        proj.setTilt(tilt);
        proj.setPsi(psi);
        proj.setShifts(xoff, yoff);
        proj.setFlip(flip);
        proj.setWeight(weight);
        proj.getTransformationMatrix(A, true);
        if (!A.isIdentity())
            selfApplyGeometry(BSPLINE3, proj(), A, IS_INV, WRAP);
        if (do_weights)
            proj() *= proj.weight();
        proj().setXmippOrigin();
        filterOneImage(proj, TSINC);
        simpleBackprojection(proj, mReconstructedVolume, diameter);

        showProgress();
    }
    if (verbose > 0)
        progress_bar(time_bar_size);

    // Symmetrize if necessary
    if (fn_sym != "")
    {
        MultidimArray<double> Vaux;
        Vaux.resize(mReconstructedVolume);
        symmetrizeVolume(SL, mReconstructedVolume, Vaux);
        mReconstructedVolume = Vaux;
        Mask mask_prm;
        mask_prm.mode = INNER_MASK;
        mask_prm.R1 = diameter / 2.;
        mask_prm.type = BINARY_CIRCULAR_MASK;
        mask_prm.generate_mask(mReconstructedVolume);
        mask_prm.apply_mask(mReconstructedVolume, mReconstructedVolume, 0.);
    }
}

