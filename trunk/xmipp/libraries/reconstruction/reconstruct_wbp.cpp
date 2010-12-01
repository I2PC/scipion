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

// Read arguments ==========================================================
void ProgRecWbp::readParams()
{
    fn_sel = getParam("-i");
    fn_doc = getParam("--doc");
    fn_out =  getParam("-o");
    fn_sym =  getParam("--sym");
    threshold = getDoubleParam("--threshold");
    diameter = 2 * getIntParam("--radius");
    sampling = getDoubleParam("--filsam");
    do_all_matrices = checkParam("--use_each_image");
    do_weights = checkParam("--weight");
    // For improved killing control
    //fn_control = getParam("--control");
}

// Show ====================================================================
void ProgRecWbp::show()
{
    if (verbose > 0)
    {
        // To screen
        std::cerr << " =================================================================" << std::endl;
        std::cerr << " Weighted-back projection (arbitrary geometry) " << std::endl;
        std::cerr << " =================================================================" << std::endl;
        std::cerr << " Input selfile             : " << fn_sel << std::endl;
        if (fn_doc != "")
            std::cerr << " Input docfile             : " << fn_doc << std::endl;
        std::cerr << " Output volume             : " << fn_out << std::endl;
        if (diameter > 0)
            std::cerr << " Reconstruction radius     : " << diameter / 2 << std::endl;
        std::cerr << " Relative filter threshold : " << threshold << std::endl;
        if (fn_sym != "")
            std::cerr << " Symmetry file:            : " << fn_sym << std::endl;
        if (do_all_matrices)
            std::cerr << " --> Use all projection directions in arbitrary geometry filter" << std::endl;
        else
            std::cerr << " --> Use sampled directions for filter, sampling = " << sampling << std::endl;
        if (do_weights)
            std::cerr << " --> Use weights stored in the image headers" << std::endl;
        std::cerr << " -----------------------------------------------------------------" << std::endl;
    }
}

// Usage ====================================================================
void ProgRecWbp::defineParams()
{
    // To screen
    addUsageLine("Generate 3D reconstruction from projections using the Weighted BackProjection algorithm.");
    addUsageLine("Example of use: Sample using projections included in g1t.sel input selection file and default options");
    addUsageLine("   reconstruct_wbp -i g1t.sel -o recons_wbp.vol");

    addParamsLine("   -i <input_selfile>          : selection file with input images ");
    addParamsLine(" [ -o <name=\"wbp.vol\"> ]     : filename for output volume ");
    addParamsLine(" [ --doc <docfile=\"\"> ]      : Ignore headers and get angles from this docfile ");
    addParamsLine(" [ --radius <int=-1> ]         : Reconstruction radius. int=-1 means radius=dim/2 ");
    addParamsLine(" [ --sym <symfile=\"\"> ]      : Enforce symmetry ");
    addParamsLine(" [ --threshold <float=0.005> ] : Lower (relative) threshold for filter values ");
    addParamsLine(" [ --filsam <float=5> ]        : Angular sampling rate for geometry filter ");
    addParamsLine(" [ --use_each_image]           : Use each image instead of sampled representatives for filter ");
    addParamsLine(" [ --weight]                   : Use weights stored in image headers ");
    //addParamsLine(" [ --control <fn=\"\">]             : For improved killing control ");
}

void ProgRecWbp::run()
{
    Image<double> vol;

    produceSideInfo();
    apply_2Dfilter_arbitrary_geometry(SF, vol());

    if (verbose > 0)
        std::cerr << "Fourier pixels for which the threshold was not reached: "
        << (float)(count_thr*100.) / (SF.size()*dim*dim) << " %" << std::endl;

    vol.write(fn_out);
}

void ProgRecWbp::produceSideInfo()
{
    // Read-in stuff
    //remove images with weight=0
    double dum, weight;

    //remove images with weight=0
    if (do_weights)
    {
        MetaData SF_aux;
        SF_aux.read(fn_sel);
        FOR_ALL_OBJECTS_IN_METADATA(SF_aux);
        {
            SF_aux.getValue(MDL_WEIGHT,weight);
            if (weight != 0)
            {
                SF.addObject();
                FileName fnImg;
                SF_aux.getValue(MDL_IMAGE,fnImg);
                SF.setValue(MDL_IMAGE,fnImg);
            }
        }
        if (SF.size() == 0)
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,"there is no input file with weight!=0");
    }
    else
        SF.read(fn_sel);

    ImgSize(SF,dim, dim);
    if (fn_sym != "")
        SL.read_sym_file(fn_sym);
    if (diameter == 0)
        diameter = dim;

    // Fill arrays of transformation matrices
    if (do_all_matrices)
        get_all_matrices(SF);
    else
        get_sampled_matrices(SF);

}


void ProgRecWbp::get_angles_for_image(const FileName &fn, double &rot, double &tilt, double &psi,
        double &xoff, double &yoff, double &flip, double &weight)
{
    if (fn_doc == "")
    {
        Image<double> I;
        I.read(fn,false);
        rot    = I.rot();
        tilt   = I.tilt();
        psi    = I.psi();
        xoff   = I.Xoff();
        yoff   = I.Yoff();
        flip   = I.flip();
        weight = I.weight();
    }
    else
    {
        unsigned long int id=SF.gotoFirstObject(MDValueEQ(MDL_IMAGE,(std::string)fn));
        if (id!=NO_OBJECT_FOUND && id!=NO_OBJECTS_STORED)
        {
            SF.getValue(MDL_ANGLEROT,rot);
            SF.getValue(MDL_ANGLETILT,tilt);
            SF.getValue(MDL_ANGLEPSI,psi);
            SF.getValue(MDL_SHIFTX,xoff);
            SF.getValue(MDL_SHIFTY,yoff);
            flip=0;
            SF.getValue(MDL_FLIP,flip);
            weight=0;
            SF.getValue(MDL_WEIGHT,weight);
        }
        else
            REPORT_ERROR(ERR_MD_NOOBJ, (std::string)"Cannot find " + fn + " in docfile " + fn_doc);
    }
}

void ProgRecWbp::get_sampled_matrices(MetaData &SF)
{
    MetaData          DFlib;
    FileName          fn_tmp;
    Matrix2D<double>  A(3, 3);
    Matrix2D<double>  L(4, 4), R(4, 4);
    double            newrot, newtilt, newpsi, rot, tilt, psi, dum, weight, totimgs = 0.;
    int               NN, dir, optdir;
    std::vector<double>    count_imgs;

    if (verbose > 0)
        std::cerr << "--> Sampling the filter ..." << std::endl;

    // Create an (asymmetric part of an) even projection direction distribution
    make_even_distribution(DFlib, sampling, SL, true);
    NN = DFlib.size();
    count_imgs.resize(NN);
    // Each experimental image contributes to the nearest of these directions
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        FileName fn_img;
        SF.getValue(MDL_IMAGE,fn_img);
        get_angles_for_image(fn_img, rot, tilt, dum, dum, dum, dum, weight);
        int idx=find_nearest_direction(rot,tilt,DFlib,SL);
        if (do_weights)
            count_imgs[idx] += weight;
        else
            count_imgs[idx] += 1.;
    }

    // Now calculate transformation matrices for all representative directions
    no_mats = 0;
    for (int i = 0; i < NN; i++)
        if (count_imgs[i] > 0.)
            no_mats += SL.SymsNo() + 1;
    mat_g = (column*)malloc(no_mats * sizeof(column));

    no_mats = 0;
    for (int i = 0; i < NN; i++)
    {
        if (count_imgs[i] > 0.)
        {
            DFlib.getValue(MDL_ANGLEROT,newrot,i);
            DFlib.getValue(MDL_ANGLETILT,newtilt,i);
            DFlib.getValue(MDL_ANGLEPSI,newpsi,i);
            Euler_angles2matrix(newrot, -newtilt, newpsi, A);
            mat_g[no_mats].zero = A(2, 0);
            mat_g[no_mats].one = A(2, 1);
            mat_g[no_mats].two = A(2, 2);
            mat_g[no_mats].count = count_imgs[i];
            totimgs += mat_g[no_mats].count;
            no_mats++;

            // Expand symmetric directions
            for (int j = 0; j < SL.SymsNo(); j++)
            {
                SL.get_matrices(j, L, R);
                L.resize(3, 3);
                R.resize(3, 3);
                Euler_apply_transf(L, R, newrot, newtilt, 0., rot, tilt, psi);
                Euler_angles2matrix(rot, -tilt, psi, A);
                mat_g[no_mats].zero = A(2, 0);
                mat_g[no_mats].one = A(2, 1);
                mat_g[no_mats].two = A(2, 2);
                mat_g[no_mats].count = count_imgs[i];
                totimgs += mat_g[no_mats].count;
                no_mats++;
            }
        }
    }

    // Adjust relative threshold
    threshold *= totimgs;
}

// Fill array with transformation matrices needed for arbitrary geometry filter
void ProgRecWbp::get_all_matrices(MetaData &SF)
{
    Matrix2D<double> A(3, 3);
    Matrix2D<double> L(4, 4), R(4, 4);
    double           rot, tilt, psi, weight, dum, newrot, newtilt, newpsi, totimgs = 0.;
    int              NN;

    no_mats = 0;

    NN = SF.size();
    NN *= (SL.SymsNo() + 1);
    mat_g = (column*)malloc(NN * sizeof(column));

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        FileName fn_img;
        SF.getValue(MDL_IMAGE,fn_img);
        get_angles_for_image(fn_img, rot, tilt, psi, dum, dum, dum, weight);
        Euler_angles2matrix(rot, -tilt, psi, A);
        mat_g[no_mats].zero = A(2, 0);
        mat_g[no_mats].one = A(2, 1);
        mat_g[no_mats].two = A(2, 2);
        if (do_weights)
            mat_g[no_mats].count = weight;
        else
            mat_g[no_mats].count = 1.;
        totimgs += mat_g[no_mats].count;
        no_mats++;
        // Also add symmetry-related projection directions
        for (int i = 0; i < SL.SymsNo(); i++)
        {
            SL.get_matrices(i, L, R);
            L.resize(3, 3);
            R.resize(3, 3);
            Euler_apply_transf(L, R, rot, -tilt, psi, newrot, newtilt, newpsi);
            Euler_angles2matrix(newrot, newtilt, newpsi, A);
            mat_g[no_mats].zero = A(2, 0);
            mat_g[no_mats].one = A(2, 1);
            mat_g[no_mats].two = A(2, 2);
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
void ProgRecWbp::simple_backprojection(Projection &img, MultidimArray<double> &vol,
        int diameter)
{
    int i, j, k, l, m;
    Matrix2D<double> A(3, 3);
    float dim2, x, y, z, xp, yp;
    float value1, value2, scalex, scaley, scale1, value;
    float radius2, x2, y2, z2, z2_plus_y2;

    // Use minus-tilt, because code copied from OldXmipp
    Euler_angles2matrix(img.rot(), -img.tilt(), img.psi(), A);
    A = A.inv();

    radius2 = diameter / 2.;
    radius2 = radius2 * radius2;
    dim2 = dim / 2;

    for (i = 0; i < dim; i++)
    {
        z = -i + dim2;   /*** Z points upwards ***/
        z2 = z * z;
        for (j = 0; j < dim; j++)
        {
            y = j - dim2;
            y2 = y * y;
            z2_plus_y2 = z2 + y2;
            x = 0 - dim2;   /***** X for k == 0 *****/
            xp = x * A(0, 0) + y * A(1, 0) + z * A(2, 0) + dim2;
            yp = x * A(0, 1) + y * A(1, 1) + z * A(2, 1) + dim2;
            for (k = 0; k < dim; k++, xp += A(0, 0), yp += A(0, 1), x++)
            {
                x2 = x * x;
                if (x2 + z2_plus_y2 > radius2)
                    continue;
                if ((xp >= (dim - 1) || xp < 0) || (yp >= (dim - 1) || yp < 0))
                    continue;

                /**** interpolation ****/
                l = (int)yp;
                m = (int)xp;
                scalex = xp - m;
                scaley = yp - l;
                scale1 = 1. - scalex;
                value1 = scalex * dAij(img(), l, m + 1) + scale1 * dAij(img(), l, m);
                value2 = scalex * dAij(img(), l + 1, m + 1) + scale1 * dAij(img(), l + 1, m);
                value  = scaley * value2 + (1. - scaley) * value1;
                dAkij(vol, i, j, k) += value;
            }
        }
    }
}

// Calculate the filter in 2D and apply ======================================
void ProgRecWbp::filter_one_image(Projection &proj)
{
    Image<double> save;
    save()=proj();
    save.write("PPPwbp2_before.xmp");

    MultidimArray< std::complex<double> > IMG;
    Matrix2D<double>  A(3, 3);
    float             factor, argum, weight, x, y;

    factor = (float)diameter;

    // Tabulated sinc
    Tabsinc TSINC(0.001, dim);

    Euler_angles2matrix(proj.rot(), -proj.tilt(), proj.psi(), A);
    A = A.inv();
    FourierTransform(proj(), IMG);
    CenterFFT(IMG, true);

    // loop over all transformation matrices
    for (int k = 0; k < no_mats; k++)
    {
        mat_f[k].zero = A(0, 0) * mat_g[k].zero +
                        A(1, 0) * mat_g[k].one +
                        A(2, 0) * mat_g[k].two;
        mat_f[k].one = A(0, 1) * mat_g[k].zero +
                       A(1, 1) * mat_g[k].one +
                       A(2, 1) * mat_g[k].two;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY2D(IMG)
    {
        y = (float)i;
        x = (float)j;
        weight = 0.;
        for (int k = 0; k < no_mats; k++)
        {
            argum = diameter / (float)dim *
                    (x * mat_f[k].zero + y * mat_f[k].one);
            // The following line is the most expensive of all...
            weight += mat_g[k].count * TSINC(argum);
        }

        if (weight < threshold)
        {
            count_thr++;
            A2D_ELEM(IMG, i, j) /= (threshold * factor);
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

// Calculate the filter in 2D and apply ======================================
void ProgRecWbp::apply_2Dfilter_arbitrary_geometry(MetaData &SF, MultidimArray<double> &vol)
{

    int               c, nn, imgno;
    double            rot, tilt, psi, newrot, newtilt, newpsi, xoff, yoff, flip, weight;
    Projection        proj;
    Matrix2D<double>  L(4, 4), R(4, 4), A;
    Mask       mask_prm;
    FileName          fn_img;

    vol.resize(dim, dim, dim);
    vol.setXmippOrigin();
    vol.initZeros();
    count_thr = 0;

    // Initialize time bar
    if (verbose > 0)
        std::cerr << "--> Back-projecting ..." << std::endl;
    nn = SF.size();
    if (verbose > 0)
        init_progress_bar(nn);
    c = XMIPP_MAX(1, nn / 60);

    mat_f = (column*)malloc(no_mats * sizeof(column));

    imgno = 0;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        // Check whether to kill job
        //exit_if_not_exists(fn_control);
        SF.getValue(MDL_IMAGE,fn_img);
        if (fn_doc == "")
        {
            proj.read(fn_img, true);
        }
        else
        {
            proj.read(fn_img, false);
            get_angles_for_image(fn_img, rot, tilt, psi, xoff, yoff, flip, weight);
            proj.setRot(rot);
            proj.setTilt(tilt);
            proj.setPsi(psi);
            proj.setShifts(xoff,yoff);
            proj.setFlip(flip);
            proj.setWeight(weight);
            A = proj.getTransformationMatrix(true);
            if (!A.isIdentity())
                selfApplyGeometry(BSPLINE3, proj(), A, IS_INV, WRAP);
        }
        if (do_weights)
            proj() *= proj.weight();
        proj().setXmippOrigin();
        filter_one_image(proj);
        simple_backprojection(proj, vol, diameter);

        if (verbose > 0)
            if (imgno % c == 0)
                progress_bar(imgno);
        imgno++;

    }
    if (verbose > 0)
        progress_bar(nn);

    // Symmetrize if necessary
    if (fn_sym != "")
    {
        MultidimArray<double> Vaux;
        Vaux.resize(vol);
        symmetrize(SL, vol, Vaux);
        vol = Vaux;
        vol *= (SL.SymsNo() + 1);
        mask_prm.mode = INNER_MASK;
        mask_prm.R1 = diameter / 2.;
        mask_prm.type = BINARY_CIRCULAR_MASK;
        mask_prm.generate_mask(vol);
        mask_prm.apply_mask(vol, vol, 0.);
    }

    // free memory
    free(mat_g);
    free(mat_f);
}

