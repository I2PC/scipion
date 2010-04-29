/***************************************************************************
 *
 * Authors:    Roberto Marabini
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

#include "crystal_angular_projection_matching.h"

// Read arguments ==========================================================
void Prog_projection_matching_crystal_prm::read(int argc, char **argv)
{

    // Read command line
    SFref.read(getParameter(argc, argv, "-ref"));
    SFref.ImgSize(dim, dim);
    SFexp.read(getParameter(argc, argv, "-exp"));
    fn_root = getParameter(argc, argv, "-o", "out");
    psi_distance  = textToFloat(getParameter(argc, argv,   "-psi_distance", "0."));
    rot_distance  = textToFloat(getParameter(argc, argv,   "-rot_distance", "0."));
    tilt_distance = textToFloat(getParameter(argc, argv,   "-tilt_distance", "0."));
    shift_distance  = textToFloat(getParameter(argc, argv, "-shift_distance", "0."));
    scale_distance  = textToFloat(getParameter(argc, argv, "-scale_distance", "0."));
    psi_sampling  = textToFloat(getParameter(argc, argv,   "-psi_sam", "1."));
    if (shift_sampling == 0.)
        REPORT_ERROR(1, " shift_sampling must be inizialized!");
    if (psi_sampling == 0.)
        REPORT_ERROR(1, " psi_sampling must be inizialized!");
    max_shift  = textToFloat(getParameter(argc, argv, "-max_shift", "0."));
    shift_sampling  = textToFloat(getParameter(argc, argv, "-shift_sam", "1."));
    scale_sampling  = textToFloat(getParameter(argc, argv, "-scale_sam", "1."));
    modify_header   = !checkParameter(argc, argv, "-dont_modify_header");

}

// Usage ===================================================================
void Prog_projection_matching_crystal_prm::usage()
{
    std::cerr << "Usage:  projection_matching [options] " << std::endl;
    std::cerr << "   -ref <selfile>                : Selfile with reference images \n"
              << "   -exp <selfile>                : Selfile with experimenal images \n"
              << "   -o filerootname               : Output file root name \n"
              << "   -psi_sam number               : number of samples \n"
              << "   -psi_distance degrees  : psi_distance to experimental \
    proj\n"
              << "   -rot_distance degrees  : rot_distance to experimental \
    proj\n"
              << "   -tilt_distance degrees : tilt_distance to experimental \
    proj\n"
              << "   -max_shift double      : maximun allowed shift in pixels \n"
              << " [ -dont_modify_header ]       : Do not store alignment parameters in the image headers \n";
    ;
}

// Show ====================================================================
void Prog_projection_matching_crystal_prm::show()
{

    std::cerr << "  Reference images            : " << SFref.name() << " (" << SFref.ImgNo() << ")" << std::endl;
    std::cerr << "  Experimental images         : " << SFexp.name() << " (" << SFexp.ImgNo() << ")" << std::endl;
    std::cerr << "  Output rootname             : " << fn_root << std::endl;
    std::cerr << "  Psi distance                : " << psi_distance << std::endl;
    std::cerr << "  Shift distance              : " << shift_distance << std::endl;
    std::cerr << "  Scale distance              : " << scale_distance << std::endl;
    std::cerr << "  Rot distance                : " << rot_distance << std::endl;
    std::cerr << "  Tilt distance               : " << tilt_distance << std::endl;
    std::cerr << "  Psi number samples          : " << psi_sampling << std::endl;
    std::cerr << "  Shift number samples        : " << shift_sampling << std::endl;
    std::cerr << "  Scale number samples        : " << scale_sampling << std::endl;
    std::cerr << "  Do not modify the image headers (only output docfile): " << !modify_header << std::endl;

    std::cerr << " =================================================================" << std::endl;
}

// Side info stuff ===================================================================
void Prog_projection_matching_crystal_prm::produce_Side_info()
{
    Projection       proj;

    SFref.go_beginning();
    double           mean_ref, stddev_ref, dummy, psi = 0.;
    // Read projections from selfile
    int nl = SFref.ImgNo();
    ref_img.clear();
    ref_rot = (double*)malloc(nl * sizeof(double));
    ref_tilt = (double*)malloc(nl * sizeof(double));
    ref_mean = (double*)malloc(nl * sizeof(double));
    ref_stddev = (double*)malloc(nl * sizeof(double));
    SFref.go_beginning();
    nr_dir = 0;

    while (!SFref.eof())
    {
        //readimage
        FileName fn_img=SFref.NextImg();
        if (fn_img=="") break;
        proj.read(fn_img);
        proj().setXmippOrigin();
        ref_rot[nr_dir] = proj.rot();
        ref_tilt[nr_dir]  = proj.tilt();
        proj().computeStats(mean_ref, stddev_ref, dummy, dummy);
        proj() -= mean_ref;
        ref_img.push_back(proj());
        ref_stddev[nr_dir] = stddev_ref;
        ref_mean[nr_dir] = mean_ref;
        nr_dir++;
    }

    //compute shift range
    double min_shift, max_shift, increment_shift, my_shift2;
    double shift_distance2 = shift_distance * shift_distance;
    shift_vector.clear();
    Matrix1D<double> v_aux(2);
    //shift already applied to the image
    //psi has not been applied
    min_shift = (-1.0) * shift_distance;
    max_shift =   1.0 * shift_distance;
    if (shift_sampling == 1) increment_shift = 0;
    else increment_shift = 2 * shift_distance / (shift_sampling - 1);

    for (int ishift_x = 0; ishift_x < shift_sampling; ishift_x++)
        for (int ishift_y = 0; ishift_y < shift_sampling; ishift_y++)
        {
            my_shift2 = (min_shift + increment_shift * (double)ishift_x) *
                        (min_shift + increment_shift * (double)ishift_x) +
                        (min_shift + increment_shift * (double)ishift_y) *
                        (min_shift + increment_shift * (double)ishift_y);
            if (my_shift2 > shift_distance2)
                continue;
            XX(v_aux) = min_shift + increment_shift * (double)ishift_x;
            YY(v_aux) = min_shift + increment_shift * (double)ishift_y;
            shift_vector.push_back(v_aux);
        }


//#define DEBUG
#ifdef DEBUG
    for (int ii = 0 ; ii < nr_dir; ii++)
        std::cout << ref_rot[ii] << " " << ref_tilt[ii] << std::endl;
    for (int jj = 0 ; jj < shift_vector.size(); jj++)
        std::cout << shift_vector[jj] << std::endl;
#endif
#undef DEBUG
}




void Prog_projection_matching_crystal_prm::PM_process_one_image(Matrix2D<double> &Mexp,
        float img_rot,
        float img_tilt,
        float img_psi,
        float img_scale,
        int &opt_dirno,
        double &opt_psi,
        double &opt_scale,
        double &opt_xoff,
        double &opt_yoff,
        double &maxCC)
{


    // Rotational search ====================================================
    Matrix2D<double> Mimg, Mref, Maux, Mcorr;
    double act_rot_range, psi, thisCC, oldCC, aveCC = 0., varCC = 0.;
    double stddev_img, mean_img, dummy, xmax, ymax;
    int c = 0, ioptpsi = 0, ioptflip = 0;
    bool search;
    std::vector<Matrix2D<double> >::iterator ipp;

    maxCC = -99.e99;
    Mimg.resize(dim, dim);
    Mimg.setXmippOrigin();
    Mref.resize(dim, dim);
    Mref.setXmippOrigin();
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    // Calculate mean_img,stddev_img and apply rotmask
    Maux = Mexp;
    Mexp.computeStats(mean_img, stddev_img, dummy, dummy);
    Maux -= mean_img;

    //compute psi range
    double min_psi, max_psi, increment_psi, my_psi;
    min_psi = img_psi - psi_distance;
    max_psi = img_psi + psi_distance;
    if (psi_sampling == 1) increment_psi = 0;
    else increment_psi = 2 * psi_distance / (psi_sampling - 1);
    double min_scale, max_scale, increment_scale, my_scale;
    min_scale = img_scale - scale_distance;
    max_scale = img_scale + scale_distance;
    if (scale_sampling == 1) increment_scale = 0;
    else increment_scale = 2 * scale_distance / (scale_sampling - 1);
    double dim2 = dim * dim;
    Matrix2D<double> my_geom_mat, shift_mat;
    //PSI LOOP
    for (int ipsi = 0; ipsi < psi_sampling; ipsi++)
    {
        my_psi = min_psi + increment_psi * (double)ipsi;
        Euler_angles2matrix(0., 0., my_psi, my_geom_mat);
        //my_geom_mat=rotation2DMatrix(my_psi);
        //Mimg=Maux.rotate(my_psi,WRAP);
        //SHIFT LOOP
        for (int ishift = 0; ishift < shift_vector.size(); ishift++)
        {
            my_geom_mat(0, 2) = -XX(shift_vector[ishift]);
            my_geom_mat(1, 2) = -YY(shift_vector[ishift]);
            //SCALE LOOP
            for (int iscale = 0; iscale < scale_sampling; iscale++)
            {
                my_scale = min_scale + increment_scale * (double)iscale;
                my_geom_mat *= my_scale;
                applyGeometryBSpline(Mimg, my_geom_mat, Maux, 3, IS_INV, true, (double)0.);
                ipp = ref_img.begin();
                //ROT-TILT LOOP (ref projections)
                for (int dir_counter = 0; dir_counter < nr_dir; dir_counter++)
                {
                    Mref = *(ipp);
                    ipp++;
                    // For some strange reason I need to access the vector via its pointer
                    // otherwise it goes 50x slower on jumilla (Alpha-Unix)

                    if ((ABS(realWRAP(ref_tilt[dir_counter] - img_tilt, -180., 180.)) <=
                         tilt_distance) &&
                        (ABS(realWRAP(ref_rot [dir_counter] - img_rot, -180., 180.)) <=
                         rot_distance))
                    {
                        thisCC = 0.;
                        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mimg)
                        {
                            thisCC += dMij(Mref, i, j) * dMij(Mimg, i, j);
                        }//FOR_ALL
                        thisCC /= ref_stddev[dir_counter] * stddev_img * dim2;
                        //#define DEBUG
#ifdef DEBUG
                        xmax = XX(shift_vector[ishift]) * my_scale;
                        ymax = YY(shift_vector[ishift]) * my_scale;
                        std::cout <<  SFref.get_file_number(dir_counter) <<
                                  " rot= "  << ref_rot[dir_counter]  <<
                                  " tilt= " << ref_tilt[dir_counter] <<
                                  " psi= "  << my_psi << " CC= " << thisCC << std::endl;
                        " xmax= "
                                << xmax << " ymax= " << ymax << std::endl;
                        " scale= "
                                << my_scale << std::endl;
#endif
#undef DEBUG
                        if (thisCC > maxCC)
                        {
                            xmax = XX(shift_vector[ishift]) * my_scale;
                            ymax = YY(shift_vector[ishift]) * my_scale;
                            maxCC = thisCC;
                            opt_psi = my_psi;
                            opt_scale = my_scale;
                            opt_dirno = dir_counter;
                            opt_xoff =   xmax;
                            opt_yoff =   ymax;
                            // opt_xoff =   xmax*COSD(opt_psi)+ymax*SIND(opt_psi);
                            // opt_yoff =  -xmax*SIND(opt_psi)+ymax*COSD(opt_psi);
                        }//if (thisCC>maxCC)

                    }//if ( (ABS(real
                }//for dir_counter
            }//for iscale
        }// for ishift
    }//for ipsi

    // Interpolated translational search for optimal angles ===================================



    //#define DEBUG
#ifdef DEBUG
    std::cout <<  SFref.get_file_number(opt_dirno) <<
              " rot= "  << ref_rot[opt_dirno]  <<
              " tilt= " << ref_tilt[opt_dirno] <<
              " psi= "  << opt_psi <<
              " opt_xoff= "  << opt_xoff <<
              " opt_yoff= "  << opt_yoff <<
              " opt_scale= "  << opt_scale <<
              " CC= " << maxCC << std::endl;
#endif
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    FileName fn_tmp;
    ImageXmipp save;
    save() = Mimg;
    fn_tmp.compose("tmp", ipsi, "xmp");
    save.write(fn_tmp);
#endif
#undef DEBUG
}
void Prog_projection_matching_crystal_prm::PM_loop_over_all_images(DocFile &DFo,
        double &sumCC)
{
    // Loop over all images
    int imgno = 0, nn;
    double opt_psi, opt_scale, opt_xoff, opt_yoff, maxCC, Zscore;
    int opt_dirno;
    ImageXmipp img;
    FileName fn_img, fn_tmp;
    sumCC = 0.;
    Matrix1D<double> dataline(8);

    // Initialize
    nn = SFexp.ImgNo();
    init_progress_bar(nn);

    SFexp.go_beginning();
    while ((!SFexp.eof()))
    {
        imgno++;
        fn_img = SFexp.NextImg();
        if (fn_img=="") break;
        //last true means shift are applied but not psi
        img.read(fn_img, false, false, false, true);
        img().setXmippOrigin();
        // Perform the projection matching for each image separately
        // NOTE (float)1. should be the scale but since
        // scale field in the image is not reliable a put 1
        PM_process_one_image(img(), img.Phi(), img.Theta(), img.psi(), (float)1.,
                             opt_dirno, opt_psi, opt_scale, opt_xoff, opt_yoff, maxCC);


        opt_xoff += img.Xoff();
        opt_yoff += img.Yoff();

        sumCC += maxCC;
        dataline(0) = ref_rot[opt_dirno];    // rot
        dataline(1) = ref_tilt[opt_dirno];   // tilt
        dataline(2) = opt_psi;               // psi
        dataline(3) = opt_scale;             // scale
        dataline(4) = opt_xoff;              // Xoff
        dataline(5) = opt_yoff;              // Yoff
        dataline(6) = opt_dirno + 1;         // optimal direction number
        dataline(7) = maxCC;                 // maximum CC
        DFo.append_comment(img.name());
        DFo.append_data_line(dataline);

        if (modify_header)
        {
            // Re-read image to get the untransformed image matrix again
            img.read(fn_img);
            img.set_eulerAngles(ref_rot[opt_dirno], ref_tilt[opt_dirno], opt_psi);
            img.set_originOffsets(opt_xoff, opt_yoff);
            img.write(fn_img);
        }

        progress_bar(imgno);
    }// while ((!SFexp.eof()))

    progress_bar(nn);

    // free memory
    free(ref_mean);
    free(ref_rot);
    free(ref_tilt);
    ref_img.clear();
}
