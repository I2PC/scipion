/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.csic.es (2008)
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

#include "angular_class_average.h"

// Read arguments ==========================================================
void Prog_angular_class_average_prm::read(int argc, char **argv)  {

    // Read command line
    DF.read(getParameter(argc, argv, "-i"));
    DFlib.read(getParameter(argc, argv, "-lib"));
    if (checkParameter(argc, argv, "-add_to"))
    {
        do_add = true;
        fn_out = getParameter(argc, argv, "-add_to");
    }
    else
    {
	do_add = false;
        fn_out = getParameter(argc, argv, "-o");
    }
    col_ref = textToInteger(getParameter(argc, argv, "-refno","6"));
    col_ref--;
    col_select = textToInteger(getParameter(argc, argv, "-select", "8"));
    col_select--;
    if (checkParameter(argc, argv, "-limitR"))
    {
	limitR = textToFloat(getParameter(argc, argv, "-limitR"));
        if (limitR < -100. || limitR > 100.)
            REPORT_ERROR(1,"limitR should be a percentage: provide values between -100 and 100.");
        if (limitR > 0.)
            do_limitR0 = true;
        else if (limitR < 0.)
        {
            limitR *= -1.;
            do_limitRF = true;
        }
    }
    do_limit0=checkParameter(argc, argv, "-limit0");
    if (do_limit0)
    {
	limit0 = textToFloat(getParameter(argc, argv, "-limit0"));
    }
    do_limitF=checkParameter(argc, argv, "-limitF");
    if (do_limitF)
    {
	limitF = textToFloat(getParameter(argc, argv, "-limitF"));
    }
    
    // Perform splitting of the data?
    do_split = checkParameter(argc, argv, "-split"); 

    // Perform Wiener filtering of average?
    fn_wien = getParameter(argc, argv, "-wien","");
    pad = XMIPP_MAX(1.,textToFloat(getParameter(argc, argv, "-pad","1.")));

    // Skip writing selfiles?
    dont_write_selfiles = checkParameter(argc, argv, "-dont_write_selfiles"); 

    // Internal re-alignment of the class averages
    Ri      = textToInteger(getParameter(argc,argv,"-Ri","-1"));
    Ro      = textToInteger(getParameter(argc,argv,"-Ro","-1"));
    nr_iter = textToInteger(getParameter(argc,argv,"-iter","0"));
    max_shift        = textToFloat(getParameter(argc, argv, "-max_shift","999."));
    max_shift_change = textToFloat(getParameter(argc, argv, "-max_shift_change","999."));
    max_psi_change   = textToFloat(getParameter(argc, argv, "-max_psi_change","360."));
}

// Show ====================================================================
void Prog_angular_class_average_prm::show() {

    std::cerr << "  Input docfile           : "<< DF.name()<<std::endl;
    std::cerr << "  Library docfile         : "<< DFlib.name()<<std::endl;
    if (do_add)
        std::cerr << "  Add class averages to   : "<< fn_out<<std::endl;
    else
        std::cerr << "  Output rootname         : "<< fn_out<<std::endl;
    if (do_split)
	std::cerr << "     -> Split data in random halves and output class averages "<<std::endl;
    if (do_mirrors)
	std::cerr << "     -> Take mirror operation into account "<<std::endl;
    if (dont_write_selfiles)
	std::cerr << "     -> Do not write class selfiles to disc "<<std::endl;
    // election
    if (do_limit0 || do_limitF || do_limitR0 || do_limitRF)
    {
        std::cerr << "  PERFORM IMAGE SELECTION BASED ON INPUT DOCFILE"<<std::endl;
        std::cerr << "    Column number to use    : "<<col_select+1<<std::endl;
        if (do_limitR0)
            std::cerr << "    Discard lowest          : "<<limitR<<" %"<<std::endl;
        else if (do_limitRF)
            std::cerr << "    Discard highest         : "<<limitR<<" %"<<std::endl;
        if (do_limit0)
            std::cerr << "    Discard images below    : "<<limit0<<std::endl;
        if (do_limitF)
            std::cerr << "    Discard images above    : "<<limitF<<std::endl;
    }
    // Realignment
    if (nr_iter > 0)
    {
        std::cerr << "  PERFORM REALIGNMENT OF CLASSES"<<std::endl;
        std::cerr << "    Number of iterations    : "<<nr_iter<<std::endl;
        std::cerr << "    Maximum shift           : "<<max_shift<<std::endl;
        std::cerr << "    Maximum shift change    : "<<max_shift_change<<std::endl;
        std::cerr << "    Maximum psi change      : "<<max_psi_change<<std::endl;
        if (Ri>0)
            std::cerr << "    Inner radius rot-search : "<<Ri<<std::endl;
        if (Ro>0)
            std::cerr << "    Outer radius rot-search : "<<Ro<<std::endl;
    }
    // Wiener filter correction
    if (fn_wien != "")
    {
        std::cerr << "  PERFORM WIENER CORRECTION ON CLASS AVERAGES"<<std::endl;
        std::cerr << "    Padding factor          : "<<pad<<std::endl;
    }

    std::cerr << " ================================================================="<<std::endl;

}

// Usage ===================================================================
void Prog_angular_class_average_prm::usage() {
    printf("Purpose:\n");
    printf(" Makes class average images and corresponding selfiles from angular_projection_matching docfiles.\n");
    printf("Usage:\n");
    printf("   angular_class_average \n");
    printf("        -i <docfile>        : docfile with assigned angles for all experimental particles\n");
    printf("        -lib <docfile>      : docfile with angles used to generate the projection matching library\n");
    printf("        -o <rootname>       : output rootname for class averages and selfiles\n");
    printf("    OR: -add_to <rootname>  : Add output to existing files\n");
    printf("       [-split ]            : Also output averages of random halves of the data\n");
    printf("       [-wien <img=\"\"> ]    : Apply this Wiener filter to the averages\n");    
    printf("       [-pad <float=1.> ]   : Padding factor for Wiener correction\n");
    printf("       [-refno  <int=6>]    : Column number in docfile to use for class number\n");
    printf("       [-dont_write_selfiles]  : Do not write class selfiles to disc\n");
    printf(" IMAGE SELECTION BASED ON INPUT DOCFILE \n");
    printf("       [-select <int=8>]    : Column number to use for image selection (limit0, limitF or limitR)\n");
    printf("       [-limit0 <float>]    : Discard images below <limit0>\n");
    printf("       [-limitF <float>]    : Discard images above <limitF>\n");
    printf("       [-limitR <float>     : if (limitR>0 && limitR< 100): discard lowest  <limitR> % in each class\n");
    printf("                            : if (limitR<0 && limitR>-100): discard highest <limitR> % in each class\n");
    printf(" REALIGNMENT OF CLASSES \n");
    printf("       [-iter <int=0>]                  : Number of iterations for re-alignment\n");
    printf("       [-Ri <int=1>]                    : Inner radius to limit rotational search\n");
    printf("       [-Ro <int=dim/2-1>]              : Outer radius to limit rotational search\n");
    printf("       [-max_shift <float=999.>]        : Maximum shift (larger shifts will be set to 0)\n");
    printf("       [-max_shift_change <float=999.>] : Discard images that change shift more in the last iteration \n");
    printf("       [-max_psi_change <float=360.>]   : Discard images that change psi more in the last iteration \n");
    exit(1);
}

// Side info stuff ===================================================================
void Prog_angular_class_average_prm::produceSideInfo() {

    
    FileName fn_tst, fn_img;

    // Initialize selfiles and docfiles for all class averages
    std::string header="Headerinfo columns: rot (1) , tilt (2), psi (3), Xoff (4), Yoff (5), Weight (6), Flip (7)";
    DFclasses.append_comment(header);
    if (do_split)
    {
        DFclasses1.append_comment(header);
        DFclasses2.append_comment(header);
    }

    // Set up output rootnames
    if (do_split)
    {
	fn_out1 = fn_out+"_split_1";
	fn_out2 = fn_out+"_split_2";
    }

    // get column numbers from NewXmipp-type docfile header
    DF.go_beginning();
    col_rot    = DF.getColNumberFromHeader("rot")  - 1;
    col_tilt   = DF.getColNumberFromHeader("tilt") - 1;
    col_psi    = DF.getColNumberFromHeader("psi")  - 1;
    col_xshift = DF.getColNumberFromHeader("Xoff") - 1;
    col_yshift = DF.getColNumberFromHeader("Yoff") - 1;
    col_mirror = DF.getColNumberFromHeader("Flip") - 1;
    if (col_mirror < 0)
        do_mirrors = false;
    else
        do_mirrors = true;

    if (DF.get_current_line().Is_comment()) fn_tst = (DF.get_current_line()).get_text();
    if (strstr(fn_tst.c_str(), "Headerinfo") == NULL)
    {
	REPORT_ERROR(1,"ERROR: Docfile is of non-NewXmipp type.");
    }

    // Read empty image with the correct dimensions
    DF.next();
    if (DF.get_current_line().Is_comment()) fn_img = ((DF.get_current_line()).get_text()).erase(0, 3);
    else  REPORT_ERROR(1, "Problem with NewXmipp-type document file");
    Iempty.read(fn_img);
    Iempty().setXmippOrigin();
    Iempty().initZeros();
    dim = XSIZE(Iempty());

    // Read Wiener filter image
    if (fn_wien!="") 
    {
        // Get padding dimensions
        paddim = ROUND(pad * dim);
        ImageXmipp It;
        It.read(fn_wien);
        Mwien=It();
        if (XSIZE(Mwien) != paddim)
        {
            std::cerr<<"image size= "<<dim<<" padding factor= "<<pad<<" padded image size= "<<paddim<<" Wiener filter size= "<<XSIZE(Mwien)<<std::endl;
            REPORT_ERROR(1,"Incompatible padding factor for this Wiener filter");
        }
    }
    // Set ring defaults
    if (Ri<1) Ri=1;
    if (Ro<0) Ro=(dim/2)-1;

    // Set limitR
    if (do_limitR0 || do_limitRF)
    {
        std::vector<double> vals;
        DF.go_beginning();
        for (int n = 0; n < DF.dataLineNo(); n++)
        {
            DF.adjust_to_data_line();
            vals.push_back(DF(col_select));
            DF.next();
        }
        int nn = vals.size();
        std::sort(vals.begin(), vals.end());
        if (do_limitR0)
        {
            double val = vals[ROUND((limitR/100.) * vals.size())];
            if (do_limit0)
                limit0 = XMIPP_MAX(limit0, val);
            else
            {
                limit0 = val;
                do_limit0 = true;
            }
        }
        else if (do_limitRF)
        {
            double val = vals[ROUND(((100. - limitR)/100.) * vals.size())];
            if (do_limitF)
                limitF = XMIPP_MIN(limitF, val);
            else
            {
                limitF = val;
                do_limitF = true;
            }
        }
    }

    // Randomization
    if (do_split) randomize_random_generator();

    // Set up FFTW transformers
    Matrix2D<double> Maux;
    Polar<double> P;
    Polar<std::complex<double> > fP;

    Iempty().produceSplineCoefficients(Maux,3);
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro);
    P.calculateFftwPlans(global_plans);
    fourierTransformRings(P,fP,global_plans,false);
    corr.resize(P.getSampleNoOuterRing());
    global_transformer.setReal(corr);
    global_transformer.FourierTransform();

}

void Prog_angular_class_average_prm::getPolar(Matrix2D<double> &img, Polar<std::complex <double> > &fP,
                                              bool conjugated, float xoff, float yoff)
{
    Matrix2D<double> Maux;
    Polar<double> P;

    // Calculate FTs of polar rings and its stddev
    img.produceSplineCoefficients(Maux,3);
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro,xoff,yoff);
    fourierTransformRings(P,fP,global_plans,conjugated);

}

void Prog_angular_class_average_prm::reAlignClass(ImageXmipp &avg1,
                                                  ImageXmipp &avg2,
                                                  SelFile    &SFclass1,
                                                  SelFile    &SFclass2,
                                                  std::vector<ImageXmipp> imgs,
                                                  std::vector<int> splits,
                                                  std::vector<int> numbers,
                                                  int dirno,
                                                  double * my_output)
{


    Polar<std::complex <double> >   fPref, fPrefm, fPimg;
    std::vector<double>             ccfs(splits.size());
    Matrix1D<double>                ang;
    Matrix2D<double>                Mimg, Mref, Maux;
    double                          maxcorr, diff_psi, diff_shift, new_xoff, new_yoff;
    double                          w1, w2, opt_flip=0., opt_psi=0., opt_xoff=0., opt_yoff=0.;
    bool                            do_discard;

    SFclass1.clear();
    SFclass2.clear();
    Mref = avg1() + avg2();
//#define DEBUG
#ifdef DEBUG
        ImageXmipp It;
        It() = Mref;
        It.write("ref.xmp");
#endif


    for (int iter = 0; iter < nr_iter; iter++)
    {
        // Initialize iteration
        getPolar(Mref,fPref,true);
        getPolar(Mref,fPrefm,false);
        avg1().initZeros();
        avg2().initZeros();
        w1 = w2 = 0.;

#ifdef DEBUG
        std::cerr<<" entering iter "<<iter<<std::endl;
#endif      
        for (int imgno = 0; imgno < imgs.size(); imgno++)
        {
            
            do_discard = false;
            maxcorr = -99.e99;
            // Rotationally align
            getPolar(imgs[imgno](),fPimg,false,(float)-imgs[imgno].Xoff(),(float)-imgs[imgno].Yoff());
            // A. Check straight image
            rotationalCorrelation(fPimg,fPref,ang,global_transformer);
	    for (int k = 0; k < XSIZE(corr); k++)
	    {
		if (corr(k)> maxcorr)
		{
		    maxcorr = corr(k);
		    opt_psi = ang(k);
                    opt_flip = 0.;
		}
	    }

            // B. Check mirrored image
            rotationalCorrelation(fPimg,fPrefm,ang,global_transformer);
	    for (int k = 0; k < XSIZE(corr); k++)
	    {
		if (corr(k)> maxcorr)
		{
		    maxcorr = corr(k);
		    opt_psi = realWRAP(360. - ang(k), -180., 180.);
                    opt_flip = 1.;
		}
	    }

            // Check max_psi_change in last iteration
            if (iter == nr_iter - 1)
            {
                diff_psi = ABS(realWRAP(imgs[imgno].psi() - opt_psi, -180., 180.));
                if (diff_psi > max_psi_change) 
                {
                    do_discard = true;
#ifdef DEBUG
                    //std::cerr<<"discard psi: "<<diff_psi<<opt_psi<<" "<<opt_flip<<" "<<imgs[imgno].psi()<<std::endl;
#endif
                }
            }

            // Translationally align
            if (!do_discard)
            {
                if (opt_flip == 1.)
                {
                    // Flip experimental image
                    Matrix2D<double> A(3,3);
                    A.initIdentity();
                    A(0, 0) *= -1.;
                    A(0, 1) *= -1.;
                    applyGeometry(Mimg, A, imgs[imgno](), IS_INV, DONT_WRAP);
                    Mimg.selfRotateBSpline(3,opt_psi,DONT_WRAP);
                }
                else
                {
                    imgs[imgno]().rotateBSpline(3,opt_psi,Mimg,DONT_WRAP);
                }
                if (max_shift > 0) 
                {
                    best_shift(Mref,Mimg,opt_xoff,opt_yoff);
                    if (opt_xoff * opt_xoff + opt_yoff * opt_yoff > max_shift * max_shift) 
                    {
                        new_xoff = new_yoff = 0.;
                    }
                    else
                    {
                        Mimg.selfTranslate(vectorR2(opt_xoff,opt_yoff),true);
                        new_xoff =  opt_xoff*COSD(opt_psi) + opt_yoff*SIND(opt_psi);
                        new_yoff = -opt_xoff*SIND(opt_psi) + opt_yoff*COSD(opt_psi);
                    }

                    // Check max_shift_change in last iteration
                    if (iter == nr_iter - 1)
                    {
                        opt_yoff = imgs[imgno].Yoff();
                        opt_xoff = imgs[imgno].Xoff();
                        if (imgs[imgno].flip() == 1.)                
                            opt_xoff *= -1.;
                        diff_shift = (new_xoff - opt_xoff) * (new_xoff - opt_xoff) +
                                     (new_yoff - opt_yoff) * (new_yoff - opt_yoff);
                        if (diff_shift > max_shift_change * max_shift_change) 
                        {
                            do_discard = true;
#ifdef DEBUG
                            std::cerr <<"discard shift: "<<diff_shift<<" "<<new_xoff<<" "<<opt_xoff<<" "<<imgs[imgno].Xoff()<<" "<<new_yoff<<" "<<opt_yoff<<" "<<imgs[imgno].Yoff()<<std::endl;
#endif
                        }
                    }
                }
            }

            if (!do_discard)
            {
                ccfs[imgno] = correlation_index(Mref,Mimg);
                imgs[imgno].set_psi(opt_psi);
                imgs[imgno].set_flip(opt_flip);
                imgs[imgno].set_Yoff(new_yoff);
                imgs[imgno].set_Xoff(new_xoff);
                if (opt_flip==1.)                
                    imgs[imgno].set_Xoff(-imgs[imgno].Xoff());

                // Check max_shift_change in last iteration
                // Add to averages
                if (splits[imgno] == 0)
                {
                    w1 += 1.;
                    avg1() += Mimg;
                }
                else if (splits[imgno] == 1)
                {
                    w2 += 1.;
                    avg2() += Mimg;
                }
         }
            else
            {
                splits[imgno] = -1;
                ccfs[imgno] = 0.; 
           }
        }
        Mref = avg1() + avg2();
#ifdef DEBUG
        //It() = Mref;
        //It.weight()=w1+w2;
        //It.write("ref.xmp");
#endif
    }

    avg1.set_weight(w1);
    avg2.set_weight(w2);
    
    // Report the new angles, offsets and selfiles
    my_output[4] = imgs.size() * AVG_OUPUT_SIZE;
    for (int imgno = 0; imgno < imgs.size(); imgno++)
    {
        if (splits[imgno] < 0)
            my_output[imgno * AVG_OUPUT_SIZE + 5]  = -numbers[imgno];
        else
            my_output[imgno * AVG_OUPUT_SIZE + 5]  = numbers[imgno];
        my_output[imgno * AVG_OUPUT_SIZE + 6]  = avg1.rot();
        my_output[imgno * AVG_OUPUT_SIZE + 7]  = avg1.tilt();
        my_output[imgno * AVG_OUPUT_SIZE + 8]  = imgs[imgno].psi();
        my_output[imgno * AVG_OUPUT_SIZE + 9]  = imgs[imgno].Xoff();
        my_output[imgno * AVG_OUPUT_SIZE + 10] = imgs[imgno].Yoff();
        my_output[imgno * AVG_OUPUT_SIZE + 11] = dirno;
        my_output[imgno * AVG_OUPUT_SIZE + 12] = imgs[imgno].flip();
        my_output[imgno * AVG_OUPUT_SIZE + 13] = ccfs[imgno];

        if (splits[imgno] == 0) 
            SFclass1.insert(imgs[imgno].name());
        else if (splits[imgno] == 1) 
            SFclass2.insert(imgs[imgno].name());
    }


}
void Prog_angular_class_average_prm::applyWienerFilter(Matrix2D<double> &img)
{
    Matrix2D<std::complex<double> > Faux;
    if (paddim > dim)
    {
        // pad real-space image
        int x0 = FIRST_XMIPP_INDEX(paddim);
        int xF = LAST_XMIPP_INDEX(paddim);
        img.window(x0, x0, xF,xF, 0.);
    }
    FourierTransform(img,Faux);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mwien)
    {
        dMij(Faux,i,j) *= dMij(Mwien,i,j);
    }
    InverseFourierTransform(Faux,img);
    if (paddim > dim)
    {
        // de-pad real-space image
        int x0 = FIRST_XMIPP_INDEX(dim);
        int xF = LAST_XMIPP_INDEX(dim);
        img.window(x0, x0, xF,xF, 0.);
    }
}


void Prog_angular_class_average_prm::processOneClass(int &dirno, 
                                                     double * my_output) {

    ImageXmipp img, avg, avg1, avg2;
    FileName   fn_img, fn_tmp;
    SelFile    SFclass, SFclass1, SFclass2;
    double     rot, tilt, psi, xshift, yshift, mirror, val, w, w1, w2, my_limitR;
    int        ref_number, this_image;
    int        isplit;
    Matrix2D<double> A(3,3);
    std::vector<int> exp_number, exp_split;
    std::vector<ImageXmipp> exp_imgs;
    
    
    // Get reference angles and preset to averages
    DFlib.locate(dirno);
    rot = DFlib(col_rot);
    tilt = DFlib(col_tilt);
    Iempty.set_eulerAngles((float)rot, (float)tilt, (float)0.);
    Iempty.set_originOffsets(0., 0.);
    Iempty.set_flip(0.);
    avg=Iempty;
    avg1=Iempty;
    avg2=Iempty;

    // Loop over all images in the input docfile
    DF.go_beginning();
    w = 0.;
    w1 = 0.;
    w2 = 0.;
    for (int n = 0; n < DF.dataLineNo(); n++)
    {
        DF.next();
        if (DF.get_current_line().Is_comment()) fn_img = ((DF.get_current_line()).get_text()).erase(0, 3);
        else  REPORT_ERROR(1, "Problem with NewXmipp-type document file");
        DF.adjust_to_data_line();
        this_image = DF.get_current_key();
        ref_number = ROUND(DF(col_ref));
        if (ref_number == dirno)
        {
            bool is_select = true;
            if (do_limit0)
                if (DF(col_select) < limit0) is_select = false;
            if (do_limitF)
                if (DF(col_select) > limitF) is_select = false;
            if (is_select)
            {
                psi    = DF(col_psi);
                xshift = DF(col_xshift);
                yshift = DF(col_yshift);
                if (do_mirrors) mirror = DF(col_mirror);
                img.read(fn_img, false, false, false, false);
                img().setXmippOrigin();
                img.set_eulerAngles((float)0., (float)0., (float)psi);
                img.set_originOffsets(xshift, yshift);
                if (do_mirrors) img.set_flip(mirror);
                
                if (do_split) isplit = ROUND(rnd_unif());
                else isplit = 0;

                // For re-alignment of class: store all images in memory
                if (nr_iter > 0)
                {
                    exp_imgs.push_back(img);
                    exp_number.push_back(this_image);
                    exp_split.push_back(isplit);
                }

                // Apply in-plane transformation
                A = img.get_transformation_matrix();
                if (!A.isIdentity())
                    img().selfApplyGeometryBSpline(A, 3, IS_INV,WRAP);
                
                // Add to average
                if (isplit==0)
                {
                    avg1() += img();
                    w1 += 1.;
                    SFclass1.insert(fn_img);
                }
                else
                {
                    avg2() += img();
                    w2 += 1.;
                    SFclass2.insert(fn_img);
                }
            }
        }
    }
	
    // Re-alignment of the class
    if (nr_iter > 0)
    {
        SFclass = SFclass1 + SFclass2;
        avg() = avg1() + avg2();
        w = w1 + w2;
        avg.set_weight(w);
        writeToDisc(avg,dirno,SFclass,fn_out,false,"ref.xmp");
        reAlignClass(avg1, avg2, SFclass1, SFclass2, 
                     exp_imgs, exp_split, exp_number, 
                     dirno, my_output);
        w1 = avg1.weight();
        w2 = avg2.weight();
    }

    // Apply Wiener filters
    if (fn_wien != "")
    {
        applyWienerFilter(avg1());
        applyWienerFilter(avg2());
    }

    // Output total and split averages and selfiles to disc
    SFclass = SFclass1 + SFclass2;
    avg() = avg1() + avg2();
    w = w1 + w2;
    avg.set_weight(w);
    avg1.set_weight(w1);
    avg2.set_weight(w2);
    writeToDisc(avg,dirno,SFclass,fn_out+"_class",!dont_write_selfiles);
    if (do_split)
    {
        writeToDisc(avg1,dirno,SFclass1,fn_out1+"_class",!dont_write_selfiles);
        writeToDisc(avg2,dirno,SFclass2,fn_out2+"_class",!dont_write_selfiles);
    }
    
    my_output[0] = (double)dirno;
    my_output[1] = w;
    my_output[2] = w1;
    my_output[3] = w2;

}

void Prog_angular_class_average_prm::writeToDisc(ImageXmipp avg,
                                                 int        dirno,
                                                 SelFile    SF,
                                                 FileName   fn,
                                                 bool       write_selfile,
                                                 FileName   oext)
{
    FileName   fn_tmp;
    double     w = avg.weight(), w_old;
    ImageXmipp old;
    SelFile    SFold;

    if (w > 0.)
    {
	avg()/=w;
        // Write class average to disc
        fn_tmp.compose(fn,dirno,oext);
        if (do_add && exists(fn_tmp) )
        {
            old.read(fn_tmp);
            w_old = old.weight();
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(old())
            {
                dMij(old(),i,j) = ( w_old * dMij(old(),i,j) + w * dMij(avg(),i,j) ) / (w_old + w);
            }
            old.set_weight(w_old + w);
            old.write(fn_tmp);
        }
        else
        {
            avg.write(fn_tmp);
        }
	// Write class selfile to disc
	if (write_selfile)
	{
	    fn_tmp.compose(fn,dirno,"sel");
            if (do_add && exists(fn_tmp) ) 
            {
                SF.merge(fn_tmp);
                SF.sort_by_filenames();
            }
            SF.write(fn_tmp);
	}
        
        if (ROUND(w) != SF.ImgNo())
        {
            std::cerr<<" w = "<<w<<" SF.ImgNo()= "<<SF.ImgNo()<<" dirno = "<<dirno<<std::endl;
            REPORT_ERROR(1,"Selfile and average weight do not correspond!");
        }
    }
}


void Prog_angular_class_average_prm::addClassAverage(int dirno,
                                                     double w,
                                                     double w1,
                                                     double w2)
{
    
    Matrix1D<double> docline(7);
    FileName fn_tmp;

    DFlib.locate(dirno);
    docline(0) = DFlib(col_rot);
    docline(1) = DFlib(col_tilt);

    if (w > 0.)
    {
        fn_tmp.compose(fn_out+"_class",dirno,"xmp");
        SFclasses.insert(fn_tmp);
        DFclasses.append_comment(fn_tmp);
        docline(5) = w;
        DFclasses.append_data_line(docline);
    }
    if (do_split)
    {
        if (w1 > 0.)
        {
            fn_tmp.compose(fn_out1+"_class",dirno,"xmp");
            SFclasses1.insert(fn_tmp);
            DFclasses1.append_comment(fn_tmp);
            docline(5) = w1;
            DFclasses1.append_data_line(docline);
        }
        if (w2 > 0.)
        {
            fn_tmp.compose(fn_out2+"_class",dirno,"xmp");
            SFclasses2.insert(fn_tmp);
            DFclasses2.append_comment(fn_tmp);
            docline(5) = w2;
            DFclasses2.append_data_line(docline);
        }
    }

}


void Prog_angular_class_average_prm::finalWriteToDisc()
{

    FileName fn_tmp;
    SelFile  auxSF;
    DocFile  auxDF;

    // Write selfiles containing all classes
    fn_tmp=fn_out+"_classes.sel";
    if (do_add && exists(fn_tmp)) SFclasses.merge(fn_tmp);
    auxSF=SFclasses.sort_by_filenames();
    auxSF.write(fn_tmp);
    if (do_split)
    {
        fn_tmp=fn_out1+"_classes.sel";
        if (do_add && exists(fn_tmp)) 
            SFclasses1.merge(fn_tmp);
        SFclasses1.write(fn_tmp);
        fn_tmp=fn_out2+"_classes.sel";
        if (do_add && exists(fn_tmp)) 
            SFclasses2.merge(fn_tmp);
        SFclasses2.write(fn_tmp);
    }

    // Write docfiles with angles and weights of all classes
    fn_tmp=fn_out+"_classes.doc";
    if (do_add && exists(fn_tmp)) 
        DFclasses.merge(fn_tmp,DOCMERGE_SUM_COLUMN,5);
    DFclasses.write(fn_tmp);
    if (do_split)
    {
        fn_tmp=fn_out1+"_classes.doc";
        if (do_add && exists(fn_tmp)) 
            DFclasses1.merge(fn_tmp,DOCMERGE_SUM_COLUMN,5);
        DFclasses1.write(fn_tmp);
        fn_tmp=fn_out2+"_classes.doc";
        if (do_add && exists(fn_tmp)) 
            DFclasses2.merge(fn_tmp,DOCMERGE_SUM_COLUMN,5);
        DFclasses2.write(fn_tmp);
    }

    // Write docfile with data for all realigned individual images
    if (nr_iter > 0)
    {
        fn_tmp=fn_out+"_realigned.doc";
        if (do_add && exists(fn_tmp)) 
        {
            // Don't do any fancy merging or sorting because those functions are really slow... 
            DF.append(fn_tmp);
            DF.remove_multiple_strings("Headerinfo");
        }
        DF.write(fn_tmp);
    }

}
