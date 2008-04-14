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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "angular_class_average.h"

// Read arguments ==========================================================
void Prog_angular_class_average_prm::read(int argc, char **argv)  {

    // Read command line
    DF.read(getParameter(argc, argv, "-i"));
    DFlib.read(getParameter(argc, argv, "-lib"));
    fn_out = getParameter(argc, argv, "-o");

    // Columns numbers
    int i;
    if ((i = paremeterPosition(argc, argv, "-columns")) != -1)
    {
	if (i + 6 >= argc)
	{
	    REPORT_ERROR(1, "Not enough integers after -columns");
	}
	col_rot = textToInteger(argv[i+1]);
	col_tilt = textToInteger(argv[i+2]);
	col_psi = textToInteger(argv[i+3]);
	col_xshift = textToInteger(argv[i+4]);
	col_yshift = textToInteger(argv[i+5]);
	col_ref  = textToInteger(argv[i+6]);
    }
    else
    {
	col_rot    = 1;
	col_tilt   = 2;
	col_psi    = 3;
	col_xshift = 4;
	col_yshift = 5;
	col_ref    = 6;
    }

    // TODO: INSTEAD OF FIXED VALUE TO DISCARD IMAGES, DISCARD BOTTOM X% FOR EACH CLASS
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
    col_select = textToInteger(getParameter(argc, argv, "-select", "8"));
    

    // Also assign weights or mirror flags?
    do_mirrors = checkParameter(argc, argv, "-mirror");
    if (do_mirrors)
	col_mirror = textToInteger(getParameter(argc, argv, "-mirror", "7"));

    // Perform splitting of the data?
    do_split = checkParameter(argc, argv, "-split"); 

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
    std::cerr << "  Output rootname         : "<< fn_out<<std::endl;
    if (do_split)
	std::cerr << "     -> Split data in random halves and output class averages "<<std::endl;
    if (do_mirrors)
	std::cerr << "     -> Take mirror operation into account "<<std::endl;
    if (do_limit0)
	std::cerr << "     -> Discard images with value in column "<<col_select<<" below "<<limit0<<std::endl;
    if (do_limitF)
	std::cerr << "     -> Discard images with value in column "<<col_select<<" above "<<limitF<<std::endl;
    if (dont_write_selfiles)
	std::cerr << "     -> Do not write class selfiles to disc "<<std::endl;
    if (nr_iter > 0)
    {
        std::cerr << "     -> Re-align class averages in "<<nr_iter<<" iterations"<<std::endl;
        std::cerr << "  Maximum shift           : "<<max_shift<<std::endl;
        std::cerr << "  Maximum shift change    : "<<max_shift_change<<std::endl;
        std::cerr << "  Maximum psi change      : "<<max_psi_change<<std::endl;
        if (Ri>0)
            std::cerr << "  Inner radius rot-search : "<<Ri<<std::endl;
        if (Ro>0)
            std::cerr << "  Outer radius rot-search : "<<Ro<<std::endl;
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
    printf("        -o <rootname=class> : output rootname for class averages and selfiles\n");
    printf("        -split              : also output averages of random halves of the data\n");
    printf("       [-columns] <rot=1> <tilt=2> <psi=3> <Xoff=4> <Yoff=5> "
           "<ref=6> \n"
           "                           : where the 6 integers are the column numbers for the \n"
           "                           : respective angles and offsets and optimal reference \n"
           "                           : number in the docfile\n"
           "                           : Note that rot & tilt are used to determine the classes \n"
           "                           : and psi, xoff & yoff are applied to calculate the class averages\n");
    printf("       [-mirror <col_m=7>] : Apply mirror operation (from docfile column col_m) (0=no-flip; 1=flip)\n");
    printf("       [-select <col_s=8>] : Column number to use for limit0/F selection\n");
    printf("       [-limit0 <limit0>]  : Values in column <col_s> below this are discarded\n");
    printf("       [-limitF <limitF>]  : Values in column <col_s> above this are discarded\n");
    printf("       [-dont_write_selfiles]  : Do not write class selfiles to disc\n");
    printf(" REALIGNMENT OF CLASSES \n");
    printf("       [-iter <int=0>]                  : Number of iterations for re-alignment\n");
    printf("       [-Ri <int=1>]                    : Inner radius to limit rotational search\n");
    printf("       [-Ro <int=dim/2-1>]              : Outer radius to limit rotational search\n");
    printf("       [-max_shift <float=999.>]        : Maximum shift (larger shifts will be set to 0)\n");
    printf("       [-max_shift_change <float=999.>] : Maximum change in shift in last iteration \n");
    printf("       [-max_psi_change <float=360.>]   : Maximum change in rotation in last iteration \n");
    exit(1);
}

// Side info stuff ===================================================================
void Prog_angular_class_average_prm::produceSideInfo() {

    
    FileName fn_tst, fn_img;

    if (do_split)
    {
	fn_out1 = fn_out+"_split_1_class";
	fn_out2 = fn_out+"_split_2_class";
    }
    fn_out += "_class";

    // Check that DF is of NewXmipp-type
    DF.go_beginning();
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

    // Set ring defaults
    if (Ri<1) Ri=1;
    if (Ro<0) Ro=(XSIZE(Iempty())/2)-1;

    // Randomization
    if (do_split) randomize_random_generator();

}

void Prog_angular_class_average_prm::getPolar(Matrix2D<double> &img, Polar<std::complex <double> > &fP, 
                                              bool conjugated, float xoff, float yoff)
{
    Matrix2D<double> Maux;
    Polar<double> P;

    // Calculate FTs of polar rings and its stddev
    img.produceSplineCoefficients(Maux,3);
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro,xoff,yoff);
    fP = P.fourierTransformRings(conjugated);
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
    Matrix1D<double>                ang, corr;
    Matrix2D<double>                Mimg, Mref;
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
            rotationalCorrelation(fPimg,fPref,ang,corr);
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
            rotationalCorrelation(fPimg,fPrefm,ang,corr);
	    for (int k = 0; k < XSIZE(corr); k++)
	    {
		if (corr(k)> maxcorr)
		{
		    maxcorr = corr(k);
		    opt_psi = realWRAP(360. - ang(k), -180., 180.);
                    opt_flip = 1.;
		}
	    }

#ifdef DEBUG
            std::cout<<" imgno= "<<imgno<<" Psi()= "<<imgs[imgno].Psi()<<" opt_psi= "<<opt_psi<<" flip()= "<<imgs[imgno].flip()<<" opt_flip= "<<opt_flip<<std::endl;
#endif      
            // Check max_psi_change in last iteration
            if (iter == nr_iter - 1)
            {
                diff_psi = ABS(realWRAP(imgs[imgno].Psi() - opt_psi, -180., 180.));
                if (diff_psi > max_psi_change) 
                    do_discard = true;
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
                    // STILL CHECK THIS!!!
                    best_shift(Mref,Mimg,opt_xoff,opt_yoff);
                    if (opt_xoff * opt_xoff + opt_yoff * opt_yoff > max_shift * max_shift) 
                    {
                        opt_xoff = opt_yoff = 0.;
                    }
                    else
                    {
                        Mimg.selfTranslate(vectorR2(opt_xoff,opt_yoff),true);
                    }
                    new_xoff =  opt_xoff*COSD(opt_psi) + opt_yoff*SIND(opt_psi);
                    new_yoff = -opt_xoff*SIND(opt_psi) + opt_yoff*COSD(opt_psi);

#ifdef DEBUG
     std::cout<<" imgno= "<<imgno<<" Xoff()= "<<imgs[imgno].Xoff()<<" opt_xoff= "<<opt_xoff<<" new_xoff= "<<new_xoff<<std::endl;
     std::cout<<" imgno= "<<imgno<<" Yoff()= "<<imgs[imgno].Yoff()<<" opt_yoff= "<<opt_yoff<<" new_yoff= "<<new_yoff<<std::endl;
#endif      
                    // Check max_shift_change in last iteration
                    if (iter == nr_iter - 1)
                    {
                        diff_shift = (new_xoff - opt_xoff) * (new_xoff - opt_xoff) +
                                     (new_yoff - opt_yoff) * (new_yoff - opt_yoff);
                        if (diff_shift > max_shift_change * max_shift_change) 
                            do_discard = true;
                    }
                }
            }

            if (!do_discard)
            {
                imgs[imgno].Psi() =  opt_psi;
                imgs[imgno].flip() = opt_flip;
                if (opt_flip==1.)                
                    imgs[imgno].Xoff() = -new_xoff;
                else                
                    imgs[imgno].Xoff() = new_xoff;
                imgs[imgno].Yoff() = new_yoff;
                ccfs[imgno] = correlation_index(Mref,Mimg);

 #ifdef DEBUG
                std::cout<<" imgno= "<<imgno<<" Psi()= "<<imgs[imgno].Psi()<<" Xoff()= "<<imgs[imgno].Xoff()<<" Yoff()= "<<imgs[imgno].Yoff()<<" flip()= "<<imgs[imgno].flip()<<" ccf= "<<ccfs[imgno]<<std::endl;
#endif      
                    // Check max_shift_change in last iteration
               // Add to averages
                if (splits[imgno]==0)
                {
                    w1 += 1.;
                    avg1() += Mimg;
                }
                else if (splits[imgno]==1)
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
        It() = Mref;
        It.weight()=w1+w2;
        It.write("ref.xmp");
#endif
    }

    avg1.weight() = w1;
    avg2.weight() = w2;
    
    // Report the new angles, offsets and selfiles
    my_output[4] = imgs.size() * AVG_OUPUT_SIZE;
    for (int imgno = 0; imgno < imgs.size(); imgno++)
    {
        my_output[imgno * AVG_OUPUT_SIZE + 5] = numbers[imgno];
        my_output[imgno * AVG_OUPUT_SIZE + 6] = avg1.rot();
        my_output[imgno * AVG_OUPUT_SIZE + 7] = avg1.tilt();
        my_output[imgno * AVG_OUPUT_SIZE + 8] = imgs[imgno].psi();
        //if (imgs[imgno].flip()==1)
        //    my_output[imgno * AVG_OUPUT_SIZE + 9] = -imgs[imgno].Xoff();
        //else
            my_output[imgno * AVG_OUPUT_SIZE + 9] = imgs[imgno].Xoff();
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

void Prog_angular_class_average_prm::processOneClass(int &dirno, 
                                                     double * my_output) {

    ImageXmipp img, avg, avg1, avg2;
    FileName   fn_img, fn_tmp;
    SelFile    SFclass, SFclass1, SFclass2;
    double     rot, tilt, psi, xshift, yshift, mirror, val, w, w1, w2;
    int        ref_number, this_image;
    int        isplit;
    Matrix2D<double> A(3,3);
    std::vector<int> exp_number, exp_split;
    std::vector<ImageXmipp> exp_imgs;
    
    
    // Get reference angles and preset to averages
    DFlib.locate(dirno);
    rot = DFlib(col_rot - 1);
    tilt = DFlib(col_tilt - 1);
    Iempty.set_eulerAngles((float)rot, (float)tilt, (float)0.);
    Iempty.set_originOffsets(0., 0.);
    Iempty.flip() = 0.;
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
        ref_number = ROUND(DF(col_ref - 1));
        if (ref_number == dirno)
        {
            bool is_select = true;
            val = DF(col_select - 1);
            if ( (do_limit0 && val < limit0) || (do_limitF && val > limitF) ) is_select = false;
            if (is_select)
            {
                psi    = DF(col_psi - 1);
                xshift = DF(col_xshift - 1);
                yshift = DF(col_yshift - 1);
                if (do_mirrors) mirror = DF(col_mirror - 1);
                img.read(fn_img, false, false, false, false);
                img().setXmippOrigin();
                img.set_eulerAngles((float)0., (float)0., (float)psi);
                img.set_originOffsets(xshift, yshift);
                if (do_mirrors) img.flip() = mirror;
                
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
        avg.weight() = w;
        writeToDisc(avg,dirno,SFclass,fn_out,false,"ref.xmp");
        reAlignClass(avg1, avg2, SFclass1, SFclass2, 
                     exp_imgs, exp_split, exp_number, 
                     dirno, my_output);
        w1 = avg1.weight();
        w2 = avg2.weight();
    }

    // Output total and split averages and selfiles to disc
    SFclass = SFclass1 + SFclass2;
    avg() = avg1() + avg2();
    w = w1 + w2;
    avg.weight() = w;
    avg1.weight() = w1;
    avg2.weight() = w2;
    writeToDisc(avg,dirno,SFclass,fn_out,!dont_write_selfiles);
    if (do_split)
    {
        writeToDisc(avg1,dirno,SFclass1,fn_out1,!dont_write_selfiles);
        writeToDisc(avg2,dirno,SFclass2,fn_out2,!dont_write_selfiles);
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
    FileName fn_tmp;
    double w = avg.weight();

    if (w > 0.)
    {
	avg()/=w;
	// Write class average to disc
	fn_tmp.compose(fn,dirno,oext);
	avg.write(fn_tmp);
	// Write class selfile to disc
	if (write_selfile)
	{
	    fn_tmp.compose(fn,dirno,"sel");
	    SF.write(fn_tmp);
	}
        
        if (ROUND(w) != SF.ImgNo())
        {
            std::cerr<<" w = "<<w<<" SF.ImgNo()= "<<SF.ImgNo()<<" dirno = "<<dirno<<std::endl;
            REPORT_ERROR(1,"Selfile and average weight do not correspond!");
        }
    }
}
