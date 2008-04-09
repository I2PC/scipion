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

}

void Prog_angular_class_average_prm::processOneClass(int &dirno, 
                                                     double &w,
                                                     double &w1,
                                                     double &w2) {

    ImageXmipp img, avg, avg1, avg2;
    FileName   fn_img, fn_tmp;
    SelFile    SFclass, SFclass1, SFclass2;
    double     rot, tilt, psi, xshift, yshift, mirror, val;
    int        ref_number;
    int        isplit;
    Matrix2D<double> A(3,3);

    w = 0., w1 = 0., w2 = 0.;
    avg=Iempty;
    SFclass.clear();
    if (do_split)
    {
	     avg1=Iempty;
	     avg2=Iempty;
	     SFclass1.clear();
	     SFclass2.clear();
             randomize_random_generator();
    }

    // Loop over all images in the input docfile
    DFlib.locate(dirno);
	rot = DFlib(col_rot - 1);
	tilt = DFlib(col_tilt - 1);
    DF.go_beginning();
    for (int n = 0; n < DF.dataLineNo(); n++)
    {
	    DF.next();
	    if (DF.get_current_line().Is_comment()) fn_img = ((DF.get_current_line()).get_text()).erase(0, 3);
	    else  REPORT_ERROR(1, "Problem with NewXmipp-type document file");
	    DF.adjust_to_data_line();
            ref_number = ROUND(DF(col_ref - 1));
	    // Check for matching rot and tilt
    //	if (ABS(rot-lib_rot) < 0.01 && ABS(tilt-lib_tilt)<0.01)
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

		         // Apply in-plane transformation
		         A = img.get_transformation_matrix();
		         if (!A.isIdentity())
		             img().selfApplyGeometryBSpline(A, 3, IS_INV,WRAP);

		         // Add to average
		         avg() += img();
		         w+= 1.;
		         SFclass.insert(fn_img);

		         // Add to split averages
		         if (do_split)
		         {
		             isplit = ROUND(rnd_unif());
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
    }
	
    // Write output files
    if (w > 0.)
    {
	FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(avg())
	{
	    MULTIDIM_ELEM(avg(), i) /= w;
	}
	avg.clear_header();
	avg.set_eulerAngles((float)rot, (float)tilt, (float)0.);
	avg.set_originOffsets(0., 0.);
	avg.flip() = 0.;
	avg.weight() = w;
	// Write class average to disc
	fn_tmp.compose(fn_out,dirno,"xmp");
	avg.write(fn_tmp);
	// Write class selfile to disc
	if (!dont_write_selfiles)
	{
	    fn_tmp.compose(fn_out,dirno,"sel");
	    SFclass.write(fn_tmp);
	}
    }

    if (do_split)
    {
	if (w1 > 0.)
	{
	    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(avg1())
	    {
		MULTIDIM_ELEM(avg1(), i) /= w1;
	    }
	    avg1.clear_header();
	    avg1.set_eulerAngles((float)rot, (float)tilt, (float)0.);
	    avg1.set_originOffsets(0., 0.);
	    avg1.flip() = 0.;
	    avg1.weight() = w1;
	    // Write class ave1rage to disc
	    fn_tmp.compose(fn_out1,dirno,"xmp");
	    avg1.write(fn_tmp);
	    if (!dont_write_selfiles)
	    {
		// Write class selfile to disc
		fn_tmp.compose(fn_out1,dirno,"sel");
		SFclass1.write(fn_tmp);
	    }
	}
	if (w2 > 0.)
	{
	    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(avg2())
	    {
		MULTIDIM_ELEM(avg2(), i) /= w2;
	    }
	    avg2.clear_header();
	    avg2.set_eulerAngles((float)rot, (float)tilt, (float)0.);
	    avg2.set_originOffsets(0., 0.);
	    avg2.flip() = 0.;
	    avg2.weight() = w2;
	    // Write class ave1rage to disc
	    fn_tmp.compose(fn_out2,dirno,"xmp");
	    avg2.write(fn_tmp);
	    if (!dont_write_selfiles)
	    {
		// Write class selfile to disc
		fn_tmp.compose(fn_out2,dirno,"sel");
		SFclass2.write(fn_tmp);
	    }
	}
    }
}

