/***************************************************************************
 *
 * Authors:    Sjors Scheres
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

#include <data/args.h>
#include <data/image.h>
#include <data/selfile.h>
#include <data/docfile.h>

void Usage();

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    int              i, nmax, col_rot = 1, col_tilt = 2, col_psi = 3;
    int              col_xshift = 4, col_yshift = 5, col_mirror = 7;
    bool             do_mirrors = false;
    float            lib_rot, lib_tilt, rot, tilt, psi, xshift, yshift, mirror;
    float            w, sumw;
    Matrix2D<double> A(3,3);
    FileName         fn_img, fn_out, fn_tmp, fn_tst;
    SelFile          SFclass, SFclasses;
    DocFile          DF, DFlib;
    ImageXmipp       img, avg;

// Check command line options ===========================================
    try
    {

        DF.read(getParameter(argc, argv, "-i"));
        DFlib.read(getParameter(argc, argv, "-lib"));
        fn_out = getParameter(argc, argv, "-o");
	fn_out += "_class";

       // Columns numbers
        if ((i = paremeterPosition(argc, argv, "-columns")) != -1)
        {
            if (i + 5 >= argc)
            {
                REPORT_ERROR(1, "Not enough integers after -columns");
            }
            col_rot = textToInteger(argv[i+1]);
            col_tilt = textToInteger(argv[i+2]);
            col_psi = textToInteger(argv[i+3]);
            col_xshift = textToInteger(argv[i+4]);
            col_yshift = textToInteger(argv[i+5]);
        }

        // Also assign weights or mirror flags?
        do_mirrors = checkParameter(argc, argv, "-mirror");
        if (do_mirrors)
            col_mirror = textToInteger(getParameter(argc, argv, "-mirror", "7"));

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage();
    }

// Making class averages ==================================================
    try
    {

	// Check that DF is of NewXmipp-type
	DF.go_beginning();
        if (DF.get_current_line().Is_comment()) fn_tst = (DF.get_current_line()).get_text();
        if (strstr(fn_tst.c_str(), "Headerinfo") == NULL)
        {
	    REPORT_ERROR(1,"ERROR: Docfile is of non-NewXmipp type.");
	}

	// Read first image to get dimensions
	DF.next();
	if (DF.get_current_line().Is_comment()) fn_img = ((DF.get_current_line()).get_text()).erase(0, 3);
	else  REPORT_ERROR(1, "Problem with NewXmipp-type document file");
	avg.read(fn_img);

	// Initialize
	int nn=DFlib.dataLineNo();
	init_progress_bar(nn);
	sumw = 0.;
	int dirno = 0;
	SFclasses.clear();
	DFlib.go_first_data_line();

	// Loop over all classes
	while (!DFlib.eof())
	{
	    dirno++;
	    DFlib.adjust_to_data_line();
	    lib_rot = DFlib(ABS(col_rot) - 1);
	    lib_tilt = DFlib(ABS(col_tilt) - 1);
	    
	    // initialize outputs
	    avg().initZeros();
	    w = 0.;
	    SFclass.clear();

	    DF.go_beginning();
	    int n = 0;
            nmax = DF.dataLineNo();
	    // Loop over all images in the input docfile
            while (n < nmax)
            {
                n++;
                DF.next();
                if (DF.get_current_line().Is_comment()) fn_img = ((DF.get_current_line()).get_text()).erase(0, 3);
                else  REPORT_ERROR(1, "Problem with NewXmipp-type document file");
		DF.adjust_to_data_line();
		rot = DF(ABS(col_rot) - 1);
		tilt = DF(ABS(col_tilt) - 1);
		if (ABS(rot-lib_rot) < 0.01 && ABS(tilt-lib_tilt)<0.01)
		{
 		    // matching rot and tilt!
		    psi    = DF(ABS(col_psi) - 1);
		    xshift = DF(ABS(col_xshift) - 1);
		    yshift = DF(ABS(col_yshift) - 1);
		    if (do_mirrors) mirror = DF(ABS(col_mirror) - 1);
		    img.read(fn_img, false, false, false, false);
		    img.set_eulerAngles((float)0., (float)0., psi);
		    img.set_originOffsets(xshift, yshift);
		    if (do_mirrors) img.flip() = mirror;

		    // Apply in-plane transformation
		    A = img.get_transformation_matrix();
		    if (!A.isIdentity())
			img().selfApplyGeometryBSpline(A, 3, IS_INV,WRAP);

		    // Add to average
		    avg() += img();
		    w+= 1.;
		    sumw += 1.;
		    SFclass.insert(fn_img);
		}
	    }
 
	    // Write output files
	    if (w > 0.)
	    {
		FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(avg())
		{
		    MULTIDIM_ELEM(avg(), i) /= w;
		}
	    }
	    avg.clear_header();
	    avg.set_eulerAngles(lib_rot, lib_tilt, (float)0.);
	    avg.set_originOffsets(0., 0.);
	    avg.flip() = 0.;
	    avg.weight() = w;
	    fn_tmp.compose(fn_out,dirno,"xmp");
	    avg.write(fn_tmp);
	    SFclasses.insert(fn_tmp);
	    fn_tmp.compose(fn_out,dirno,"sel");
	    SFclass.write(fn_tmp);
	    DFlib.next_data_line();
	    progress_bar(dirno);
	}
	progress_bar(nn);

	fn_tmp=fn_out+"es.sel";
	SFclasses.write(fn_tmp);

	if (ROUND(sumw) != nmax)
	{
	    cerr<<"sumw= "<<sumw<<" nmax= "<<nmax<<endl;
	    REPORT_ERROR(1,"ERROR: there were images in the input docfile with rot and tilt angles that were not in the library!! Dont trust the averages and selfiles!");
	}
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }

}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    printf("Purpose:\n");
    printf(" Makes class average images and corresponding selfiles from angular_projection_matching docfiles.\n");
    printf("Usage:\n");
    printf("   angular_class_averages \n");
    printf("        -i <docfile>        : docfile with assigned angles for all experimental particles\n");
    printf("        -lib <docfile>      : docfile with angles used to generate the projection matching library\n");
    printf("        -o <rootname=class> : output rootname for class averages and selfiles\n");
    printf("       [-columns] <rot=1> <tilt=2> <psi=3> <Xoff=4> <Yoff=5> \n"
           "                           : where the 5 integers are the column numbers for the \n"
           "                           : respective angles and offsets in the docfile\n"
           "                           : Note that rot & tilt are used to determine the classes \n"
           "                           : and psi, xoff & yoff are applied to calculate the class averages\n");
    printf("       [-mirror <col_m=7>] : Apply mirror operation (from docfile column col_m) (0=no-flip; 1=flip)\n");
    exit(1);
}
