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
#include <reconstruction/angular_class_average.h>

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    Prog_angular_class_average_prm  prm;

    int              i, nmax, nr_ref;
    double           ref_number, rot, tilt, psi, xshift, yshift, mirror;
    double           w, w1, w2;
    SelFile          SFclasses, SFclasses1, SFclasses2;
    FileName         fn_tmp;

    // Get input parameters
    try
    {
        // Read command line & produce side info
        prm.read(argc, argv); //all nodes
        prm.show(); // only rank =0

        // Project reference volume etc.
        prm.produceSideInfo();// prerun, ONCE PER WORKING NODE

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage(); //only rank = 0 
        exit(0);
    }

    // Making class averages
    try
    {

	    // Initialize
	    SFclasses.clear(); 
	    SFclasses1.clear();
	    SFclasses2.clear();

	    nr_ref = prm.DFlib.dataLineNo();// master//library docfile
	    init_progress_bar(nr_ref);// master

	    // Loop over all classes
	    //prm.DFlib.go_beginning();// master
	    for (int dirno = 1; dirno <= nr_ref; dirno++)
	    {   //col_rot colum numbers
	        //prm.DFlib.adjust_to_data_line();
	        prm.processOneClass(dirno, w, w1, w2);

	        if (w > 0.)
	        {
		        fn_tmp.compose(prm.fn_out,dirno,"xmp");
		        SFclasses.insert(fn_tmp);
	        }
	        if (prm.do_split)
	        {
		        if (w1 > 0.)
		        {
		            fn_tmp.compose(prm.fn_out1,dirno,"xmp");
		            SFclasses1.insert(fn_tmp);
		        }
		        if (w2 > 0.)
		        {
		            fn_tmp.compose(prm.fn_out2,dirno,"xmp");
		            SFclasses2.insert(fn_tmp);
		        }
	        }
                progress_bar(dirno);
	        //prm.DFlib.next();
            }
            progress_bar(nr_ref);
 
        //maSTER SHOULD DO THE FOLLOWING...
	    // Write selfile with all classes
	    fn_tmp=prm.fn_out+"es.sel";
	    SFclasses.write(fn_tmp);
	    if (prm.do_split)
	    {
	        fn_tmp=prm.fn_out1+"es.sel";
	        SFclasses1.write(fn_tmp);
	        fn_tmp=prm.fn_out2+"es.sel";
	        SFclasses2.write(fn_tmp);
	    }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }

}
