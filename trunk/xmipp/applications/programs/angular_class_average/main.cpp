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

    int              i, nmax, nr_ref, nr_images, reserve;
    double           ref_number, rot, tilt, psi, xshift, yshift, mirror;
    double           w, w1, w2;
    SelFile          SFclasses, SFclasses1, SFclasses2;
    DocFile          DFclasses, DFclasses1, DFclasses2;
    FileName         fn_tmp;
    Matrix1D<double> dataline(8), docline(7);

    // Get input parameters
    try
    {
        // Read command line & produce side info
        prm.read(argc, argv);
        prm.produceSideInfo();
        prm.show();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

    // Making class averages
    try
    {
        // Initialize
        SFclasses.clear(); 
        SFclasses1.clear();
        SFclasses2.clear();
        
        // Write docfiles with all classes
        std::string header="Headerinfo columns: rot (1) , tilt (2), psi (3), Xoff (4), Yoff (5), Weight (6), Flip (7)";
        DFclasses.append_comment(header);
        if (prm.do_split)
        {
            DFclasses1.append_comment(header);
            DFclasses2.append_comment(header);
        }

        // Reserve memory for output from class realignment
        if (prm.nr_iter > 0) reserve = prm.DF.dataLineNo();
        else reserve = 0;
        double output_values[AVG_OUPUT_SIZE*reserve+1];

        nr_ref = prm.DFlib.dataLineNo();
        init_progress_bar(nr_ref);

        // Loop over all classes
        
        for (int dirno = 1; dirno <= nr_ref; dirno++)
        {   
            // Do the actual work
            prm.processOneClass(dirno, output_values);
            
            // Output classes sel and doc files
            prm.DFlib.locate(dirno);
            docline(0) = prm.DFlib(prm.col_rot);
            docline(1) = prm.DFlib(prm.col_tilt);
            w = output_values[1];
            w1 = output_values[2];
            w2 = output_values[3];
            if (w > 0.)
            {
                fn_tmp.compose(prm.fn_out,dirno,"xmp");
                SFclasses.insert(fn_tmp);
                DFclasses.append_comment(fn_tmp);
                docline(5) = w;
                DFclasses.append_data_line(docline);
            }
            if (prm.do_split)
            {
                if (w1 > 0.)
                {
                    fn_tmp.compose(prm.fn_out1,dirno,"xmp");
                    SFclasses1.insert(fn_tmp);
                    DFclasses1.append_comment(fn_tmp);
                    docline(5) = w1;
                    DFclasses1.append_data_line(docline);
                }
                if (w2 > 0.)
                {
                    fn_tmp.compose(prm.fn_out2,dirno,"xmp");
                    SFclasses2.insert(fn_tmp);
                    DFclasses2.append_comment(fn_tmp);
                    docline(5) = w2;
                    DFclasses2.append_data_line(docline);
                }
            }

            // Fill new docfile (with params after realignment)
            if (prm.nr_iter > 0)
            {
                nr_images = ROUND(output_values[4] / AVG_OUPUT_SIZE);
                for (int i = 0; i < nr_images; i++)
                {
                    int this_image = ROUND(output_values[i*AVG_OUPUT_SIZE+5]);
                    if (!(this_image < 0))
                    {
                        prm.DF.locate(this_image);
                        prm.DF.set(0,output_values[i*AVG_OUPUT_SIZE+6]);
                        prm.DF.set(1,output_values[i*AVG_OUPUT_SIZE+7]);
                        prm.DF.set(2,output_values[i*AVG_OUPUT_SIZE+8]);
                        prm.DF.set(3,output_values[i*AVG_OUPUT_SIZE+9]);
                        prm.DF.set(4,output_values[i*AVG_OUPUT_SIZE+10]);
                        prm.DF.set(5,output_values[i*AVG_OUPUT_SIZE+11]);
                        prm.DF.set(6,output_values[i*AVG_OUPUT_SIZE+12]);
                        prm.DF.set(7,output_values[i*AVG_OUPUT_SIZE+13]);
                    }
                }
            }
            
            progress_bar(dirno);
            
        }
        progress_bar(nr_ref);
        
        // Write new document file
        if (prm.nr_iter > 0)
        {
            fn_tmp=prm.fn_out+"es_realigned.doc";
            prm.DF.write(fn_tmp);
        }
        
        // Write selfile and docfile with all classes
        fn_tmp=prm.fn_out+"es.sel";
        SFclasses.write(fn_tmp);
        fn_tmp=prm.fn_out+"es.doc";
        DFclasses.write(fn_tmp);
        if (prm.do_split)
        {
            fn_tmp=prm.fn_out1+"es.sel";
            SFclasses1.write(fn_tmp);
            fn_tmp=prm.fn_out2+"es.sel";
            SFclasses2.write(fn_tmp);
            fn_tmp=prm.fn_out1+"es.doc";
            DFclasses1.write(fn_tmp);
            fn_tmp=prm.fn_out2+"es.doc";
            DFclasses2.write(fn_tmp);
        }

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    
}
