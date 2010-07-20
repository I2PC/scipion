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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/image.h>
#include <data/metadata.h>
#include <reconstruction/angular_class_average.h>

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    Prog_angular_class_average_prm  prm;

    int              i, nmax, nr_ref, nr_images, reserve;
    double           ref_number, rot, tilt, psi, xshift, yshift, mirror;
    double           w, w1, w2;
    Matrix1D<double> dataline(8);

    // Get input parameters
    try
    {
        // Read command line & produce side info
        prm.read(argc, argv);
        prm.produceSideInfo();
        prm.show();

        // Only for do_add: append input docfile to add_to docfile
        if (prm.do_add)
        {
            FileName fn_tmp=prm.fn_out+".doc";
            if (exists(fn_tmp))
            {
                MetaData DFaux = prm.DF;
                // Don't do any fancy merging or sorting because those functions are really slow...
                DFaux.merge(fn_tmp);
                DFaux.write(fn_tmp);
            }
            else
            {
                prm.DF.write(fn_tmp);
            }
        }
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

        // Reserve memory for output from class realignment
        if (prm.nr_iter > 0)
            reserve = prm.DF.size();
        else
            reserve = 0;
        double output_values[AVG_OUPUT_SIZE*reserve+1];

        nr_ref = prm.DFlib.size();
        init_progress_bar(nr_ref);

        // Loop over all classes

        for (int dirno = 1; dirno <= nr_ref; dirno++)
        {
            // Do the actual work
            prm.processOneClass(dirno, output_values);

            // Output classes sel and doc files
            w = output_values[1];
            w1 = output_values[2];
            w2 = output_values[3];
            prm.addClassAverage(dirno,w,w1,w2);

            // Fill new docfile (with params after realignment)
            if (prm.nr_iter > 0)
            {
                nr_images = ROUND(output_values[4] / AVG_OUPUT_SIZE);
                for (int i = 0; i < nr_images; i++)
                {
                    int this_image = ROUND(output_values[i*AVG_OUPUT_SIZE+5]);
                    if (!(this_image < 0))
                    {
                        //FIXME: The next line has no sense since the MDL_IMAGE is string
                        // and 'this_image' is of type int...
                        REPORT_ERROR(-99, "The next line has no sense since the MDL_IMAGE is string \
                                     and 'this_image' is of type int...");
                        prm.DF.gotoFirstObject(MDValueEQ(MDL_IMAGE,this_image));

                        prm.DF.setValue(MDL_ANGLEROT,output_values[i*AVG_OUPUT_SIZE+6]);
                        prm.DF.setValue(MDL_ANGLETILT,output_values[i*AVG_OUPUT_SIZE+7]);
                        prm.DF.setValue(MDL_ANGLEPSI,output_values[i*AVG_OUPUT_SIZE+8]);
                        prm.DF.setValue(MDL_SHIFTX,output_values[i*AVG_OUPUT_SIZE+9]);
                        prm.DF.setValue(MDL_SHIFTY,output_values[i*AVG_OUPUT_SIZE+10]);
                        prm.DF.setValue(MDL_REF,output_values[i*AVG_OUPUT_SIZE+11]);
                        prm.DF.setValue(MDL_FLIP,output_values[i*AVG_OUPUT_SIZE+12]);
                        prm.DF.setValue(MDL_MAXCC,output_values[i*AVG_OUPUT_SIZE+13]);
                    }
                }
            }

            progress_bar(dirno);

        }
        progress_bar(nr_ref);

        // Write selfiles and docfiles with all class averages
        prm.finalWriteToDisc();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
}
