/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.csic.es)
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

#include <reconstruction/angular_projection_matching.h>
 

int main(int argc, char **argv)
{

    MetaData                         DFo;
    Prog_angular_projection_matching_prm prm;
    FileName                         fn_tmp;
    Matrix1D<double>                 dataline(8);

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

    int nr_images = prm.DFexp.size();
    int input_images[nr_images+1];
    double output_values[MY_OUPUT_SIZE*nr_images+1];

    try
    {
        // Process all images
	    input_images[0]=nr_images;
	    for (int i = 0; i < nr_images; i++)
	    {
	        input_images[i+1]=i;
	    }
            prm.processSomeImages(input_images,output_values);

	    // Fill output docfile
	    for (int i = 0; i < nr_images; i++)
	    {
	        prm.DFexp.goToObject(ROUND(output_values[i*MY_OUPUT_SIZE+1]));
            prm.DFexp.getValue(MDL_IMAGE,fn_tmp);

            DFo.addObject();
            DFo.setValue(MDL_IMAGE,fn_tmp);
            DFo.setValue(MDL_ANGLEROT, output_values[i*MY_OUPUT_SIZE+2]);
            DFo.setValue(MDL_ANGLETILT,output_values[i*MY_OUPUT_SIZE+3]);
            DFo.setValue(MDL_ANGLEPSI, output_values[i*MY_OUPUT_SIZE+4]);
            DFo.setValue(MDL_SHIFTX,   output_values[i*MY_OUPUT_SIZE+5]);
            DFo.setValue(MDL_SHIFTY,   output_values[i*MY_OUPUT_SIZE+6]);
            DFo.setValue(MDL_REF,(int)(output_values[i*MY_OUPUT_SIZE+7]));
           	DFo.setValue(MDL_FLIP,    (output_values[i*MY_OUPUT_SIZE+8]>0));
            DFo.setValue(MDL_MAXCC,    output_values[i*MY_OUPUT_SIZE+9]);
	    }

	    fn_tmp=prm.fn_root + ".doc";
	    DFo.write(fn_tmp);
	    std::cerr<<"done!"<<std::endl;
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

}




