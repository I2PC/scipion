/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)
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

#include <reconstruction/angular_projection_matching.h>
 

int main(int argc, char **argv)
{

    DocFile                          DFo;
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

    int nr_images = prm.DFexp.dataLineNo();
    int                              input_images[nr_images+1];
    double                           output_values[MY_OUPUT_SIZE*nr_images+1];

    try
    {

        DFo.clear();
        DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Refno (6), Flip (7), maxCC (8)");

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

	    prm.DFexp.locate(round(output_values[i*MY_OUPUT_SIZE+1]+1));
	    prm.DFexp.previous();
	    fn_tmp = ((prm.DFexp.get_current_line()).get_text()).erase(0, 3);

	    dataline(0)=output_values[i*MY_OUPUT_SIZE+2];
	    dataline(1)=output_values[i*MY_OUPUT_SIZE+3];
	    dataline(2)=output_values[i*MY_OUPUT_SIZE+4];
	    dataline(3)=output_values[i*MY_OUPUT_SIZE+5];
	    dataline(4)=output_values[i*MY_OUPUT_SIZE+6];
	    dataline(5)=output_values[i*MY_OUPUT_SIZE+7] + 1;
	    dataline(6)=output_values[i*MY_OUPUT_SIZE+8];
	    dataline(7)=output_values[i*MY_OUPUT_SIZE+9];
	    DFo.append_comment(fn_tmp);
	    DFo.append_data_line(dataline);
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




