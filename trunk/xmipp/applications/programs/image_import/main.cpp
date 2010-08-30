/***************************************************************************
 *
 * Authors:    Jose Ramón Macias & Roberto Marabini
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

void Usage();

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    FileName        fn_in, fn_out;
    Image<double>   img;
    std::string taskNumber, itemNumber;
    std::ofstream filestr;

    // Check command line options ===========================================
    try
    {
        fn_in            = getParameter(argc, argv, "-inputFileName");
        taskNumber       = getParameter(argc, argv, "-taskNumber");
        itemNumber       = getParameter(argc, argv, "-itemNumber");
        fn_out           = getParameter(argc, argv, "-outputFileName");
        int x,y,z,n;
        double sampling_rate;
        img.read(fn_in,false);
        img.getDimensions(x,y,z,n);
        sampling_rate = (double) img.samplingRateX();
        filestr.open (fn_out.c_str());
        filestr  << "<ITEM id=\"" << itemNumber <<"\">\n"
        << "<MDL_IMAGE>" << fn_in << "</MDL_IMAGE>\n"
        << "<MDL_X>" << x << "</MDL_X>\n"
        << "<MDL_Y>" << y << "</MDL_Y>\n"
        << "<MDL_Z>" << z << "</MDL_Z>\n"
        << "<MDL_N>" << n << "</MDL_N>\n"
        << "</ITEM>\n";
        filestr.close();
    }
    catch (XmippError XE)
    {
        filestr.open (fn_out.c_str());
        filestr  << "<ITEM id=\"" << itemNumber <<"\">\n"
        << "<ERROR>" << XE << "</ERROR>\n"
        << "</ITEM>\n";
        filestr.close();
        std::cout << XE;
        Usage();
    }
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cout << " Purpose:\n";
    std::cout << " Extract size and sampling rate from an image file.\n";
    std::cout << " Usage:\n";
    std::cout << "    image_import \n";
    std::cout << "         -inputFileName inputFileName\n";
    std::cout << "         -outputFileName outputFileName\n";
    std::cout << "         -taskNumber taskNumber\n";
    std::cout << "         -itemNumber itemNumber\n";
    exit(1);
}
