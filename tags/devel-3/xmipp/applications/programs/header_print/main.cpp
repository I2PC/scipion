/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
 *             Sjors Scheres
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
#include <data/selfile.h>
#include <data/image.h>
#include <cstdio>

void Usage();

int main(int argc, char **argv)
{
    FileName        fn_input;
    SelFile         SF;
    Image<double>  image;
    int             show_old_rot;     // True if old rot is to be shown

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = getParameter(argc, argv, "-i");
        if (image.isImage(fn_input) )
        {
            SF.insert(fn_input, SelLine::ACTIVE);
        }
        else
            SF.read(fn_input);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }

    try
    {

        // Process each file -----------------------------------------------------
        SF.go_beginning();
        while (!SF.eof())
        {
            FileName file_name = SF.NextImg();
            if (file_name=="") break;

            image.read(file_name, false, -1, false);
			std::cout << "FileName     : " << file_name << std::endl;
			std::cout << image;
			std::cout << std::endl;
        } // while

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }

    exit(0);

} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cerr << "Purpose:\n";
    std::cerr << "    Prints to screen some of the information stored in the header\n";

    std::cerr << "Usage: header_print " << std::endl
    << "    -i               : Selfile with images/volumes \n"
    << "                        or individual image or volume \n"
    <<std::endl;
}

