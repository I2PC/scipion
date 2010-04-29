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
#include <data/volume.h>
#include <cstdio>
#include <data/image.h>
#include <data/metadata.h>
void Usage();

int main(int argc, char **argv)
{
    FileName        fn_input;
    headerXmipp     header;
    MetaData SF;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = getParameter(argc, argv, "-i");
        if (Is_VolumeXmipp(fn_input) || Is_ImageXmipp(fn_input) ||
            Is_FourierVolumeXmipp(fn_input) || Is_FourierImageXmipp(fn_input))
        {
            SF.addObject();
            SF.setValue( MDL_IMAGE, fn_input);
            SF.setValue( MDL_ENABLED, 1);
        }
        else
        {
            SF.read( fn_input ,NULL);
        }
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
        long int ret=SF.firstObject();
        if(ret==MetaData::NO_OBJECTS_STORED)
        {
            std::cerr << "Empty inputFile File\n";
            exit(1);
        }
        do
        {
            FileName file_name;
            SF.getValue( MDL_IMAGE, file_name);

            if (file_name=="") break;

            // For volumes ........................................................
            if (Is_VolumeXmipp(file_name))
            {
                header.read(file_name);
                std::cout << "FileName     : " << file_name << std::endl;
                std::cout << header;

                // For images .........................................................
            }
            else if (Is_ImageXmipp(file_name))
            {

                header.read(file_name);
                std::cout << "FileName     : " << file_name << std::endl;
                std::cout << header;
                std::cout << std::endl;

                // For fourier volumes .................................................
            }
            else if (Is_FourierVolumeXmipp(file_name))
            {
                header.read(file_name);
                std::cout << "FileName     : " << file_name << std::endl;
                std::cout << header;

                // For fourier images .................................................
            }
            else if (Is_FourierImageXmipp(file_name))
            {
                header.read(file_name);
                std::cout << "FileName     : " << file_name << std::endl;
                std::cout << header;

                // Is not an Spider file ..............................................
            }
            else
                std::cout << file_name << " is not a Spider File";

            // Finish information .................................................
            std::cout << std::endl;

        } // while
        while (SF.nextObject()!= MetaData::NO_MORE_OBJECTS);

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
              << "    -i               : metaDataFile with images/volumes \n"
              << "                        or individual image or volume \n"
              ;
}

