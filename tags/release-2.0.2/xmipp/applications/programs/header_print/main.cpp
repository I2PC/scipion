/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (1999)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/selfile.h>
#include <data/volume.h>
#include <cstdio>

void Usage();

int main(int argc, char **argv)
{
    FileName        fn_input;
    SelFile         SF;
    headerXmipp     header;
    int             show_old_rot;     // True if old rot is to be shown

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = getParameter(argc, argv, "-i");
        if (Is_VolumeXmipp(fn_input) || Is_ImageXmipp(fn_input) ||
            Is_FourierVolumeXmipp(fn_input) || Is_FourierImageXmipp(fn_input))
        {
            SF.insert(fn_input, SelLine::ACTIVE);
        }
        else
            SF.read(fn_input);

        show_old_rot = checkParameter(argc, argv, "-show_old_rot");
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
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

            // For volumes ........................................................
            if (Is_VolumeXmipp(file_name))
            {
                header.read(file_name);
                cout << "FileName     : " << file_name << endl;
                cout << header;

                // For images .........................................................
            }
            else if (Is_ImageXmipp(file_name))
            {

                header.read(file_name);
                cout << "FileName     : " << file_name << endl;
                cout << header;
                if (show_old_rot)
                    cout << "Old rot      : " << header.old_rot() << endl;
                cout << endl;

                // For fourier volumes .................................................
            }
            else if (Is_FourierVolumeXmipp(file_name))
            {
                header.read(file_name);
                cout << "FileName     : " << file_name << endl;
                cout << header;

                // For fourier images .................................................
            }
            else if (Is_FourierImageXmipp(file_name))
            {
                header.read(file_name);
                cout << "FileName     : " << file_name << endl;
                cout << header;

                // Is not an Spider file ..............................................
            }
            else
                cout << file_name << " is not a Spider File";

            // Finish information .................................................
            cout << endl;

        } // while

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }

    exit(0);

} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    cerr << "Purpose:\n";
    cerr << "    Prints to screen some of the information stored in the header\n";

    cerr << "Usage: header_print " << endl
    << "    -i               : Selfile with images/volumes \n"
    << "                        or individual image or volume \n"
    << "   [-show_old_rot]   : also show old rotational angle\n";
}

