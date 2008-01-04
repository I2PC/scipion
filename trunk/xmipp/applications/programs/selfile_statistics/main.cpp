/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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

#include <data/selfile.h>

#include <cstdlib>

/* Prototypes -============================================================= */

void Usage(char **argv);


int main(int argc, char *argv[])
{

    /* Input Parameters ======================================================== */
    FileName       sel_file;   // selection file


    /* Parameters ============================================================== */
    try
    {
        if (argc != 2)
        {
            Usage(argv);
            exit(0);
        }
        else
        {
            sel_file = argv[1];
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage(argv);
    }

    try
    {

        SelFile SF(sel_file);
        int total = 0;
        int active = 0;
        int discarded = 0;
        int commented = 0;
        while (!SF.eof())
        {
            // Get file
            if (SF.Is_ACTIVE())
                active++;
            else if (SF.Is_DISCARDED())
                discarded++;
            else if (SF.Is_COMMENT())
                commented++;
            SF.next();
            total++;
        }  // while

        std::cout << "RESULTS: " << std::endl;
        std::cout << "total number of images : " << total << " (100%)" << std::endl;
        std::cout << "active images : " << active << " (" << (float) active*100.0 / (float) total << "%)" << std::endl;
        std::cout << "discarded images : " << discarded << " (" << (float) discarded*100.0 / (float) total << "%)" << std::endl;
        std::cout << "commented images : " << commented << " (" << (float) commented*100.0 / (float) total << "%)" << std::endl;


    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    exit(0);
} //main

/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf(
        "Usage: %s [Purpose and Parameters]"
        "\nPurpose: Sel file statistics"
        "\nParameter Values: (note space before value)"
        "\nI/O parameters"
        "\n    input_file    input sel file"
        "\n  "
        "\nExample: "
        "\n    selstats c3u.sel "

        "\n"
        , argv[0]);
}

