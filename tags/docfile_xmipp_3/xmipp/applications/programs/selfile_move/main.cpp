/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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

#include <cstdlib>


void Usage(char **argv);


int main(int argc, char *argv[])
{

    FileName       sel_file;   // selection file
    std::string    dest_path;  // extension for output files in selection file.


    /* Parameters ============================================================== */
    try
    {
        if (argc != 3)
        {
            Usage(argv);
            exit(0);
        }
        else
        {
            sel_file = argv[1];
            dest_path = argv[2];
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage(argv);
    }

    try
    {

        /* Perform copy or move =================================================== */

        // Finds last slash in path (to avoid moving files to another file)

        int break_point = -1;
        for (int i = dest_path.size() - 1; i >= 0; i--)
        {
            if (dest_path[i] == '/')
            {
                break_point = i;
                break;
            }
        }

        if (break_point < 0)
        {
            std::cout << std::endl << "Error, destination path is not a valid directory name " << std::endl;
            Usage(argv);
            exit(0);
        }


        // Finds last slash in sel name
        std::string org_path;
        break_point = -1;
        for (int i = sel_file.size() - 1; i >= 0; i--)
        {
            if (sel_file[i] == '/')
            {
                break_point = i;
                break;
            }
        }

        // Copy only the path
        if (break_point >= 0)
        {
            org_path.resize(break_point + 1);
            for (int j = 0; j <= break_point; j++)
                org_path[j] = sel_file[j];
        }

        SelFile SF(sel_file);
        std::string comStr;
        while (!SF.eof())
        {
            // Get file
            SelLine line = SF.current();
            if (line.Is_data())
            {   //The SelLine is not a comment
                FileName in_name = line.get_text();
                comStr = "mv " + org_path + in_name + " " + dest_path;

                if (!system(comStr.c_str()))
                    std::cout << " file " << org_path << in_name << " moved to " << dest_path << std::endl;
            }
            SF.next();
        }  // while

        // now move sel file
        comStr = "mv " + sel_file + " " + dest_path;
        system(comStr.c_str());

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
        "\nPurpose: Moves the images in a sel file (and the sel file) to a destination directory"
        "\nParameter Values: (note space before value)"
        "\nI/O parameters"
        "\n    input_file    input sel file"
        "\n    path path     path of destination"
        "\n  "
        "\nExample: "
        "\n    mvsel c3u.sel New_Images/ "
        "\n    (will move all images in c3u.sel file to New_Images directory) "

        "\n"
        , argv[0]);
}

/* ------------------------------------------------------------------------- */
/* Menu                                                                      */
/* ------------------------------------------------------------------------- */
/*Colimate:
   PROGRAM MvSel {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/MvSel/Help/mvsel.html";
      help="Move all files in a selection file to another path";
      OPEN MENU menu_selfile;
      COMMAND LINES {
        + usual: xmipp_mvsel $SELFILE_IN $PATH
      }
   }
*/
