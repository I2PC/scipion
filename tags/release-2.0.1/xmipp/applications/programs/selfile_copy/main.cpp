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

#include <data/args.h>
#include <data/selfile.h>

#include <cstdlib>

void Usage(char **argv);

int main(int argc, char *argv[])
{

    FileName       sel_file;   // selection file
    string         dest_path;  // extension for output files in selection file.


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
        cout << XE;
        Usage(argv);
    }

    try
    {

        /* Perform copy or move =================================================== */

        // Finds last slash in path (to avoid copying files to another file)

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
            cout << endl << "Error, destination path is not a valid directory name " << endl;
            Usage(argv);
            exit(0);
        }


        // Finds last slash
        string org_path;
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
        SelFile SF_out;
        string comStr;
        while (!SF.eof())
        {
            // Get file
            SelLine line = SF.current();
            if (line.Is_data())
            {   //The SelLine is not a comment
                if (SF.Is_ACTIVE())
                {
                    FileName in_name = line.get_text();
                    string pathname = line.get_text();
                    string myfilename = pathname.substr(pathname.rfind("/") + 1, string::npos);
                    SF_out.insert(myfilename, SelLine::ACTIVE);
                    comStr = "cp " + org_path + in_name + " " + dest_path;
                    system(comStr.c_str()) ;
                }//if getlabel
            }//if isdata
            SF.next();
        }  // while

        // now copy sel file
        SF_out.write(sel_file.remove_directories());
        exit(0);
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }


} //main

/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf(
        "Usage: %s [Purpose and Parameters]"
        "\nPurpose: copy the images in a sel file (and the sel file) to a destination directory"
        "\nParameter Values: (note space before value)"
        "\nI/O parameters"
        "\n    input_file    input sel file"
        "\n    path path     path of destination"
        "\n  "
        "\nExample: "
        "\n    cpsel c3u.sel New_Images/ "
        "\n    (will copy all images in c3u.sel file to New_Images directory) "

        "\n"
        , argv[0]);
}

/* ------------------------------------------------------------------------- */
/* Menu                                                                      */
/* ------------------------------------------------------------------------- */
/*Colimate:
   PROGRAM CpSel {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/CpSel/Help/cpsel.html";
      help="Copy all files in a selection file to another path";
      OPEN MENU menu_selfile;
      COMMAND LINES {
        + usual: xmipp_cpsel $SELFILE_IN $PATH
      }
      PARAMETER DEFINITIONS {
         $SELFILE_IN {
            label="Sel File In";
            help="Input Selfile";
            type=file existing;
         }
         $PATH {
            label="Destination path";
            help="Where to perform operation";
            type=file;
         }
      }
   }
*/
