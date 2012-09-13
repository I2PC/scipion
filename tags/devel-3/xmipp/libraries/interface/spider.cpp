/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
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

#include "spider.h"

#include <data/args.h>
#include <data/image.h>
#include <data/geometry.h>

// Generate Count File -----------------------------------------------------
void generate_Spider_count(int imax, DocFile &DF_out)
{
    MultidimArray<double>   aux(1);

    DF_out.clear();
    DF_out.append_comment((std::string)"Count for Spider up to " + integerToString(imax));

    for (aux(0) = 1; aux(0) <= imax; aux(0)++)
        DF_out.append_data_line(aux);
}

// Translate to Spider selfile ---------------------------------------------
void translate_to_Spider_sel(SelFile &SF_in, DocFile &DF_out, bool new_style)
{
    MultidimArray<double>   aux(1);
    int               selline = 1;

    DF_out.clear();
    DF_out.append_comment((std::string)"Translation for Spider of " + SF_in.name());

    SF_in.go_beginning();
    while (!SF_in.eof())
    {
        bool store = true;
        if (!SF_in.Is_COMMENT())
        {
            if (SF_in.Is_ACTIVE())
            {
                if (!new_style) aux(0) = 1;
                else            aux(0) = ((FileName)SF_in.get_current_file()).get_number();
            }
            else
            {
                if (!new_style) aux(0) = 0;
                else            store = false;
            }
            if (store) DF_out.append_data_line(aux);
        }
        SF_in.next();
    }
}

// Extract angles ----------------------------------------------------------
void extract_angles(SelFile &SF_in, DocFile &DF_out,
                    const std::string &ang1, const std::string &ang2,
                    const std::string &ang3)
{

    checkAngle(ang1);
    checkAngle(ang2);
    checkAngle(ang3);

    DF_out.clear();
    DF_out.append_comment((std::string)"Angles for " + SF_in.name() +
                          ".   Angle order: " + ang1 + " " + ang2 + " " + ang3);

    int i = 0;
    time_config();
    std::cerr << "Extracting angles ...\n";
    init_progress_bar(SF_in.ImgNo());
    while (!SF_in.eof())
    {
        // Read image
        Image<double> P;
        FileName fn_img=SF_in.NextImg();
        if (fn_img=="") break;
        P.read(fn_img);
        DF_out.append_angles(P.rot(), P.tilt(), P.psi(),
        		ang1, ang2, ang3);
        i++;
        if (i % 10 == 0) progress_bar(i);
    }
    progress_bar(SF_in.ImgNo());
}
#ifdef NEVERDEFINED
// write_angles
void write_angles(SelFile &SF_in, DocFile &DF_in,
                  const string &ang1 = "rot", const string &ang2 = "tilt",
                  const string &ang3 = "psi")
{
    double rot, tilt, psi;
    int FirstLine_colNumber;

    checkAngle(ang1);
    checkAngle(ang2);
    checkAngle(ang3);

//   std::cout << "FirstLine_colNumber" << DF_in.FirstLine_colNumber();

    int i = 0;
    time_config();
    std::cerr << "Writting new headers ...\n";
    init_progress_bar(SF_in.ImgNo());
    FirstLine_colNumber = DF_in.FirstLine_colNumber();
    while (!SF_in.eof())
    {
        // Read image
        ImageXmipp P;
        FileName fn_img=SF_in.NextImg();
        if (fn_img=="") break;
        P.read(fn_img);
        P.clear_fFlag_flag();
        if (FirstLine_colNumber >= 3)
        {
            DF_in.get_angles(i + 1, rot, tilt, psi, ang1, ang2, ang3);
            P.set_eulerAngles(rot, tilt, psi);
        }
        if (FirstLine_colNumber >= 6)
        {
            DF_in.get_angles1(i + 1, rot, tilt, psi, ang1, ang2, ang3);
            P.set_eulerAngles1(rot, tilt, psi);
        }
        if (FirstLine_colNumber >= 9)
        {
            DF_in.get_angles2(i + 1, rot, tilt, psi, ang1, ang2, ang3);
            P.set_eulerAngles2(rot, tilt, psi);
        }

        P.write(P.name());
        i++;
        if (i % 10 == 0) progress_bar(i);
    }
    progress_bar(SF_in.ImgNo());

}
#endif
// Rename for Spider -------------------------------------------------------
void rename_for_Spider(SelFile &SF_in, SelFile &SF_out, const FileName &fn_root,
                       const FileName &out_ext)
{
    FileName fn_in, fn_out;
    int counter = 1;

    SF_out.clear();
    while (!SF_in.eof())
    {
        fn_in = SF_in.NextImg();
        if (fn_in=="") break;
        fn_out = fn_root + integerToString(counter, 5);
        if (out_ext == "") fn_out = fn_out.add_extension(fn_in.get_extension());
        else             fn_out = fn_out.add_extension(out_ext);
        SF_out.insert(fn_out);

        std::cout << "Renaming " << fn_in << " as " << fn_out << std::endl;
        std::string command = (std::string)"cp " + fn_in + " " + fn_out;
        system(command.c_str());

        counter++;
    }
}


