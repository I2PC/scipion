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
#include <data/geometry.h>

// Generate Count File -----------------------------------------------------
void generate_Spider_count(int imax, DocFile &DF_out)
{
    Matrix1D<double>   aux(1);

    DF_out.clear();
    DF_out.append_comment((std::string)"Count for Spider up to " + integerToString(imax));

    for (aux(0) = 1; aux(0) <= imax; aux(0)++)
        DF_out.append_data_line(aux);
}

// Translate to Spider selfile ---------------------------------------------
void translate_to_Spider_sel(MetaData &SF_in, DocFile &DF_out, bool new_style)
{
    Matrix1D<double>   aux(1);
    //int               selline = 1;

    DF_out.clear();
    DF_out.append_comment((std::string)"Translation for Spider of " + SF_in.getFilename());
    int i=1;
    FOR_ALL_OBJECTS_IN_METADATA(SF_in)
    {
        bool store = true;
        //if (!SF_in.Is_COMMENT())
        {
            int enabled;

            SF_in.getValue( MDL_ENABLED, enabled ,__iter.objId);

            if ( enabled ==1)
            {
                if (!new_style)
                    aux(0) = 1;
                //else            aux(0) = ((FileName)SF_in.get_current_file()).get_number();
                else
                    aux(0) = i++;
            }
            else
            {
                if (!new_style)
                    aux(0) = 0;
                else
                {
                    store = false;
                    i++;
                }
            }
            if (store)
                DF_out.append_data_line(aux);
        }
    }

}

// Extract angles ----------------------------------------------------------
void extract_angles(MetaData &SF_in, DocFile &DF_out,
                    const std::string &ang1, const std::string &ang2,
                    const std::string &ang3, bool fromMetadata)
{
    checkAngle(ang1);
    checkAngle(ang2);
    checkAngle(ang3);

    DF_out.clear();

    FileName auxFn;

    SF_in.getValue( MDL_IMAGE, auxFn, SF_in.firstObject() );

    DF_out.append_comment((std::string)"Angles for " + auxFn +
                          ".   Angle order: " + ang1 + " " + ang2 + " " + ang3);

    int i = 0;
    time_config();
    std::cerr << "Extracting angles ...\n";
    init_progress_bar(SF_in.size());
    Image<double> P;
    FileName fn_img;
    ApplyGeoParams params;
    params.datamode = HEADER;

    FOR_ALL_OBJECTS_IN_METADATA(SF_in)
    {
        if (fromMetadata)
        {
            double rot;
            SF_in.getValue(MDL_ANGLE_ROT,rot, __iter.objId);
            double tilt;
            SF_in.getValue(MDL_ANGLE_TILT,tilt, __iter.objId);
            double psi;
            SF_in.getValue(MDL_ANGLE_PSI,psi, __iter.objId);
            DF_out.append_angles(rot, tilt, psi,
                                 ang1, ang2, ang3);
        }
        else
        {
            // Read image
            SF_in.getValue( MDL_IMAGE, fn_img, __iter.objId);
            if (fn_img=="")
                break;
            P.readApplyGeo(fn_img,SF_in,__iter.objId, params);
            DF_out.append_angles(P.rot(), P.tilt(), P.psi(),
                                 ang1, ang2, ang3);
        }
        i++;
        if (i % 10 == 0)
            progress_bar(i);
    }

    progress_bar(SF_in.size());
}

// Rename for Spider -------------------------------------------------------
void rename_for_Spider(MetaData &SF_in, MetaData &SF_out, const FileName &fn_root,
                       const FileName &out_ext)
{
    FileName fn_in, fn_out;
    int counter = 1;
    size_t id;

    FOR_ALL_OBJECTS_IN_METADATA(SF_out)
    {
        SF_in.getValue( MDL_IMAGE, fn_in, __iter.objId);
        if (fn_in=="")
            break;
        fn_out = fn_root + integerToString(counter, 5);
        if (out_ext == "")
            fn_out = fn_out.addExtension(fn_in.getExtension());
        else
            fn_out = fn_out.addExtension(out_ext);
        id = SF_out.addObject();
        SF_out.setValue( MDL_IMAGE, fn_out, id);
        SF_out.setValue( MDL_ENABLED, 1, id);

        std::cout << "Renaming " << fn_in << " as " << fn_out << std::endl;
        std::string command = (std::string)"cp " + fn_in + " " + fn_out;
        system(command.c_str());

        counter++;
    }

}
