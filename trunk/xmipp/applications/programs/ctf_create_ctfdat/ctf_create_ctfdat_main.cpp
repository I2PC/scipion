/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <data/metadata.h>
#include <data/ctf.h>

int main(int argc, char **argv)
{
    FileName fn_sel, fn_doc, fn_ctf, fn_param, fn_out;
    bool do_mode2;

    try
    {
        fn_sel = getParameter(argc, argv, "-i");
        fn_out = getParameter(argc, argv, "-o","out");
        do_mode2 = checkParameter(argc, argv, "-ctfs");
        if (do_mode2)
        {
            fn_ctf = getParameter(argc, argv, "-ctfs");
        }
        else
        {
            fn_doc = getParameter(argc, argv, "-doc");
            fn_param = getParameter(argc, argv, "-param");
        }
    }
    catch (XmippError XE)
    {
        std::cerr << XE << std::endl;
        std::cerr << "Usage: docfile_select_subset\n"
        << "   -i <selfile>          : Input selfile of selfiles for each micrograph\n"
        << "   -o <rootname=\"out\">   : Root name for output files \n"
        << " MODE 1: \n"
        << "   -param <param file>   : CTFparam file with common parameters \n"
        << "   -doc <docfile>        : Docfile with defocus values for each micrograph \n "
        << "                           this file may have either a single column (defocus)\n"
        << "                           or three columns (defocusU, defocusV, azimuth_angle)\n"
        << " MODE 2: \n"
        << "   -ctfs <selfile>       : Selfile of CTF param files for each micropgrah \n"
        ;
        return 1;
    }

    try
    {
        MetaData SFsel, SFind, SFctf, ctfdat;
        FileName fnsel, fnimg, fnctf;

        SFsel.read(fn_sel);

        if (do_mode2)
        {
            SFctf.read(fn_ctf);
        }
        else
        {
            // Write param files for each micrograph to disc and make an internal SFctf
            CTFDescription ctf;
            MetaData DFdef;
            double defU, defV, azi;
            int ii = 0;
            DFdef.read(fn_doc);
            ctf.read(fn_param);
            ctf.enable_CTF = true;
            ctf.enable_CTFnoise = false;

            if (SFsel.size() != DFdef.size())
                REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                             "Selfile -1 and docfile -doc have unequal number of entries! ");

            FOR_ALL_OBJECTS_IN_METADATA(DFdef)
            {
                ii++;
                DFdef.getValue(MDL_CTF_DEFOCUSU,defU);
                if (!DFdef.getValue(MDL_CTF_DEFOCUSU,defV))
                    defV=defU;
                if (!DFdef.getValue(MDL_CTF_DEFOCUS_ANGLE,azi))
                    azi=0;

                ctf.DeltafU=defU;
                ctf.DeltafV=defV;
                ctf.azimuthal_angle=azi;
                fnctf.compose(fn_out,ii,"ctfparam");
                ctf.write(fnctf);
                SFctf.addObject();
                SFctf.setValue(MDL_CTFMODEL,fnctf);
                std::cerr<<" Saved CTF parameter file "<<fnctf<<" for micrograph number "
                <<ii<<std::endl;
            }
        }

        // For both modes
        if (SFsel.size() != SFctf.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "Selfiles of options -i and -ctfs have unequal number of entries! ");
        SFctf.firstObject();
        FOR_ALL_OBJECTS_IN_METADATA(SFsel)
        {
            SFsel.getValue(MDL_SELFILE,fnsel);
            SFctf.getValue(MDL_IMAGE,fnctf);
            SFind.read(fnsel);
            FOR_ALL_OBJECTS_IN_METADATA(SFind)
            {
                SFind.getValue(MDL_IMAGE,fnimg);
                ctfdat.addObject();
                ctfdat.setValue(MDL_IMAGE,fnimg);
                ctfdat.setValue(MDL_CTFMODEL,fnctf);
            }
        }

        ctfdat.write(fn_out+".ctfdat");
        std::cerr<<" Saved CTFdat file as "<<fn_out+".ctfdat"<<std::endl;
        std::cerr<< " Done! "<<std::endl;
    }
    catch (XmippError XE)
    {
        std::cerr << XE << std::endl;
        return 2;
    }
    return 0;
}
