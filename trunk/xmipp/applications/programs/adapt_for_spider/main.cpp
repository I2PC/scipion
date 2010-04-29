/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
#include <interface/spider.h>

void Usage();

int main(int argc, char **argv)
{
    std::string       command;
    FileName          fn_sel;
    FileName          fn_out;
    FileName          out_root;
    FileName          out_ext;
    std::string       ang1="rot",ang2="tilt",ang3="psi";
    int               maxcount;

    MetaData          SF;
    MetaData          SF_out;
    DocFile           DF;
    Matrix1D<float>   aux;
    int               selline;
    bool              newsel_style;
    bool              from_metadata;

// Check command line ------------------------------------------------------
    try
    {
        if (argc<2) REPORT_ERROR(1,"Adapt for Spider: Not enough parameters");
        command=argv[1];

        fn_out=getParameter(argc, argv, "-o");
        if        (command=="rename")
        {
            fn_sel=getParameter(argc, argv, "-i");
            out_root=getParameter(argc, argv, "-oroot");
            out_ext=getParameter(argc, argv, "-oext","");
        }
        else if (command=="translate_sel")
        {
            fn_sel=getParameter(argc, argv, "-i");
            newsel_style=checkParameter(argc,argv,"-new_style");
        }
        else if (command=="generate_count")
        {
            maxcount=textToInteger(getParameter(argc, argv, "-max"));
        }
        else if (command=="extract_angles")
        {
            int i;
            fn_sel=getParameter(argc, argv, "-i");
            if ((i=paremeterPosition(argc,argv,"-order"))!=-1)
            {
                if (i+3>=argc)
                    REPORT_ERROR(1,"Adapt for Spider: Not enough parameters behind -ang\n");
                ang1=argv[i+1];
                ang2=argv[i+2];
                ang3=argv[i+3];
            }
            from_metadata=checkParameter(argc,argv,"-from_metadata");
        }
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE;
        Usage();
        exit(1);
    }

    try
    {

// Perform operation
        if (command=="rename")
        {
            SF.read(fn_sel,NULL);
            rename_for_Spider(SF,SF_out,out_root,out_ext);
            SF_out.write(fn_out);
        }
        else if (command=="generate_count")
        {
            //saves old doc file for spider
            generate_Spider_count(maxcount,DF);
            DF.write(fn_out);
        }
        else if (command=="translate_sel")
        {
            SF.read(fn_sel,NULL);
            translate_to_Spider_sel(SF,DF,newsel_style);
            DF.write(fn_out);
        }
        else if (command=="extract_angles")
        {
            SF.read(fn_sel,NULL);
            extract_angles(SF,DF,ang1,ang2,ang3,from_metadata);
            DF.write(fn_out);
        }
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE;
    }
    exit(0);
}

// Usage -------------------------------------------------------------------
void Usage()
{
    std::cerr << "Usage: adapt_for_Spider rename            : Generate correlative names\n"
              << "            -i metaData File              : Input Xmipp selection file\n"
              << "            -o metaData file              : Output Xmipp SelFile\n"
              << "            -oroot <root_name>            : root name for output images\n"
              << "            [-oext <output extension="">  : if nothing is provided the same\n"
              << "                                            as the original images' one is used\n"
              << "       adapt_for_Spider generate_count    : Generate corresponding Spider SelFile\n"
              << "            -max <max_count>              : number images\n"
              << "            -o <DocFile>                  : Output Spider CountFile\n"
              << "       adapt_for_Spider translate_sel     : Generate corresponding Spider SelFile\n"
              << "            -i metadata file              : Input Xmipp selection file\n"
              << "            -o metadata file              : Output Spider SelFile\n"
              << "           [-new_style]                   : Generate new Spider Selfile style\n"
              << "       adapt_for_Spider extract_angles    : Generate a Docfile with angles\n"
              << "            -i metadata                   : Input Xmipp selection file\n"
              << "            -o <Ang DocFile>              : Output Docfile\n"
              << "           [-order <ang1> <ang2> <ang3>   : order of the angles\n"
              << "                                            by default, psi, tilt, rot\n"
              << "           [-from_metadata]               : get angles from metadata\n"
              << "                                            instead of from the header\n"
              ;
}

