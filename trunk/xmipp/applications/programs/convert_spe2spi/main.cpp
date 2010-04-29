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
#include <data/selfile.h>
#include <data/image.h>

void Usage(char **argv);
void spe2spi(const FileName &, const FileName &, bool);

int main(int argc, char *argv[])
{
    FileName       fn_in;    // input file
    FileName       fn_out;   // output file

    FileName       fn_sel;   // input selfile
    FileName       fn_oext;  // output extension

    bool           reverse_endian;

    /* Parameters ============================================================== */
    try
    {
        if (argc == 1) Usage(argv);
        if (checkParameter(argc, argv, "-i"))
        {
            fn_in = getParameter(argc, argv, "-i");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "spe2spi: -i option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "-o"))
        {
            fn_out = getParameter(argc, argv, "-o");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "spe2spi: -o option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "-sel"))
        {
            fn_sel = getParameter(argc, argv, "-sel");
            fn_oext  = getParameter(argc, argv, "-oext", "xmp");
            if (checkParameter(argc, argv, "-i") || checkParameter(argc, argv, "-o"))
                EXIT_ERROR(1, "spe2spi: -sel option is not compatible with -i or -o");
        }

        reverse_endian = checkParameter(argc, argv, "-reverse_endian");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage(argv);
    }

    try
    {
        /* Perform conversion ====================================================== */

        /* input is a sel file*/
        if (fn_sel!="")
        {
            SelFile SF(fn_sel), SF_out;
            std::cerr << "Converting from SPE to SPI ...\n";
            init_progress_bar(SF.ImgNo());
            int i=0;
            while (!SF.eof())
            {
                FileName in_name = SF.NextImg();
                FileName out_name = in_name.without_extension()+"."+fn_oext;
                SF_out.insert(out_name);
                spe2spi(in_name, out_name, reverse_endian);
                progress_bar(i++);
            }
            progress_bar(SF.ImgNo());
            SF_out.write(fn_sel.without_extension()+"_spider.sel");
        }

        /* input/output are single files */
        else if (fn_in!="" && fn_out!="")
            spe2spi(fn_in, fn_out, reverse_endian);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    exit(0);
}

/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf(
        "Usage: %s [Purpose and Parameters]"
        "\nPurpose: Convert from SPE to Spider format "
        "\nParameter Values: (note space before value)"
        "\nESPECIFIC PARAMETERS FOR SINGLE-FILE CONVERSION"
        "\n    -i    file_in        input spe file"
        "\n    -o    file_out       output Spider file"
        "\nESPECIFIC PARAMETERS FOR SEL-FILE CONVERSION"
        "\n    -sel  input_file     input sel file"
        "\n    -oext input_file     extension for the output files"
        "\nGENERAL PARAMETERS"
        "\n   [-reverse_endian]     by default, output has the same endiannes as input"
        "\n"
        , argv[0]);
}

void spe2spi(const FileName &fn_in, const FileName &fn_out,
             bool reverse_endian)
{
    // Get the image size
    int Xdim, Ydim;
    short int aux;
    FILE *fh_in=fopen(fn_in.c_str(), "rb");
    if (!fh_in)
        REPORT_ERROR(1,(std::string)"Cannot find "+fn_in);
    fseek(fh_in,42,SEEK_SET);
    FREAD(&aux,2,1,fh_in,reverse_endian);
    Xdim=aux;
    fseek(fh_in,656,SEEK_SET);
    FREAD(&aux,2,1,fh_in,reverse_endian);
    Ydim=aux;
    fclose(fh_in);

    ImageXmipp Ix;
    Image *I = &Ix; // This is a trick for the compiler
    I->read(fn_in, 0, Ydim, Xdim, reverse_endian, I16, 4100);
    Ix.write(fn_out);
}
