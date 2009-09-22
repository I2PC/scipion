/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

void Usage(char **argv);
void raw22spi(const FileName &, const FileName &, char, bool,
              int, int, int, int, bool, bool);

int main(int argc, char *argv[])
{
    FileName       fn_in;    // input file
    FileName       fn_out;   // output file
    FileName       sel_file; // selection file
    std::string    sel_ext;  // extension for output files in selection file.
    char           raw_type = 'b';
    bool           force_byte=false;
    int            Zdim, Ydim, Xdim;
    int            size_arg;
    int            header_size;
    bool           generate_inf;
    bool           reverse_endian;

    Zdim = Ydim = Xdim = 0;
    /* Parameters ============================================================== */
    try
    {
        if (checkParameter(argc, argv, "-i"))
        {
            fn_in = getParameter(argc, argv, "-i");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
            {
                EXIT_ERROR(1, "Raw22spi: -i option is not compatible with -sel or -oext");
            }
        }
        if (checkParameter(argc, argv, "-o"))
        {
            fn_out = getParameter(argc, argv, "-o");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
            {
                EXIT_ERROR(1, "Raw22spi: -o option is not compatible with -sel or -oext");
            }
        }
        if (checkParameter(argc, argv, "-sel"))
        {
            sel_file = getParameter(argc, argv, "-sel");
            sel_ext  = getParameter(argc, argv, "-oext");
            if (checkParameter(argc, argv, "-i") || checkParameter(argc, argv, "-o"))
            {
                /*error cause -sel is not compatible with -i -o*/
                EXIT_ERROR(1, "Raw22spi: -sel option is not compatible with -i or -o");
            }
        }

        if (checkParameter(argc, argv, "-f"))  raw_type = 'f';
        if (checkParameter(argc, argv, "-16")) raw_type = 'h';
        force_byte=checkParameter(argc, argv, "-8");

        if (checkParameter(argc, argv, "-s"))
        {
            size_arg = paremeterPosition(argc, argv, "-s");
            if (size_arg + 3 >= argc) EXIT_ERROR(1, "Not enough parameters behind -s");
            Zdim = textToInteger(argv[size_arg+1]);
            Ydim = textToInteger(argv[size_arg+2]);
            Xdim = textToInteger(argv[size_arg+3]);
        }
        generate_inf = checkParameter(argc, argv, "-generate_inf");
        reverse_endian = checkParameter(argc, argv, "-reverse_endian");
        if (checkParameter(argc, argv, "-is_micrograph"))
        {
            FileName fn_inf = fn_in + ".inf";
            FILE *fh_inf = fopen(fn_inf.c_str(), "r");
            if (!fh_inf)
                REPORT_ERROR(1, (std::string)"Raw22Spi:Cannot find " + fn_inf);
            Zdim = 1;
            Ydim = textToInteger(getParameter(fh_inf, "Ydim"));
            Xdim = textToInteger(getParameter(fh_inf, "Xdim"));
            fclose(fh_inf);
        }
        header_size = textToInteger(getParameter(argc, argv, "-header", "0"));
        if (argc == 1)
        {
            Usage(argv);
        }
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
        if (checkParameter(argc, argv, "-sel") && checkParameter(argc, argv, "-oext"))
        {
            SelFile SF(sel_file);
            FileName SF_out_name;
            SF_out_name = sel_file.without_extension().add_prefix("out_");
            SF_out_name += (std::string)"." + sel_file.get_extension();
            SelFile SF_out;
            SF_out.clear();
            while (!SF.eof())
            {
                SelLine line = SF.current();
                if (line.Is_data())
                { //The SelLine is not a comment
                    FileName in_name = line.get_text();
                    short label = line.get_label();
                    FileName out_name = in_name.without_extension();
                    out_name = out_name.add_extension(sel_ext);
                    SF_out.insert(out_name, (SelLine::Label)label);
                    raw22spi(in_name, out_name, raw_type, force_byte,
                        header_size, Zdim, Ydim, Xdim,
                        generate_inf, reverse_endian);
                }
                else if (line.Is_comment())
                {
                    SF_out.insert(line);
                }
                SF.next();
            }  // while
            SF_out.write(SF_out_name); //write output sel file
        }
//input/output are single files

        else if (checkParameter(argc, argv, "-i") && checkParameter(argc, argv, "-o"))
        {
            raw22spi(fn_in, fn_out, raw_type, force_byte,
                header_size, Zdim, Ydim, Xdim, generate_inf, reverse_endian);
        }

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
        "\nPurpose: Convert from a 2d/3d raw images to Xmipp ones (and viceversa)"
        "\n        Input/Output can be either a single file or a set of them "
        "\n        (specified in a 'sel' file)"
        "\nParameter Values: (note space before value)"
        "\nI/O parameters"
        "\nESPECIFIC PARAMETERS FOR SINGLE-FILE CONVERSION"
        "\n    -i    file_in        input raw or Xmipp file"
        "\n    -o    file_out       output Xmipp or raw file"
        "\n   [-generate_inf]       only valid if the output is a raw image"
        "\nESPECIFIC PARAMETERS FOR SEL-FILE CONVERSION"
        "\n    -sel  input_file     input sel file"
        "\n    -oext input_file     extension for the output files if the input"
        "\n    files were specified in a sel file"
        "\nGENERAL PARAMETERS"
        "\n   [-f]           raw file is read/written in float format "
        "\n                         (byte by default)."
        "\n   [-16]                 raw file is read/written in 16-bit integers"
        "\n   [-8]                  raw file is read/written in 8-bit integers"
        "\n   [-reverse_endian]     by default, output has the same endiannes as input"
        "\n                         use this option to change endianness\n"
        "\n   [-header size=0]      Valid for raw to Xmipp conversions\n"
        "\n   [-s Zdim Ydim Xdim]   Z,Y,X dimensions for input files."
        "\n    For 2D raw images set the Zdim to 1"
        "\n   [-is_micrograph]      If this flag is provided the size is taken"
        "\n                         from <file_in>.inf\n"
        , argv[0]);
}

void raw22spi(const FileName &fn_in, const FileName &fn_out,
              char raw_type, bool force_byte,
              int header_size, int Zdim,  int Ydim, int Xdim,
              bool generate_inf, bool reverse_endian)
{

    VolumeXmipp    Vx;
    ImageXmipp     Ix;

    // Volume Xmipp --> Raw Volume
    if (Is_VolumeXmipp(fn_in))
    {
        Vx.read(fn_in);
        bool endianness = (reverse_endian) ? !Vx.reversed() : Vx.reversed();
        if (!force_byte) raw_type='f';
        switch (raw_type)
        {
        case 'b':
            ((Volume)Vx).write(fn_out, endianness, VBYTE);
            break;
        case 'h':
            ((Volume)Vx).write(fn_out, endianness, V16);
            break;
        case 'f':
            ((Volume)Vx).write(fn_out, endianness, VFLOAT);
            break;
        }
        // Image Xmipp --> Raw Image
    }
    else if (Is_ImageXmipp(fn_in))
    {
        Ix.read(fn_in);
        bool endianness = (reverse_endian) ? !Ix.reversed() : Ix.reversed();
        if (!force_byte) raw_type='f';
        int bits;
        switch (raw_type)
        {
        case 'b':
            ((Image)Ix).write(fn_out, endianness, IBYTE);
            bits = 8;
            break;
        case 'h':
            ((Image)Ix).write(fn_out, endianness, I16);
            bits = 16;
            break;
        case 'f':
            ((Image)Ix).write(fn_out, endianness, IFLOAT);
            bits = 32;
            break;
        }
        if (generate_inf)
        {
            std::ofstream fh((fn_out + ".inf").c_str());
            if (!fh)
                REPORT_ERROR(1, "Cannot create output .inf file");
            fh << "# Image width\n";
            fh << "Xdim=" << XSIZE(Ix()) << std::endl;
            fh << "# Image length\n";
            fh << "Ydim=" << YSIZE(Ix()) << std::endl;
            fh << "# Pixel depth\n";
            fh << "bitspersample=" << bits << std::endl;
            fh.close();
        }
    }
    // Raw image --> Xmipp Image
    else if (Zdim == 1)
    {
        Image *I = &Ix; // This is a trick for the compiler
        switch (raw_type)
        {
        case 'b':
            I->read(fn_in, 0, Ydim, Xdim, reverse_endian, IBYTE, header_size);
            break;
        case 'h':
            I->read(fn_in, 0, Ydim, Xdim, reverse_endian, I16, header_size);
            break;
        case 'f':
            I->read(fn_in, 0, Ydim, Xdim, reverse_endian, IFLOAT, header_size);
            break;
        }
        Ix.write(fn_out);
        // Raw Volume --> Xmipp Volume
    }
    else
    {
        Volume *V = &Vx; // This is a trick for the compiler
        switch (raw_type)
        {
        case 'b':
            V->read(fn_in, Zdim, Ydim, Xdim, reverse_endian, VBYTE, header_size);
            break;
        case 'h':
            V->read(fn_in, Zdim, Ydim, Xdim, reverse_endian, V16, header_size);
            break;
        case 'f':
            V->read(fn_in, Zdim, Ydim, Xdim, reverse_endian, VFLOAT, header_size);
            break;
        }
        Vx.write(fn_out);
    }
}
