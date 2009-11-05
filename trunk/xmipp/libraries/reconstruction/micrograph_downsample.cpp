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

#include "micrograph_downsample.h"
#include <data/args.h>
#include <data/mask.h>

// Read --------------------------------------------------------------------
void Prog_downsample_prm::read(int argc, char **argv, bool do_not_read_files)
{
    if (!do_not_read_files)
    {
        fn_micrograph  = getParameter(argc, argv, "-i");
        fn_downsampled = getParameter(argc, argv, "-o");
    }
    bitsMp         = textToInteger(getParameter(argc, argv, "-output_bits", "32"));
    if (bitsMp != 8 && bitsMp != 16 && bitsMp != 32)
        REPORT_ERROR(1, "Downsample: you must specify 8, 16 or 32 bits only");

    if (checkParameter(argc, argv, "-fourier"))
    {
        do_fourier=true;
        scale          = textToFloat(getParameter(argc, argv, "-fourier"));
    }
    else
    {
        do_fourier=false;
        Xstep          = textToInteger(getParameter(argc, argv, "-Xstep"));
        if (checkParameter(argc, argv, "-Ystep"))
            Ystep     = textToInteger(getParameter(argc, argv, "-Ystep"));
        else
            Ystep     = Xstep;
    }
    if (checkParameter(argc, argv, "-kernel"))
    {
        std::string aux = getParameter(argc, argv, "-kernel");
        int i = paremeterPosition(argc, argv, "-kernel");
        if (aux == "rectangle")
        {
            kernel_mode = KER_RECTANGLE;
            if (i + 2 >= argc)
                REPORT_ERROR(1, "Downsample: Not enough parameters after rectangle");
            Yrect = textToInteger(argv[i+2]);
            Xrect = textToInteger(argv[i+3]);
        }
        else if (aux == "circle")
        {
            kernel_mode = KER_CIRCLE;
            if (i + 2 >= argc)
                REPORT_ERROR(1, "Downsample: Not enough parameters after circle");
            r = textToFloat(argv[i+2]);
        }
        else if (aux == "gaussian")
        {
            kernel_mode = KER_GAUSSIAN;
            if (i + 3 >= argc)
                REPORT_ERROR(1, "Downsample: Not enough parameters after gaussian");
            if (Xstep != Ystep)
                REPORT_ERROR(1, "Downsample: You cannot apply different steps in this mode");
            r = textToFloat(argv[i+2]);
            sigma = textToFloat(argv[i+3]);
        }
        else if (aux == "pick")
        {
            kernel_mode = KER_PICK;
        }
        else if (aux == "sinc")
        {
            kernel_mode = KER_SINC;
            if (i + 3 >= argc)
                REPORT_ERROR(1, "Downsample: Not enough parameters after sinc");
            if (Xstep != Ystep)
                REPORT_ERROR(1, "Downsample: You cannot apply different steps in this mode");
            delta = textToFloat(argv[i+2]);
            Deltaw = textToFloat(argv[i+3]);
        }
        else
            REPORT_ERROR(1, "Downsample: Unknown kernel mode");
    }
    else
    {
        if(!do_fourier)
        {
            do_fourier=true;
            scale          = (double)1./Xstep;
        }
        if(checkParameter(argc, argv, "-old"))//undocument option
        {
            kernel_mode = KER_SINC;
            delta = 0.02;
            Deltaw = 1.0 / 10.0;
            do_fourier=false;
        }    
    }
    reversed = checkParameter(argc, argv, "-reverse_endian");
}

// Usage -------------------------------------------------------------------
void Prog_downsample_prm::usage() const
{
    std::cerr << "  [-output_bits <bits=32>]   : Must be 8, 16 or 32 bits\n"
              << "   -Xstep <xstep>            : Look at the documentation\n"
              << "  [-Ystep <ystep=xstep>]\n"
              << "  [-fourier <factor=0.3333>] : work in fourier space, factor < 1\n"
              << "  [-kernel rectangle <Ydim> <Xdim>]\n"
              << "  [-kernel circle    <r>]\n"
              << "  [-kernel gaussian  <r> <sigma>]\n"
              << "  [-kernel pick]\n"
              << "  [-kernel sinc <delta=0.02> <Deltaw=0.1>]\n"
              << "  [-reverse_endian]       : Reverse endian\n"
    ;
}
#ifdef NEVERDEFINED

// Produce command line ----------------------------------------------------
std::string Prog_downsample_prm::command_linezz() const
{
    std::string retval;
    retval += (std::string)"-i " + fn_micrograph + " ";
    retval += (std::string)"-o " + fn_downsampled + " ";
    retval += (std::string)"-output_bits " + integerToString(bitsMp) + " ";
    retval += (std::string)"-Xstep " + integerToString(Xstep) + " ";
    retval += (std::string)"-Ystep " + integerToString(Ystep) + " ";
    retval += (std::string)"-kernel ";
    switch (kernel_mode)
    {
    case KER_RECTANGLE:
        retval += (std::string)"rectangle" + integerToString(Yrect) + " " + integerToString(Xrect) + " ";
        break;
    case KER_CIRCLE:
        retval += (std::string)"circle" + floatToString(r) + " ";
        break;
    case KER_GAUSSIAN:
        retval += (std::string)"gaussian" + floatToString(r) + " " + floatToString(sigma) + " ";
        break;
    case KER_PICK:
        retval += (std::string)"pick";
        break;
    case KER_SINC:
        retval += (std::string)"sinc" + floatToString(delta) + " " + floatToString(Deltaw) + " ";
        break;
    }
    return retval;
}
#endif
// Generate kernel ---------------------------------------------------------
void Prog_downsample_prm::generate_kernel()
{
    // Integer Kernel
    Matrix2D<int>    ikernel;

    switch (kernel_mode)
    {
    case KER_RECTANGLE:
        kernel.resize(Yrect, Xrect);
        kernel.initConstant(1);
        break;
    case KER_CIRCLE:
        ikernel.resize(CEIL(2*r) + 1, CEIL(2*r) + 1);
        ikernel.setXmippOrigin();
        BinaryCircularMask(ikernel, r);
        typeCast(ikernel, kernel);
        break;
    case KER_GAUSSIAN:
        kernel.resize(CEIL(2*r) + 1, CEIL(2*r) + 1);
        kernel.setXmippOrigin();
        GaussianMask(kernel, sigma);
        break;
    case KER_PICK:
        kernel.resize(1, 1);
        kernel.initConstant(1);
        break;
    case KER_SINC:
        SeparableSincKaiserMask(kernel, (1.0 / Xstep),
                                delta, Deltaw);
        break;
    }
    kernel.setXmippOrigin();
    // Keep energy constant
    // kernel /=sqrt(kernel.sum2());
    // Keep average value constant
    kernel /= kernel.sum();
    //ImageXmipp save;
    //save() = kernel;
    //save.write("PPPkernel.xmp");
}

// Create output inf file --------------------------------------------------
void Prog_downsample_prm::create_empty_output_file()
{
    std::cerr << "Creating empty downsampled file ...\n";
    if (do_fourier) 
    {
        Ypdim = FLOOR((double)Ydim *scale);
        Xpdim = FLOOR((double)Xdim *scale);
    }
    else
    {
        Ypdim = FLOOR(Ydim / Ystep);
        Xpdim = FLOOR(Xdim / Xstep);
    }
    create_empty_file(fn_downsampled, ((unsigned long long)Ypdim)*
      	 Xpdim*bitsMp / 8);

    std::ofstream fh_downsample_inf;
    fh_downsample_inf.open((fn_downsampled + ".inf").c_str());
    if (!fh_downsample_inf)
        REPORT_ERROR(1, (std::string)"Downsample: Cannot open " + fn_downsampled +
                     ".inf");
    fh_downsample_inf << "# Generated by Downsample\n";
    fh_downsample_inf << "# Original file: " << fn_micrograph << std::endl;
    if(!do_fourier)
        fh_downsample_inf << "# Ystep x Xstep: " << Ystep << " x "
                          << Xstep << std::endl;
    if(do_fourier)
        fh_downsample_inf << "# scale : " << scale << std::endl;
    fh_downsample_inf << "# Image width\n";
    fh_downsample_inf << "Xdim= " << Xpdim << std::endl;
    fh_downsample_inf << "# Image length\n";
    fh_downsample_inf << "Ydim= " << Ypdim << std::endl;
    fh_downsample_inf << "# Pixel depth\n";
    fh_downsample_inf << "bitspersample= " << bitsMp << std::endl;
    fh_downsample_inf.close();
}

// Open input micrograph ---------------------------------------------------
void Prog_downsample_prm::open_input_micrograph()
{
    M.open_micrograph(fn_micrograph, reversed);
    bitsM = M.depth();
    M.size(Xdim, Ydim);
}

// Close input micrograph --------------------------------------------------
void Prog_downsample_prm::close_input_micrograph()
{
    M.close_micrograph();
}

// Downsample micrograph ---------------------------------------------------
void Prog_downsample_prm::Downsample() const
{
    Micrograph Mp;
    Mp.open_micrograph(fn_downsampled, reversed);
    downsample(M, Xstep, Ystep, kernel, Mp,do_fourier);
    Mp.close_micrograph();
}
