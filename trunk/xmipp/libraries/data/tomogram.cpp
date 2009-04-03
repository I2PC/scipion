/***************************************************************************
 *
 * Authors: Carlos Oscar (coss@cnb.uam.es)
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

#include "tomogram.h"
#include "args.h"
#include "volume.h"

#include <fstream>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#ifdef LINUX
#include <unistd.h>
#endif


/* Clear ------------------------------------------------------------------- */
void Tomogram::clear()
{
    fn_tomogram = "";
    fh_tomogram = -1;
    Zdim = Xdim = Ydim = -1;
    __depth = -1;
    m8 = NULL;
    m16 = NULL;
    um16 = NULL;
    m32 = NULL;
}

/* Open tomogram --------------------------------------------------------- */
void Tomogram::open_tomogram(const FileName &_fn_tomogram,
                             bool reversed)
{
    struct stat info;

    // tomogram name
    fn_tomogram = _fn_tomogram;

    // Look for tomogram dimensions
    // check if the file format is spider
    if (Is_VolumeXmipp(fn_tomogram))
    {
        headerXmipp     header;
        header.read(fn_tomogram);
        Xdim = header.iXdim();
        Ydim = header.iYdim();
        Zdim = header.iZdim();
        __offset = header.get_header_size();
        __depth = 32;
        reversed = header.reversed();
    }
    else
    {
        fn_inf = fn_tomogram.add_extension("inf");
        FILE *fh_inf = fopen(fn_inf.c_str(), "r");
        if (!fh_inf)
            REPORT_ERROR(1, (std::string)"Tomogram::open_tomogram: Cannot find " +
                         fn_inf);
        Xdim = textToInteger(getParameter(fh_inf, "Xdim"));
        Ydim = textToInteger(getParameter(fh_inf, "Ydim"));
        Zdim = textToInteger(getParameter(fh_inf, "Zdim"));
        __depth = textToInteger(getParameter(fh_inf, "bitspersample"));
        if (checkParameter(fh_inf, "offset"))
            __offset = textToInteger(getParameter(fh_inf, "offset"));
        else
            __offset = 0;
        if (checkParameter(fh_inf, "is_signed"))
            __is_signed = (getParameter(fh_inf, "is_signed") == "true" ||
                           getParameter(fh_inf, "is_signed") == "TRUE");
        else
            __is_signed = false;
        fclose(fh_inf);
    }
    XYdim = Xdim * Ydim;

    // Open tomogram and map
    fh_tomogram = open(fn_tomogram.c_str(), O_RDWR, S_IREAD | S_IWRITE);
    if (fh_tomogram == -1)
        REPORT_ERROR(1, (std::string)"Tomogram::open_tomogram: There is a "
                     "problem opening " + fn_tomogram);
    char *aux_ptr;
    switch (__depth)
    {
    case 8:
        m8 = (unsigned char *) mmap(0, (__depth / 8) * Zdim * Ydim * Xdim + __offset,
                                    PROT_READ | PROT_WRITE, MAP_SHARED, fh_tomogram, 0);
        if (m8 == MAP_FAILED)
            REPORT_ERROR(1, (std::string)"Tomogram::open_tomogram: cannot map " +
                         _fn_tomogram + " in memory");
        m8 += __offset;
        break;
    case 16:
        if (__is_signed)
        {
            m16 = (short int *) mmap(0, (__depth / 8) * Zdim * Ydim * Xdim + __offset,
                                     PROT_READ | PROT_WRITE, MAP_SHARED, fh_tomogram, 0);
            if (m16 == MAP_FAILED)
            {
                /*
                switch (errno) {
                case EACCES:    std::cout << "EACCES:   \n"; break;
                  case EAGAIN:    std::cout << "EAGAIN:   \n"; break;
                case EBADF:     std::cout << "EBADF:    \n"; break;
                case EINVAL:    std::cout << "EINVAL:   \n"; break;
                case EMFILE:    std::cout << "EMFILE:   \n"; break;
                case ENODEV:    std::cout << "ENODEV:   \n"; break;
                case ENOMEM:    std::cout << "ENOMEM:   \n"; break;
                case ENOTSUP:   std::cout << "ENOTSUP:  \n"; break;
                case ENXIO:     std::cout << "ENXIO:    \n"; break;
                case EOVERFLOW: std::cout << "EOVERFLOW:\n"; break;
                }
                */
                REPORT_ERROR(1, (std::string)"Tomogram::open_tomogram: cannot map " +
                             _fn_tomogram + " in memory");
            }
            aux_ptr = (char *)m16;

            aux_ptr += __offset;
            m16 = (short int *) aux_ptr;
        }
        else
        {
            um16 = (unsigned short int *) mmap(0, (__depth / 8) * Zdim * Ydim * Xdim + __offset,
                                               PROT_READ | PROT_WRITE, MAP_SHARED, fh_tomogram, 0);
            if (um16 == MAP_FAILED)
            {
                /*
                switch (errno) {
                case EACCES:    std::cout << "EACCES:   \n"; break;
                  case EAGAIN:    std::cout << "EAGAIN:   \n"; break;
                case EBADF:     std::cout << "EBADF:    \n"; break;
                case EINVAL:    std::cout << "EINVAL:   \n"; break;
                case EMFILE:    std::cout << "EMFILE:   \n"; break;
                case ENODEV:    std::cout << "ENODEV:   \n"; break;
                case ENOMEM:    std::cout << "ENOMEM:   \n"; break;
                case ENOTSUP:   std::cout << "ENOTSUP:  \n"; break;
                case ENXIO:     std::cout << "ENXIO:    \n"; break;
                case EOVERFLOW: std::cout << "EOVERFLOW:\n"; break;
                }
                */
                REPORT_ERROR(1, (std::string)"Tomogram::open_tomogram: cannot map " +
                             _fn_tomogram + " in memory");
            }
            aux_ptr = (char *)um16;
            aux_ptr += __offset;
            um16 = (unsigned short int *) aux_ptr;
        }
        break;
    case 32:
        // Map file in memory
        m32 = (float*) mmap(0, (__depth / 8) * Zdim * Ydim * Xdim + __offset,
                            PROT_READ | PROT_WRITE, MAP_SHARED, fh_tomogram, 0);
        if (m32 == MAP_FAILED)
            REPORT_ERROR(1, (std::string)"Tomogram::open_tomogram: cannot map " +
                         _fn_tomogram + " in memory");
        aux_ptr = (char *)m32;
        aux_ptr += __offset;
        m32 = (float *) aux_ptr;
        break;
    default:
        REPORT_ERROR(1, "Tomogram::open_tomogram: depth is not 8, 16 nor 32");
    }
    __reversed = reversed;
}

/* Close tomogram -------------------------------------------------------- */
void Tomogram::close_tomogram()
{
    if (fh_tomogram != -1)
    {
        close(fh_tomogram);
        if (__depth == 8)
        {
            m8 -= __offset;
            munmap((char *)m8, Zdim*Ydim*Xdim*__depth / 8 + __offset);
        }
        else if (__depth == 16)
        {
            if (__is_signed)
            {
                char *aux_ptr = (char *)m16;
                aux_ptr -= __offset;
                m16 = (short int *)aux_ptr;
                munmap((char *)m16, Zdim*Ydim*Xdim*__depth / 8 + __offset);
            }
            else
            {
                char *aux_ptr = (char *)um16;
                aux_ptr -= __offset;
                um16 = (unsigned short int *)aux_ptr;
                munmap((char *)um16, Zdim*Ydim*Xdim*__depth / 8 + __offset);
            }
        }
        else if (__depth == 32)
        {
            char *aux_ptr = (char *)m32;
            aux_ptr -= __offset;
            m32 = (float *)aux_ptr;
            munmap((char *)m32, Zdim*Ydim*Xdim*__depth / 8 + __offset);
        }
    }
}

/* Get piece --------------------------------------------------------------- */
void Tomogram::get_piece(Matrix1D<int> &r0, Matrix1D<int> &length,
                         Matrix3D<double> &piece)
{
    Matrix1D<int> rF = r0 + length - 1;
    std::cout << r0.transpose() << std::endl;
    std::cout << length.transpose() << std::endl;
    std::cout << Xdim << " " << Ydim << " " << Zdim << std::endl;
    if (XX(rF) >= Xdim)
        XX(r0) = Xdim - XX(length) - 1;
    if (YY(rF) >= Ydim)
        YY(r0) = Ydim - YY(length) - 1;
    if (ZZ(rF) >= Zdim)
        ZZ(r0) = Zdim - ZZ(length) - 1;
    std::cout << r0.transpose() << std::endl;
    if (XX(r0) < 0 || YY(r0) < 0 || ZZ(r0) < 0)
        REPORT_ERROR(1, "Tomogram::get_piece: piece does not fit into tomogram");

    piece.resize(ZZ(length), YY(length), XX(length));
    int k, i, j, kp, ip, jp;
    for (kp = ZZ(r0), k = 0; k < ZSIZE(piece); k++, kp++)
        for (ip = YY(r0), i = 0; i < YSIZE(piece); i++, ip++)
            for (jp = XX(r0), j = 0; j < XSIZE(piece); j++, jp++)
                DIRECT_VOL_ELEM(piece, k, i, j) = (*this)(jp, ip, kp);
}

/* Set piece --------------------------------------------------------------- */
void Tomogram::set_piece(Matrix1D<int> &r0, Matrix1D<int> &length,
                         Matrix3D<double> &piece)
{
    int k, i, j, kp, ip, jp;
    for (kp = ZZ(r0), k = 0; k < ZSIZE(piece); k++, kp++)
        for (ip = YY(r0), i = 0; i < YSIZE(piece); i++, ip++)
            for (jp = XX(r0), j = 0; j < XSIZE(piece); j++, jp++)
                set_val(jp, ip, kp, DIRECT_VOL_ELEM(piece, k, i, j));
}
