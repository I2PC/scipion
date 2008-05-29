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

#include "micrograph.h"
#include "args.h"
#include "selfile.h"
#include "mask.h"
#include "geometry.h"

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
void Micrograph::clear()
{
    single_particle.clear();
    coords.clear();
    fn_coords = fn_micrograph = "";
    X_window_size = Y_window_size = -1;
    fh_micrograph = -1;
    Xdim = Ydim = -1;
    __depth = -1;
    m8 = NULL;
    m16 = NULL;
    um16 = NULL;
    m32 = NULL;
    compute_transmitance = false;
    compute_inverse = false;
    __scaling_valid = false;
    /* __in_core=FALSE;*/
}

/* Open micrograph --------------------------------------------------------- */
void Micrograph::open_micrograph(const FileName &_fn_micrograph,
                                 /*bool in_core,*/ bool reversed)
{
    struct stat info;

    // Micrograph name
    fn_micrograph = _fn_micrograph;

    // Look for micrograph dimensions
    // check if the file format is spider
    if (Is_ImageXmipp(fn_micrograph))
    {
        headerXmipp     header;
        header.read(fn_micrograph);
        float fXdim, fYdim;
        header.getDimensionension(fYdim, fXdim);
        Xdim = (int) fXdim;
        Ydim = (int)fYdim;
        __offset = header.get_header_size();
        __depth = 32;
        reversed = header.reversed();
    }
    else
    {
        fn_inf = fn_micrograph.add_extension("inf");
        FILE *fh_inf = fopen(fn_inf.c_str(), "r");
        if (!fh_inf)
            REPORT_ERROR(1, (std::string)"Micrograph::open_micrograph: Cannot find " +
                         fn_inf);
        Xdim = textToInteger(getParameter(fh_inf, "Xdim"));
        Ydim = textToInteger(getParameter(fh_inf, "Ydim"));
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
    // Open micrograph and map
    fh_micrograph = open(fn_micrograph.c_str(), O_RDWR, S_IREAD | S_IWRITE);
    if (fh_micrograph == -1)
        REPORT_ERROR(1, (std::string)"Micrograph::open_micrograph: There is a "
                     "problem opening " + fn_micrograph +
                     "\nCheck that the file has write permission");
    char *aux_ptr;
    switch (__depth)
    {
    case 8:
        /* if (!in_core) { */
        m8 = (unsigned char *) mmap(0, (__depth / 8) * Ydim * Xdim + __offset,
                                    PROT_READ | PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
        if (m8 == MAP_FAILED)
            REPORT_ERROR(1, (std::string)"Micrograph::open_micrograph: cannot map " +
                         _fn_micrograph + " in memory");
        m8 += __offset;
        /* } else {
           m8=new unsigned char (Ydim*Xdim*__depth/8);
           int length=Ydim*Xdim*__depth/8;
           int read_length=read(fh_micrograph,m8,length);
           std::cout << SSIZE_MAX << std::endl;
           std::cout << read_length << std::endl;
           if (read_length!=length)
              REPORT_ERROR(1,(std::string)"Micrograph::open_micrograph: cannot read "+
                 _fn_micrograph+" in memory");
        } */
        break;
    case 16:
        /* if (!in_core) { */
        if (__is_signed)
        {
            m16 = (short int *) mmap(0, (__depth / 8) * Ydim * Xdim + __offset,
                                     PROT_READ | PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
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
                REPORT_ERROR(1, (std::string)"Micrograph::open_micrograph: cannot map " +
                             _fn_micrograph + " in memory");
            }
            aux_ptr = (char *)m16;

            aux_ptr += __offset;
            m16 = (short int *) aux_ptr;
        }
        else
        {
            um16 = (unsigned short int *) mmap(0, (__depth / 8) * Ydim * Xdim + __offset,
                                               PROT_READ | PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
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
                REPORT_ERROR(1, (std::string)"Micrograph::open_micrograph: cannot map " +
                             _fn_micrograph + " in memory");
            }
            aux_ptr = (char *)um16;
            aux_ptr += __offset;
            um16 = (unsigned short int *) aux_ptr;
        }
        /* } else {
           m16=new short int (Ydim*Xdim*__depth/8);
           int length=Ydim*Xdim*__depth/8;
           if (read(fh_micrograph,m16,length)!=length)
              REPORT_ERROR(1,(std::string)"Micrograph::open_micrograph: cannot read "+
                 _fn_micrograph+" in memory");
        } */
        break;
    case 32:
        // Map file in memory
        m32 = (float*) mmap(0, (__depth / 8) * Ydim * Xdim + __offset,
                            PROT_READ | PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
        if (m32 == MAP_FAILED)
            REPORT_ERROR(1, (std::string)"Micrograph::open_micrograph: cannot map " +
                         _fn_micrograph + " in memory");
        aux_ptr = (char *)m32;
        aux_ptr += __offset;
        m32 = (float *) aux_ptr;
        break;
    default:
        REPORT_ERROR(1, "Micrograph::open_micrograph: depth is not 8, 16 nor 32");
    }
    /*__in_core=in_core; */
    __reversed = reversed;
}

/* Close micrograph -------------------------------------------------------- */
void Micrograph::close_micrograph()
{
    if (fh_micrograph != -1)
    {
        close(fh_micrograph);
        /* if (!__in_core) { */
        if (__depth == 8)
        {
            m8 -= __offset;
            munmap((char *)m8, Ydim*Xdim*__depth / 8 + __offset);
        }
        else if (__depth == 16)
        {
            if (__is_signed)
            {
                char *aux_ptr = (char *)m16;
                aux_ptr -= __offset;
                m16 = (short int *)aux_ptr;
                munmap((char *)m16, Ydim*Xdim*__depth / 8 + __offset);
            }
            else
            {
                char *aux_ptr = (char *)um16;
                aux_ptr -= __offset;
                um16 = (unsigned short int *)aux_ptr;
                munmap((char *)um16, Ydim*Xdim*__depth / 8 + __offset);
            }
        }
        else if (__depth == 32)
        {
            char *aux_ptr = (char *)m32;
            aux_ptr -= __offset;
            m32 = (float *)aux_ptr;
            munmap((char *)m32, Ydim*Xdim*__depth / 8 + __offset);
        }
        /* } else {
           if      (__depth== 8) delete m8;
           else if (__depth==16) delete m16;
        } */
    }
}

/* Compute 8 bit scaling --------------------------------------------------- */
void Micrograph::compute_8_bit_scaling()
{
    std::cerr << "Computing 8 bit scaling ...\n";

    // Compute minimum and maximum value
    float minval, maxval;
    minval = maxval = (*this)(0, 0);
    for (int i = 0; i < Ydim; i++)
    {
        for (int j = 0; j < Xdim; j++)
        {
            float tmp = (*this)(j, i);
            if (tmp < minval)
                minval = tmp;
            else if (tmp > maxval)
                maxval = tmp;
            /*
            if(maxval > 32000)
              std::cout << "(i,j) max min valuefloat value" << i << " " << j
                 << " " << maxval << " " << minval << " "<< tmp
                 << " " << (*this)(j,i) << std::endl;
            */
        }
    }
    // Compute output range
    float minF, maxF;
    if (minval < 0)
    {
        minF = 0;
        maxF = XMIPP_MIN(255, maxval - minval);
    }
    else if (maxval > 255)
    {
        minF = XMIPP_MAX(0, minval - (maxval - 255));
        maxF = 255;
    }
    else if (maxval - minval < 32)
    {
        minF = 0;
        maxF = 255;
    }
    else
    {
        minF = minval;
        maxF = maxval;
    }

    // Compute scaling
    __a = (maxF - minF) / (maxval - minval);
    __b = minF - __a * minval;
    __scaling_valid = true;
    //std::cerr <<  "__a  " << __a  << "__b" << __b << std::endl;

}

/* Write as 8 bits --------------------------------------------------------- */
void Micrograph::write_as_8_bits(const FileName &fn8bits)
{
    if (!__scaling_valid) compute_8_bit_scaling();

    // Create empty output file
    create_empty_file(fn8bits, ((unsigned long long)Ydim)*Xdim);

    std::ofstream fh8bits_inf;
    fh8bits_inf.open((fn8bits + ".inf").c_str());
    if (!fh8bits_inf)
        REPORT_ERROR(1, (std::string)"write_as_8_bits: Cannot open " + fn8bits +
                     ".inf");
    fh8bits_inf << "# Generated by write_as_8_bits\n";
    fh8bits_inf << "# Original file: " << fn_micrograph << std::endl;
    fh8bits_inf << "# Image width\n";
    fh8bits_inf << "Xdim= " << Xdim << std::endl;
    fh8bits_inf << "# Image length\n";
    fh8bits_inf << "Ydim= " << Ydim << std::endl;
    fh8bits_inf << "# Pixel depth\n";
    fh8bits_inf << "bitspersample=8\n";
    fh8bits_inf.close();

    // Open micrograph
    Micrograph Mp;
    Mp.open_micrograph(fn8bits, false);
    for (int y=0; y<Ydim; y++)
        for (int x=0; x<Xdim; x++)
           Mp.set_val(x,y,val8(x,y));
    Mp.close_micrograph();
}

/* Save coordinates to disk ------------------------------------------------ */
void Micrograph::write_coordinates(int label, const FileName &_fn_coords)
{
    std::ofstream fh;
    if (_fn_coords != "")
        fn_coords = _fn_coords;
    fh.open(fn_coords.c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR(1, (std::string)"Micrograph::write: File " + fn_coords +
                     " cannot be openned for output");
    int imax = coords.size();
    fh << "# <X position> <Y position>\n";
    for (int i = 0; i < imax; i++)
        if (coords[i].valid && coords[i].label == label)
            fh << coords[i].X << " " << coords[i].Y << std::endl;
    fh.close();
}

/* Read coordinates from disk ---------------------------------------------- */
void Micrograph::read_coordinates(int label, const FileName &_fn_coords)
{
    std::ifstream  fh;
    int            line_no = 0;
    std::string    line;

    fn_coords = _fn_coords;
    fh.open(fn_coords.c_str(), std::ios::in);
    if (!fh)
        REPORT_ERROR(1, (std::string)"Micrograph::read: File " + fn_coords + " not found");

    // Count the number of lines
    fh.peek();
    while (!fh.eof())
    {
        getline(fh, line);
        if (line.length() > 0 && line[0] != '#' && line[0] != ';')
            line_no++;
        fh.peek();
    }
    fh.close();
    fh.clear();

    // Resize coordinate list and read
    fh.open(fn_coords.c_str(), std::ios::in);
    coords.reserve(line_no);
    struct Particle_coords aux;
    aux.valid = true;
    aux.label = label;
    fh.peek();
    while (!fh.eof())
    {
        getline(fh, line);
        if (line.length() > 0 && line[0] != '#' && line[0] != ';')
        {
            int converted_elements = sscanf(line.c_str(), "%d %d",
                                            &aux.X, &aux.Y);
            if (converted_elements != 2)
                std::cerr << "Ignoring line: " << line << std::endl;
            else
                coords.push_back(aux);
        }
        fh.peek();
    }
    fh.close();
}

/* Transform all coordinates ---------------------------------------------- */
void Micrograph::transform_coordinates(const Matrix2D<double> &M)
{
    Matrix1D<double> m(3);
    SPEED_UP_temps;
    
    int imax = coords.size();
    for (int i = 0; i < imax; i++)
    {
	if (coords[i].valid)
	{
	    VECTOR_R3(m,coords[i].X, coords[i].Y,1);
	    M3x3_BY_V3x1(m, M, m);
	    coords[i].X=(int)XX(m);
	    coords[i].Y=(int)YY(m);
	}
    }
}

/* Multiply coordinates by a constant -------------------------------------- */
void Micrograph::scale_coordinates(const double &c)
{
    int imax = coords.size();
    for (int i = 0; i < imax; i++)
    {
	if (coords[i].valid)
	{
	    coords[i].X = (int)coords[i].X*c;
	    coords[i].Y = (int)coords[i].Y*c;
	}
    }

}

/* Scissor ----------------------------------------------------------------- */
int Micrograph::scissor(const Particle_coords &P, ImageT<double> &result,
                        double Dmin, double Dmax, double scaleX, double scaleY,
                        bool only_check)
{
    if (X_window_size == -1 || Y_window_size == -1)
        REPORT_ERROR(1, "Micrograph::scissor: window size not set");

    result().resize(Y_window_size, X_window_size);
    int i0 = ROUND(scaleY * P.Y) + FIRST_XMIPP_INDEX(Y_window_size);
    int iF = ROUND(scaleY * P.Y) + LAST_XMIPP_INDEX(Y_window_size);
    int j0 = ROUND(scaleX * P.X) + FIRST_XMIPP_INDEX(X_window_size);
    int jF = ROUND(scaleX * P.X) + LAST_XMIPP_INDEX(X_window_size);
    int retval = 1;
    double range, temp;
    range = Dmax - Dmin;
    if (i0 < 0 || iF >= Ydim || j0 < 0 || jF >= Xdim)
    {
        result().initZeros();
        retval = 0;
    }
    else
        if (!only_check)
        {
            for (int i = i0; i <= iF; i++)
                for (int j = j0; j <= jF; j++)
                {
                    if (compute_transmitance)
                    {
                        if ((*this)(j, i) < 1)
                            temp = (*this)(j, i);
                        else
                            temp = log10((double)(*this)(j, i));
                        if (compute_inverse)
                            result(i - i0, j - j0) = (Dmax - temp) / range;
                        else
                            result(i - i0, j - j0) = (temp - Dmin) / range;
                    }
                    else
                    {
                        if (compute_inverse)
                            result(i - i0, j - j0) = (Dmax - (*this)(j, i)) / range;
                        else
                            result(i - i0, j - j0) = (*this)(j, i);
                    }
                }
        }
    return retval;
}

/* Produce all images ------------------------------------------------------ */
void Micrograph::produce_all_images(int label, const FileName &fn_root,
                                    int starting_index, const FileName &fn_image, double ang, double tilt,
                                    double psi)
{
    SelFile SF;
    FileName fn_out;
    ImageXmipp I;
    Micrograph *M;

    // Set Source image
    if (fn_image == "")
        M = this;
    else
    {
        M = new Micrograph;
        M->open_micrograph(fn_image, __reversed);
        M->set_window_size(X_window_size, Y_window_size);
        M->set_transmitance_flag(compute_transmitance);
        M->set_inverse_flag(compute_inverse);
    }

    // Set scale for particles
    int MXdim, MYdim, thisXdim, thisYdim;
    M->size(MXdim, MYdim);
    this->size(thisXdim, thisYdim);
    double scaleX = (double)MXdim / thisXdim;
    double scaleY = (double)MYdim / thisYdim;

    // Compute max and minimum if compute_transmitance
    // or compute_inverse flags are ON
    double Dmax, Dmin;
    if (compute_transmitance || compute_inverse)
    {
        (*this).computeDoubleMinMax(Dmin, Dmax);
        //#define DEBUG66
#ifdef DEBUG66

        std::cout << "Min= " << Dmin << " Dmax" << Dmax << std::endl;
#endif
#undef DEBUG66

        if (compute_transmitance)
        {
            if (Dmin > 1)
                Dmin = log10(Dmin);
            if (Dmax > 1)
                Dmax = log10(Dmax);
        }
    }
    // Scissor all particles
    if (ang != 0)
        std::cout << "Angle from Y axis to tilt axis " << ang << std::endl
                  << "   applying apropriate rotation\n";
    int i = starting_index;
    int nmax = ParticleNo();
    for (int n = 0; n < nmax; n++)
        if (coords[n].valid && coords[n].label == label)
        {
            fn_out.compose(fn_root, i++, "xmp");
            if (!M->scissor(coords[n], (Image &) I, Dmin, Dmax, scaleX, scaleY))
            {
                std::cout << "Particle " << fn_out << " is very near the border, "
                          << "corresponding image is set to blank\n";
                SF.insert(fn_out, SelLine::DISCARDED);
            }
            else
                SF.insert(fn_out);
            //  if (ang!=0) I().rotate(-ang);
            I.set_rot((float)ang);
            I.set_tilt((float)tilt);
            I.set_psi((float)psi);
            I.write(fn_out);
        }
    if (labels[label] != "")
    {
        SF.write(fn_micrograph.remove_directories() + "." + labels[label] + ".sel");
        write_coordinates(label, fn_micrograph + "." + labels[label] + ".pos");
    }
    else
    {
        SF.write(fn_micrograph.remove_directories() + ".sel");
        write_coordinates(label, fn_micrograph + ".pos");
    }

    // Free source image??
    if (fn_image != "")
    {
        M->close_micrograph();
        delete M;
    }
}

/* Search coordinate near a position --------------------------------------- */
int Micrograph::search_coord_near(int x, int y, int prec) const
{
    int imax = coords.size();
    int prec2 = prec * prec;
    for (int i = 0; i < imax; i++)
        if ((coords[i].X - x)*(coords[i].X - x) + (coords[i].Y - y)*(coords[i].Y - y) < prec2
            && coords[i].valid)
            return i;
    return -1;
}

/* Invalidate a coordinate ------------------------------------------------- */
void Micrograph::invalidate_coord(int n)
{
    if (n < 0 || n >= ParticleNo())
        REPORT_ERROR(1, "Micrograph::invalidate_coord: Index out of range");
    coords[n].valid = false;
}

/* Add coordinate ---------------------------------------------------------- */
int Micrograph::add_coord(int x, int y, int label)
{
    struct Particle_coords aux;
    aux.valid = true;
    aux.X = x;
    aux.Y = y;
    aux.label = label;
    coords.push_back(aux);
    return coords.size() - 1;
}

/* Move last coordinate ---------------------------------------------------- */
void Micrograph::move_last_coord_to(int x, int y)
{
    if (coords.size() > 0)
    {
        coords.back().X = x;
        coords.back().Y = y;
    }
}

/* Downsample -------------------------------------------------------------- */
void downsample(const Micrograph &M, int Xstep, int Ystep,
                const Matrix2D<double> &kernel, Micrograph &Mp)
{
    // Find first and last indexes in each direction
    // it is assumed that (0,0) is within the kernel
    int Ydim, Xdim, Ypdim, Xpdim;
    M.size(Xdim, Ydim);
    Mp.size(Xpdim, Ypdim);
    int x0 = 0;
    int y0 = 0;
    int xF = Xdim;
    int yF = Ydim;

    double pixval;
    int ii, y, jj, x;
    time_config();

    // Look for input/output ranges
    double a = 1;
    double b = 0;
    double scale = 1;

    if (Mp.depth() != 32)
    {
        double imin, imax;
        double omin, omax;
        bool ifirst = true, ofirst = true;

        if (M.depth() != 32)
            scale = (pow(2.0, Mp.depth()) - 1.0) / (pow(2.0, M.depth()) - 1.0);
        else if (M.depth() == 32)
            scale = 1;

        init_progress_bar(yF / Ystep);
        for (ii = 0, y = y0; y < yF; y += Ystep, ii++)
        {
            for (jj = 0, x = x0; x < xF; x += Xstep, jj++)
            {
                pixval = 0;
                FOR_ALL_ELEMENTS_IN_MATRIX2D(kernel)
                {
                    int j2 = intWRAP(j + x, 0, xF - 1);
                    int i2 = intWRAP(i + y, 0, yF - 1);
                    if (ifirst)
                    {
                        imin = imax = M(j2, i2);
                        ifirst = false;
                    }
                    else
                    {
                        imin = XMIPP_MIN(imin, M(j2, i2));
                        imax = XMIPP_MAX(imax, M(j2, i2));
                    }
                    pixval += kernel(i, j) * M(j2, i2);
                }
                pixval *= scale;
                if (ii < Ypdim && jj < Xpdim)
                {
                    if (ofirst)
                    {
                        omin = omax = pixval;
                        ofirst = false;
                    }
                    else
                    {
                        omin = XMIPP_MIN(omin, pixval);
                        omax = XMIPP_MAX(omax, pixval);
                    }
                }
            }
            if (ii % 50 == 0)
                progress_bar(ii);
        }
        progress_bar(yF / Ystep);

        // Compute range transformation
        double irange = imax - imin;
        double orange = omax - omin;

        if (M.depth() != 32)
        {
            a = scale * irange / orange;
            b = -omin;
        }
        else
        {
            a = (pow(2.0, Mp.depth()) - 1.0) / orange;
            scale = 1;
            b = -omin;
        }
    }

    // Really downsample
    init_progress_bar(yF / Ystep);
    for (ii = 0, y = y0; y < yF; y += Ystep, ii++)
    {
        for (jj = 0, x = x0; x < xF; x += Xstep, jj++)
        {
            pixval = 0;
            for (int i=STARTINGY(kernel); i<=FINISHINGY(kernel); i++)
            {
                int i2 = intWRAP(i + y, 0, yF - 1);
                for (int j=STARTINGX(kernel); j<=FINISHINGX(kernel); j++)
                {
                    int j2 = intWRAP(j + x, 0, xF - 1);
                    pixval += kernel(i, j) * M(j2, i2);
                }
            }

            if (ii < Ypdim && jj < Xpdim)
                if (Mp.depth() != 32)
                    Mp.set_val(jj, ii, FLOOR(a*(pixval*scale + b)));
                else
                    Mp.set_val(jj, ii, pixval);
        }
        if (ii % 50 == 0)
            progress_bar(ii);
    }
    progress_bar(yF / Ystep);
}

/* Normalizations ---------------------------------------------------------- */
void normalize_OldXmipp(Image *I)
{
    double avg, stddev, min, max;
    (*I)().computeStats(avg, stddev, min, max);
    (*I)() -= avg;
    (*I)() /= stddev;
}

void normalize_Near_OldXmipp(Image *I, const Matrix2D<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    (*I)().computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
                                     stddevbg);
    (*I)() -= avg;
    (*I)() /= stddevbg;
}

void normalize_OldXmipp_decomposition(Image *I, const Matrix2D<int> &bg_mask,
                                     const Matrix2D<double> *mask)
{
    double avgbg, stddevbg, minbg, maxbg;
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
                                     stddevbg);
    (*I)() -= avgbg;
    (*I)() /= stddevbg;
    if (mask != NULL)
        (*I)() *= *mask;
    normalize_OldXmipp(I);
}

void normalize_Michael(Image *I, const Matrix2D<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    (*I)().computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
                                     stddevbg);
    if (avgbg > 0)
    {
        (*I)() -= avgbg;
        (*I)() /= avgbg;
    }
    else
    { // To avoid the contrast inversion
        (*I)() -= (avgbg - min);
        (*I)() /= (avgbg - min);
    }
}

void normalize_NewXmipp(Image *I, const Matrix2D<int> &bg_mask)
{
    double avgbg, stddevbg, minbg, maxbg;
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
                                     stddevbg);
    (*I)() -= avgbg;
    (*I)() /= stddevbg;
}

void normalize_NewXmipp2(Image *I, const Matrix2D<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    (*I)().computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
                                     stddevbg);
    (*I)() -= avgbg;
    (*I)() /= ABS(avg - avgbg);
}

void normalize_ramp(Image *I, const Matrix2D<int> &bg_mask)
{
    fit_point          onepoint;
    std::vector<fit_point>  allpoints;
    double             pA, pB, pC;
    double             avgbg, stddevbg, minbg, maxbg;

    // Fit a least squares plane through the background pixels
    allpoints.clear();
    (*I)().setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)())
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
            onepoint.x = j;
            onepoint.y = i;
            onepoint.z = MAT_ELEM((*I)(), i, j);
            onepoint.w = 1.;
            allpoints.push_back(onepoint);
        }
    }
    least_squares_plane_fit(allpoints, pA, pB, pC);
    // Substract the plane from the image
    FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)())
    {
        MAT_ELEM((*I)(), i, j) -= pA * j + pB * i + pC;
    }
    // Divide by the remaining std.dev. in the background region
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
                                     stddevbg);
    (*I)() /= stddevbg;

}

void normalize_remove_neighbours(Image *I, 
				 const Matrix2D<int> &bg_mask,
                                 const double &threshold)
{
    fit_point          onepoint;
    std::vector<fit_point>  allpoints;
    double             pA, pB, pC;
    double             avgbg, stddevbg, minbg, maxbg, aux, newstddev;
    double             sum1 = 0.;
    double             sum2 = 0;
    int                N = 0;

    // Fit a least squares plane through the background pixels
    allpoints.clear();
    (*I)().setXmippOrigin();
    
    // Get initial statistics
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,stddevbg);

    // Fit plane through those pixels within +/- threshold*sigma
    FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)())
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
	    if ( ABS(avgbg - MAT_ELEM((*I)(), i, j)) < threshold * stddevbg)
	    {
		onepoint.x = j;
		onepoint.y = i;
		onepoint.z = MAT_ELEM((*I)(), i, j);
		onepoint.w = 1.;
		allpoints.push_back(onepoint);
	    }
	}
    }
    least_squares_plane_fit(allpoints, pA, pB, pC);

    // Substract the plane from the image
    FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)())
    {
        MAT_ELEM((*I)(), i, j) -= pA * j + pB * i + pC;
    }

    // Get std.dev. of the background pixels within +/- threshold*sigma 
    FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)())
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
	    if ( ABS(MAT_ELEM((*I)(), i, j)) < threshold * stddevbg)
	    {
		N++;
		sum1 +=  (double) MAT_ELEM((*I)(), i, j);
		sum2 += ((double) MAT_ELEM((*I)(), i, j)) * 
		    ((double) MAT_ELEM((*I)(), i, j));
	    }
	}
    }
    // average and standard deviation
    aux = sum1 / (double) N;
    newstddev = sqrt(ABS(sum2 / N - aux*aux) * N / (N - 1));

    // Replace pixels outside +/- threshold*sigma by samples from 
    // a gaussian with avg-plane and newstddev
    FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)())
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
	    if ( ABS(MAT_ELEM((*I)(), i, j)) > threshold * stddevbg)
	    {
		// get local average
		aux = pA * j + pB * i + pC;
		MAT_ELEM((*I)(), i, j)=rnd_gaus(aux, newstddev );
	    }
	}
    }

    // Divide the entire image by the new background
    (*I)() /= newstddev;

}

#ifdef NEVER_DEFINED
// This version doesn't work because of the high variance of avg-avg_bg
void normalize_NewXmipp(Image *I, const Matrix2D<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    (*I)().computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
                                     stddevbg);
    (*I)() -= avgbg;
    (*I)() /= avg - avgbg;
}
#endif
