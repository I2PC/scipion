/***************************************************************************
 *
 * Authors:     Sjors Scheres
 *              Roberto Marabini
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
 *                                      <
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "em.h"

#include <fstream>
#include <iomanip>

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
#include <sstream>
#else
#include <strstream.h>
#endif

#define VERBOSE

void EM::write(const FileName &fn_out, const VolumeXmipp &V, bool reversed)
{
    FILE *fp;

    //fill em header and reverse if needed
    fill_header_from_xmippvolume(V, reversed);

    //open file
    if ((fp = fopen(fn_out.c_str(), "wb")) == NULL)
        REPORT_ERROR(1503, "EM::write: File " + fn_out + " cannot be saved");

    //write header. note that FWRITE can not be used because
    //floats and longs are invoved
    if (fwrite(&my_em_header, sizeof(char), SIZEOF_EM_HEADER, fp) !=
        SIZEOF_EM_HEADER)
        REPORT_ERROR(1503, "EM::write: Header of file " + fn_out + " cannot be saved");

    //data,
    float f; //only floats are suported
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(V())
    {
        f = (float) MULTIDIM_ELEM(V(), i);
        FWRITE(&f, sizeof(float), 1, fp, reversed);
    }

    fclose(fp);
}

/* ------------------------------------------------------------------------- */
void EM::read(const FileName &fn_in, VolumeXmipp &V, bool reversed)
{

    FILE *fp;
    // Check reversed
    reversed = read_header_from_file(fn_in, reversed);

    //open file
    if ((fp = fopen(fn_in.c_str(), "rb")) == NULL)
        REPORT_ERROR(1503, "EM::read: File " + fn_in + " cannot be read");

    V().resize((int)my_em_header.nz, (int)my_em_header.ny, (int)my_em_header.nx);
    //EM data is made by
    int mode_size;
    switch (my_em_header.mode)
    {
    case MODE_BYTE:
        mode_size = sizeof(unsigned char);
        break;
    case MODE_SHORT:
        mode_size = sizeof(short int);
        break;
    case MODE_LONG_INT:
        mode_size = sizeof(long int);
        break;
    case MODE_FLOAT:
        mode_size = sizeof(float);
        break;
    case MODE_COMPLEX:
        REPORT_ERROR(1503, "EM::read: Complex-like map not supported!");
        break;
    case MODE_DOUBLE:
        mode_size = sizeof(double);
        break;
    default:
        REPORT_ERROR(1503, "EM::read: I do not know how to read this em file \
                     format, try the reverse_endian flag");
        break;
    }

    // Get header size
    struct stat info;
    if (fstat(fileno(fp), &info))
        EXIT_ERROR(1, (std::string)"EM: Cannot get size of " + fn_in);
    int header_size = info.st_size - my_em_header.nx *
                      my_em_header.ny *
                      my_em_header.nz * mode_size;

    // Skip header
    fseek(fp, header_size, SEEK_SET);

    switch (my_em_header.mode)
    {
    case MODE_BYTE:
        unsigned char c;
        FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(V())
        {
            fread(&c, sizeof(unsigned char), 1, fp);
            MULTIDIM_ELEM(V(), i) = (double)c;
        }
        break;
    case MODE_SHORT:
        short int si;
        FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(V())
        {
            FREAD(&si, sizeof(short int), 1, fp, reversed);
            MULTIDIM_ELEM(V(), i) = (double)si;
        }
        break;
    case MODE_LONG_INT:
        long int li;
        FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(V())
        {
            FREAD(&li, sizeof(long int), 1, fp, reversed);
            MULTIDIM_ELEM(V(), i) = (double)li;
        }
        break;
    case MODE_FLOAT:
        float f;
        FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(V())
        {
            FREAD(&f, sizeof(float), 1, fp, reversed);
            MULTIDIM_ELEM(V(), i) = (double)f;
        }
        break;
    case MODE_DOUBLE:
        double d;
        FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(V())
        {
            FREAD(&d, sizeof(float), 1, fp, reversed);
            MULTIDIM_ELEM(V(), i) = (double)d;
        }
        break;
    default:
        REPORT_ERROR(1503, "EM::read: I do not know how to read this em file format");
        break;
    }
    fclose(fp);
}


/* ------------------------------------------------------------------------- */
void EM::clear()
{
    memset(&my_em_header, '\0', SIZEOF_EM_HEADER);
}


/* ------------------------------------------------------------------------- */
/** Fill em header from xmipp image. */
void EM::fill_header_from_xmippvolume(VolumeXmipp V, bool reversed)
{

    clear();
    my_em_header.mode  = MODE_FLOAT;
    if (IsLittleEndian())
        my_em_header.machine = 6;
    else
        my_em_header.machine = 3;
    if (reversed == false)
    {
        my_em_header.nx = V().colNumber();
        my_em_header.ny = V().rowNumber();
        my_em_header.nz = V().sliceNumber();
    }
    else
    {
        my_em_header.nx    = V().colNumber();
        little22bigendian(my_em_header.nx);
        my_em_header.ny    = V().rowNumber();
        little22bigendian(my_em_header.ny);
        my_em_header.nz = V().sliceNumber();
        little22bigendian(my_em_header.nz);
    }
}

/* ------------------------------------------------------------------------- */
/** Fill em header from em file. */
bool EM::read_header_from_file(const FileName &fn_in, bool reversed)
{
    clear();
    FILE *fp;
    //open file
    if ((fp = fopen(fn_in.c_str(), "rb")) == NULL)
        REPORT_ERROR(1503, "EM::read_header_from_file: File " + fn_in + " cannot be read");

    FREAD(&(my_em_header.machine), 1, 1, fp, false);
    FREAD(&(my_em_header.general_use), 1, 1, fp, false);
    FREAD(&(my_em_header.not_used), 1, 1, fp, false);
    FREAD(&(my_em_header.mode), 1, 1, fp, false);

    if (my_em_header.machine != 6 && IsLittleEndian() && reversed == false)
        reversed = true;
    else if (my_em_header.machine == 6 && IsBigEndian() && reversed == false)
        reversed = true;

    FREAD(&(my_em_header.nx), 4, 1, fp, reversed);
    FREAD(&(my_em_header.ny), 4, 1, fp, reversed);
    FREAD(&(my_em_header.nz), 4, 1, fp, reversed);
    fclose(fp);

    return(reversed);
}

