/***************************************************************************
 *
 * Authors:     Debora gil
 *              Roberto Marabini
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
 *                                      <
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "ccp4.h"

#include <data/args.h>

#include <fstream>
#include <iomanip>

#define VERBOSE
//#define DEBUG

/* ------------------------------------------------------------------------- */
void CCP4::write(const FileName &fn_out, const ImageXmipp &I, bool reversed,
                 double x_length, double y_length, double z_length)
{
    FILE *fp;

    //fill mrc header and reverse if needed
    fill_header_from_xmippimage(I, reversed, x_length, y_length, z_length);

    //open file
    if ((fp = fopen(fn_out.c_str(), "wb")) == NULL)
        REPORT_ERROR(1503, "CCP4::write: File " + fn_out + " cannot be saved");

    //write header. note that FWRITE can not be used because
    //floats and longs are invoved
    if (fwrite(&my_mrc_header, sizeof(char), SIZEOF_MRC_HEADER, fp) !=
        SIZEOF_MRC_HEADER)
        REPORT_ERROR(1503, "CCP4::write: Header of file " + fn_out + " cannot be saved");

    //data,
    float f; //only float are suported
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
    {
        f = (float) MAT_ELEM(I(), i, j);
        FWRITE(&f, sizeof(float), 1, fp, reversed);
    }

    fclose(fp);
}

void CCP4::write(const FileName &fn_out, const Tomogram &V, bool reversed,
    double x_length, double y_length, double z_length, bool isStack)
{
    FILE *fp;

    //fill mrc header and reverse if needed
    int Xdim, Ydim, Zdim;
    V.size(Xdim, Ydim, Zdim);
    float minval, maxval, avg, stddev;
    V.computeStats(minval, maxval, avg, stddev);
    fill_header3D(Xdim, Ydim, Zdim, minval, maxval, avg,
        reversed, x_length, y_length, z_length, isStack);

    //open file
    if ((fp = fopen(fn_out.c_str(), "wb")) == NULL)
        REPORT_ERROR(1503, "CCP4::write: File " + fn_out + " cannot be saved");

    //write header. note that FWRITE can not be used because
    //floats and longs are involved
    if (fwrite(&my_mrc_header, sizeof(char), SIZEOF_MRC_HEADER, fp) !=
        SIZEOF_MRC_HEADER)
        REPORT_ERROR(1503, "CCP4::write: Header of file " + fn_out + " cannot be saved");

    //data,
    for (int k=0; k<Zdim; k++)
        for (int i=0; i<Ydim; i++)
            for (int j=0; j<Xdim; j++)
            {
                float f = V(j,i,k);
                FWRITE(&f, sizeof(float), 1, fp, reversed);
            }

    fclose(fp);
}

void CCP4::write(const FileName &fn_out, SelFile &SF, bool reversed,
    double x_length, double y_length, double z_length,
    const FileName &fn_tilt)
{
    // Get angles and stack min, max, avg
    int Zdim, Ydim, Xdim;
    SF.ImgSize(Ydim,Xdim);
    Zdim=SF.ImgNo();
    int k=0;
    Matrix1D<double> tiltAngles;
    tiltAngles.resize(Zdim);
    double minval, maxval, avg;
    while (!SF.eof())
    {
        ImageXmipp I;
        I.read(SF.NextImg());
        tiltAngles(k)=I.tilt();
        double minvali, maxvali, avgi, stddevi;
        I().computeStats(avgi, stddevi, minvali, maxvali);

        if (k==0)
        {
            minval=minvali;
            maxval=maxvali;
            avg=avgi;
        }
        minval=XMIPP_MIN(minval,minvali);
        maxval=XMIPP_MAX(maxval,maxvali);
        avg+=avgi;

        k++;
    }
    avg/=Zdim;

    if (fn_tilt!="")
        tiltAngles.write(fn_tilt);

    // Now write

    //fill mrc header and reverse if needed
    fill_header3D(Xdim, Ydim, Zdim, minval, maxval, avg,
        reversed, x_length, y_length, z_length, true);

    //open file
    FILE *fp;
    if ((fp = fopen(fn_out.c_str(), "wb")) == NULL)
        REPORT_ERROR(1503, "CCP4::write: File " + fn_out + " cannot be saved");

    //write header. note that FWRITE can not be used because
    //floats and longs are involved
    if (fwrite(&my_mrc_header, sizeof(char), SIZEOF_MRC_HEADER, fp) !=
        SIZEOF_MRC_HEADER)
        REPORT_ERROR(1503, "CCP4::write: Header of file " + fn_out + " cannot be saved");

    // Write the data
    SF.go_first_ACTIVE();

    //data,
    while (!SF.eof())
    {
        ImageXmipp I;
        I.read(SF.NextImg());
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
        {
            float f=(float)I(i,j);
            FWRITE(&f, sizeof(float), 1, fp, reversed);
        }
    }
    fclose(fp);
}

/* ------------------------------------------------------------------------- */
void CCP4::read(const FileName &fn_in,
                ImageXmipp &I, bool reversed)
{

    FILE *fp;
    reversed = read_header_from_file(fn_in, reversed);
    //
    //open file
    if ((fp = fopen(fn_in.c_str(), "rb")) == NULL)
        REPORT_ERROR(1503, "CCP4::read: File " + fn_in + " cannot be read");
    if (my_mrc_header.mode != 2)
        REPORT_ERROR(1503, "Only float mode is implemented");
    I().resize(my_mrc_header.ny, my_mrc_header.nx);
    //MRC data is made by
    int mode_size;
    switch (my_mrc_header.mode)
    {
    case MODE_BYTE:
        mode_size = sizeof(unsigned char);
        break;
    case MODE_SHORT:
        mode_size = sizeof(short int);
        break;
    case MODE_FLOAT:
        mode_size = sizeof(float);
        break;
    default:
        REPORT_ERROR(1503, "CCP4::read: I do not know how to read this mrc file \
                     format, try the reverse_endian flag");
        break;
    }


    // Get header size
    struct stat info;
    if (fstat(fileno(fp), &info))
        EXIT_ERROR(1, (std::string)"CCP4: Cannot get size of " + fn_in);
    int header_size = info.st_size - my_mrc_header.nx *
                      my_mrc_header.ny *
                      my_mrc_header.nz * mode_size;

    // Skip header
    fseek(fp, header_size, SEEK_SET);

    switch (my_mrc_header.mode)
    {
    case MODE_BYTE:
        unsigned char c;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
        {
            fread(&c, sizeof(unsigned char), 1, fp);
            MAT_ELEM(I(), i, j) = (double)c;
        }
        break;
    case MODE_SHORT:
        short int si;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
        {
            FREAD(&si, sizeof(short int), 1, fp, reversed);
            MAT_ELEM(I(), i, j) = (double)si;
        }
        break;
    case MODE_FLOAT:
        float f;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
        {
            FREAD(&f, sizeof(float), 1, fp, reversed);
            MAT_ELEM(I(), i, j) = (double)f;
        }
        break;
    default:
        REPORT_ERROR(1503, "CCP4::read: I do not know how to read this mrc file format");
        break;
    }
    fclose(fp);
}

/* ------------------------------------------------------------------------- */
void CCP4::read(const FileName &fn_in,
                VolumeXmipp &V, bool reversed)
{

    FILE *fp;
    reversed = read_header_from_file(fn_in, reversed);
    //
    //open file
    if ((fp = fopen(fn_in.c_str(), "rb")) == NULL)
        REPORT_ERROR(1503, "CCP4::read: File " + fn_in + " cannot be read");

    V().resize(my_mrc_header.nz, my_mrc_header.ny, my_mrc_header.nx);
    //MRC data is made by
    int mode_size;
    switch (my_mrc_header.mode)
    {
    case MODE_BYTE:
        mode_size = sizeof(unsigned char);
        break;
    case MODE_SHORT:
        mode_size = sizeof(short int);
        break;
    case MODE_FLOAT:
        mode_size = sizeof(float);
        break;
    default:
        REPORT_ERROR(1503, "CCP4::read: I do not know how to read this mrc file \
                     format, try the reverse_endian flag");
        break;
    }


    // Get header size
    struct stat info;
    if (fstat(fileno(fp), &info))
        EXIT_ERROR(1, (std::string)"CCP4: Cannot get size of " + fn_in);
    int header_size = info.st_size - my_mrc_header.nx *
                      my_mrc_header.ny *
                      my_mrc_header.nz * mode_size;

    // Skip header
    fseek(fp, header_size, SEEK_SET);

    switch (my_mrc_header.mode)
    {
    case MODE_BYTE:
        unsigned char c;
        FOR_ALL_ELEMENTS_IN_MATRIX3D(V())
        {
            fread(&c, sizeof(unsigned char), 1, fp);
            VOL_ELEM(V(), k, i, j) = (double)c;
        }
        break;
    case MODE_SHORT:
        short int si;
        FOR_ALL_ELEMENTS_IN_MATRIX3D(V())
        {
            FREAD(&si, sizeof(short int), 1, fp, reversed);
            VOL_ELEM(V(), k, i, j) = (double)si;
        }
        break;
    case MODE_FLOAT:
        float f;
        FOR_ALL_ELEMENTS_IN_MATRIX3D(V())
        {
            FREAD(&f, sizeof(float), 1, fp, reversed);
            VOL_ELEM(V(), k, i, j) = (double)f;
        }
        break;
    default:
        REPORT_ERROR(1503, "CCP4::read: I do not know how to read this mrc file format");
        break;
    }
    fclose(fp);
}


/* ------------------------------------------------------------------------- */
void CCP4::clear()
{
    memset(&my_mrc_header, '\0', SIZEOF_MRC_HEADER);
} /*clear*/


/* ------------------------------------------------------------------------- */
/** Fill mrc header from xmipp image. */
void CCP4::fill_header_from_xmippimage(const ImageXmipp &I, bool reversed,
    double x_length, double y_length, double z_length)
{
    clear();
    if (IsLittleEndian())
    {
        (my_mrc_header.machst)[0]  = 0x4 << 2;
        (my_mrc_header.machst)[0] += 0x4 ;
        (my_mrc_header.machst)[1]  = 0x4 << 2;
        (my_mrc_header.machst)[1] += 0x4 ;
    }
    else /*Big Endian*/
    {
        (my_mrc_header.machst)[0]  = 0x1 << 2;
        (my_mrc_header.machst)[0] += 0x1 ;
        (my_mrc_header.machst)[1]  = 0x1 << 2;
        (my_mrc_header.machst)[1] += 0x1 ;
    }
    (my_mrc_header.map)[0] = 'M';
    (my_mrc_header.map)[1] = 'A';
    (my_mrc_header.map)[2] = 'P';
    (my_mrc_header.map)[3] = ' ';
    if (reversed == false)
    {
        my_mrc_header.xlen = my_mrc_header.nx = my_mrc_header.mx = I().colNumber();
        my_mrc_header.ylen = my_mrc_header.ny = my_mrc_header.my = I().rowNumber();
        my_mrc_header.zlen = my_mrc_header.nz = my_mrc_header.mz = 1;
        if (x_length != 0) my_mrc_header.xlen = x_length;
        if (y_length != 0) my_mrc_header.ylen = y_length;
        if (z_length != 0) my_mrc_header.zlen = z_length;

        my_mrc_header.mode  = MODE_FLOAT;
        my_mrc_header.mapc  = X_AXIS;
        my_mrc_header.mapr  = Y_AXIS;
        my_mrc_header.maps  = Z_AXIS;
        my_mrc_header.nxstart = -1 * int(my_mrc_header.nx / 2);
        my_mrc_header.nystart = -1 * int(my_mrc_header.ny / 2);
        my_mrc_header.nzstart = -1 * int(my_mrc_header.nz / 2);
        my_mrc_header.amin  = (float)(I().computeMin());
        my_mrc_header.amax  = (float)(I().computeMax());
        my_mrc_header.amean = (float)(I().computeAvg());
    }
    else
    {
        my_mrc_header.nx    = I().colNumber();
        little22bigendian(my_mrc_header.nx);
        my_mrc_header.xlen = my_mrc_header.mx = my_mrc_header.nx;

        my_mrc_header.ny    = I().rowNumber();
        little22bigendian(my_mrc_header.ny);
        my_mrc_header.ylen = my_mrc_header.my = my_mrc_header.ny;

        my_mrc_header.nz = 1;
        little22bigendian(my_mrc_header.nz);
        my_mrc_header.zlen = my_mrc_header.mz = my_mrc_header.nz;

        if (x_length != 0)
        {
            my_mrc_header.xlen = x_length;
            little22bigendian(my_mrc_header.xlen);
        }
        if (y_length != 0)
        {
            my_mrc_header.ylen = y_length;
            little22bigendian(my_mrc_header.ylen);
        }
        if (z_length != 0)
        {
            my_mrc_header.zlen = z_length;
            little22bigendian(my_mrc_header.zlen);
        }

        my_mrc_header.nxstart = -1 * int(my_mrc_header.nx / 2);
        little22bigendian(my_mrc_header.nxstart);
        my_mrc_header.nystart = -1 * int(my_mrc_header.ny / 2);
        little22bigendian(my_mrc_header.nystart);
        my_mrc_header.nzstart = -1 * int(my_mrc_header.nz / 2);
        little22bigendian(my_mrc_header.nzstart);

        my_mrc_header.mode  = MODE_FLOAT;
        little22bigendian((my_mrc_header.mode));

        my_mrc_header.mapc  = X_AXIS;
        little22bigendian(my_mrc_header.mapc);
        my_mrc_header.mapr  = Y_AXIS;
        little22bigendian(my_mrc_header.mapr);
        my_mrc_header.maps  = Z_AXIS;
        little22bigendian(my_mrc_header.maps);

        my_mrc_header.amin = I().computeMin();
        little22bigendian(my_mrc_header.amin);
        my_mrc_header.amax = I().computeMax();
        little22bigendian(my_mrc_header.amax);
        my_mrc_header.amean = I().computeAvg();
        little22bigendian(my_mrc_header.amean);

    }
    my_mrc_header.alpha = my_mrc_header.beta = my_mrc_header.gamma = 90;
}

/* ------------------------------------------------------------------------- */
/** Fill mrc header from tomogram. */
void CCP4::fill_header3D(int Xdim, int Ydim, int Zdim,
    float minval, float maxval, float avg, bool reversed,
    double x_length, double y_length, double z_length, bool isStack)
{
    clear();
    (my_mrc_header.map)[0] = 'M';
    (my_mrc_header.map)[1] = 'A';
    (my_mrc_header.map)[2] = 'P';
    (my_mrc_header.map)[3] = ' ';

    if (reversed == false)
    {
        my_mrc_header.xlen = my_mrc_header.nx = my_mrc_header.mx = Xdim;
        my_mrc_header.ylen = my_mrc_header.ny = my_mrc_header.my = Ydim;
        my_mrc_header.zlen = my_mrc_header.nz = my_mrc_header.mz = Zdim;
        if (isStack) my_mrc_header.mz = 1;
        if (x_length != 0) my_mrc_header.xlen = x_length;
        if (y_length != 0) my_mrc_header.ylen = y_length;
        if (z_length != 0) my_mrc_header.zlen = z_length;
        my_mrc_header.mode  = MODE_FLOAT;
        my_mrc_header.mapc  = X_AXIS;
        my_mrc_header.mapr  = Y_AXIS;
        my_mrc_header.maps  = Z_AXIS;
        my_mrc_header.alpha =
            my_mrc_header.beta  =
                my_mrc_header.gamma = 90.;
        if (isStack)
        {
            my_mrc_header.nxstart = my_mrc_header.nystart =
                my_mrc_header.nzstart = 0;
        }
        else
        {
            my_mrc_header.nxstart = -1 * int(my_mrc_header.nx / 2);
            my_mrc_header.nystart = -1 * int(my_mrc_header.ny / 2);
            my_mrc_header.nzstart = -1 * int(my_mrc_header.nz / 2);
        }
        my_mrc_header.amin  = minval;
        my_mrc_header.amax  = maxval;
        my_mrc_header.amean = avg;
        if (IsLittleEndian())
        {
            (my_mrc_header.machst)[0]  = 0x4 << 4;
            (my_mrc_header.machst)[0] += 0x4 ;
            (my_mrc_header.machst)[1]  = 0x4 << 4;
            (my_mrc_header.machst)[1] += 0x4 ;
        }
        else /*Big Endian*/
        {
            (my_mrc_header.machst)[0]  = 0x1 << 4;
            (my_mrc_header.machst)[0] += 0x1 ;
            (my_mrc_header.machst)[1]  = 0x1 << 4;
            (my_mrc_header.machst)[1] += 0x1 ;
        }
    }
    else
    {
        my_mrc_header.nx    = Xdim;
        little22bigendian(my_mrc_header.nx);
        my_mrc_header.xlen = my_mrc_header.mx = my_mrc_header.nx;

        my_mrc_header.ny    = Ydim;
        little22bigendian(my_mrc_header.ny);
        my_mrc_header.ylen = my_mrc_header.my = my_mrc_header.ny;

        my_mrc_header.nz = Zdim;
        little22bigendian(my_mrc_header.nz);
        my_mrc_header.zlen = my_mrc_header.mz = my_mrc_header.nz;

        if (x_length != 0)
        {
            my_mrc_header.xlen = x_length;
            little22bigendian(my_mrc_header.xlen);
        }
        if (y_length != 0)
        {
            my_mrc_header.ylen = y_length;
            little22bigendian(my_mrc_header.ylen);
        }
        if (z_length != 0)
        {
            my_mrc_header.zlen = z_length;
            little22bigendian(my_mrc_header.zlen);
        }

        my_mrc_header.nxstart = -1 * int(my_mrc_header.nx / 2);
        little22bigendian(my_mrc_header.nxstart);
        my_mrc_header.nystart = -1 * int(my_mrc_header.ny / 2);
        little22bigendian(my_mrc_header.nystart);
        my_mrc_header.nzstart = -1 * int(my_mrc_header.nz / 2);
        little22bigendian(my_mrc_header.nzstart);

        my_mrc_header.mode  = MODE_FLOAT;
        little22bigendian((my_mrc_header.mode));

        my_mrc_header.mapc  = X_AXIS;
        little22bigendian(my_mrc_header.mapc);
        my_mrc_header.mapr  = Y_AXIS;
        little22bigendian(my_mrc_header.mapr);
        my_mrc_header.maps  = Z_AXIS;
        little22bigendian(my_mrc_header.maps);

        my_mrc_header.alpha =
            my_mrc_header.beta  =
                my_mrc_header.gamma = 90.;
        little22bigendian(my_mrc_header.alpha);
        little22bigendian(my_mrc_header.beta);
        little22bigendian(my_mrc_header.gamma);

        my_mrc_header.amin = minval;
        little22bigendian(my_mrc_header.amin);
        my_mrc_header.amax = maxval;
        little22bigendian(my_mrc_header.amax);
        my_mrc_header.amean = avg;
        little22bigendian(my_mrc_header.amean);

        if (IsLittleEndian())
        {
            (my_mrc_header.machst)[0]  = 0x1 << 4;
            (my_mrc_header.machst)[0] += 0x1 ;
            (my_mrc_header.machst)[1]  = 0x1 << 4;
            (my_mrc_header.machst)[1] += 0x1 ;
        }
        else /*Big Endian*/
        {
            (my_mrc_header.machst)[0]  = 0x4 << 4;
            (my_mrc_header.machst)[0] += 0x4 ;
            (my_mrc_header.machst)[1]  = 0x4 << 4;
            (my_mrc_header.machst)[1] += 0x4 ;
        }

    }
}

/* ------------------------------------------------------------------------- */
/** Fill mrc header from mrc file. */
bool CCP4::read_header_from_file(const FileName &fn_in, bool reversed)
{
    clear();
    FILE *fp;
    //open file
    if ((fp = fopen(fn_in.c_str(), "rb")) == NULL)
        REPORT_ERROR(1503, "CCP4::read_header_from_file: File " + fn_in + " cannot be saved");

    fseek(fp, 0xd0, SEEK_SET);
    FREAD(&(my_mrc_header.map), sizeof(unsigned char), 4, fp, reversed);
    FREAD(&(my_mrc_header.machst), sizeof(unsigned char), 4, fp, reversed);

    // Get  size
    if ((my_mrc_header.machst)[0] == ((0x1 << 4) + 0x1) &&
        IsLittleEndian() &&
        reversed == false)
        reversed = true;
    else if ((my_mrc_header.machst)[0] == ((0x4 << 4) + 0x4) &&
             IsBigEndian() &&
             reversed == false)
        reversed = true;
    fseek(fp, 0x00, SEEK_SET);
    FREAD(&(my_mrc_header.nx), sizeof(int), 1, fp, reversed);
    FREAD(&(my_mrc_header.ny), sizeof(int), 1, fp, reversed);
    FREAD(&(my_mrc_header.nz), sizeof(int), 1, fp, reversed);
    FREAD(&(my_mrc_header.mode), sizeof(int), 1, fp, reversed);

    fclose(fp);
    return(reversed);
}

