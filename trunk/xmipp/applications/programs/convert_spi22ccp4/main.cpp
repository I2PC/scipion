/***************************************************************************
 *
 * Authors:     Debora Gil
                Roberto Marabini
                Carlos Oscar Sánchez Sorzano
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
#include <data/volume.h>
#include <data/image.h>
#include <interface/ccp4.h>

#include <sys/stat.h>
#include <fstream>

/* Prototypes -============================================================= */
void Usage();

int main(int argc, char *argv[])
{
    /* Input Parameters ======================================================== */
    FileName       fn_in;    // input file
    FileName       fn_out;   // output file
    FileName       fn_tilt;  // File with tilt angles in case of a series
    double         x_length = 0;  // Cell Dimensions (Angstroms) for x-axis
    double         y_length = 0;  // Cell Dimensions (Angstroms) for y-axis
    double         z_length = 0;  // Cell Dimensions (Angstroms) for z-axis
    bool           reverse_endian;
    CCP4           mrcimage;

    /* Parameters ============================================================== */
    try
    {
        fn_in   = getParameter(argc, argv, "-i");
        fn_out  = getParameter(argc, argv, "-o");
        fn_tilt = getParameter(argc, argv, "-tilt_angles","");
        if (checkParameter(argc, argv, "-x_length"))
        {
            y_length = x_length = textToFloat(getParameter(argc, argv, "-x_length"));
            z_length = 1;
        }
        if (checkParameter(argc, argv, "-y_length"))
            y_length = textToFloat(getParameter(argc, argv, "-y_length"));
        if (checkParameter(argc, argv, "-z_length"))
            z_length = textToFloat(getParameter(argc, argv, "-z_length"));
        reverse_endian = checkParameter(argc, argv, "-reverse_endian");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }
    ImageXmipp   I;
    VolumeXmipp  V;

    try
    {
        if (Is_ImageXmipp(fn_in))
        {//is this a spider image?
            I.read(fn_in);
            mrcimage.write(fn_out, I, reverse_endian, x_length, y_length, z_length);
        }
        else if (Is_VolumeXmipp(fn_in))
        {//is this a spider volume?
            V.read(fn_in);
            mrcimage.write(fn_out, V, reverse_endian, x_length, y_length, z_length);
        }
        else if (fn_in.get_extension()=="sel")
        {// is this a selfile?
            SelFile SF;
            SF.read(fn_in);
            int Zdim, Ydim, Xdim;
            SF.ImgSize(Ydim,Xdim);
            Zdim=SF.ImgNo();
            V().initZeros(Zdim,Ydim,Xdim);
            int k=0;
            Matrix1D<double> tiltAngles;
            tiltAngles.resize(Zdim);
            while (!SF.eof())
            {
                I.read(SF.NextImg());
                V().setSlice(k,I());
                tiltAngles(k)=I.tilt();
                k++;
            }
            mrcimage.write(fn_out, V, reverse_endian, Xdim, Ydim, Zdim, true);
            if (fn_tilt!="")
                tiltAngles.write(fn_tilt);
        }
        else
        {
            mrcimage.read_header_from_file(fn_in, reverse_endian);
            if (mrcimage.my_mrc_header.nz > 1)
            {
                mrcimage.read(fn_in, V, reverse_endian);
                if (fn_tilt!="")
                {
                    std::ifstream fh_tilt;
                    fh_tilt.open(fn_tilt.c_str());
                    if (!fh_tilt)
                        REPORT_ERROR(1,(std::string)"Cannot read "+fn_tilt);
                    SelFile SF;
                    for (int k=0; k<ZSIZE(V()); k++)
                    {
                        float tilt;
                        if (fh_tilt.eof())
                            REPORT_ERROR(1,(std::string)"Not enough angles in "+
                                fn_tilt);
                        fh_tilt >> tilt;
                        ImageXmipp I;
                        V().getSlice(k,I());
                        I.set_tilt(tilt);
                        FileName fn_final;
                        fn_final.compose(fn_out, k, "xmp");
                        I.write(fn_final);
                        SF.insert(fn_final);
                    }
                    fh_tilt.close();
                    SF.write(fn_out+".sel");
                }
                else
                    V.write(fn_out);
            }
            else
            {
                mrcimage.read(fn_in, I, reverse_endian);
                I.write(fn_out);
            }
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
void Usage()
{
    printf(
        "Usage: ccp42spi [Purpose and Parameters]"
        "\nPurpose: Convert between CCP4 (map) and Spider/Xmipp format"
        "\n    -i   <file_in>           input CCP4/Xmipp file (2D or 3D)"
        "\n    -o   <file_out>          output Xmipp/CCP4 file. If CCP4 tilt series, rootname."
        "\n   [-tilt_angles <filename>] Filename with angles for a CCP4 tilt series"
        "\n   [-x_length  <length>]     Cell Dimensions (Angstroms) for x-axis"
        "\n   [-y_length  <length>]     Cell Dimensions (Angstroms) for y-axis. By default y_length=x_length"
        "\n   [-z_length  <length>]     Cell Dimensions (Angstroms) for z-axis. By default z_length=1"
        "\n   [-reverse_endian]         by default, output has the same endiannes as input"
        "\n                             use this option to change endianness\n");
}
