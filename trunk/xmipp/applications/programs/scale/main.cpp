/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2000)
               Roberto Marabini (added fourier option)
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

#include <data/image.h>
#include <data/args.h>
#include <data/metadata.h>
#include <data/fftw.h>

void Usage();

int main(int argc, char **argv)
{
    FileName        fn_input, fn_output, fn_oext, fn_in, fn_out;
    MetaData        SF, SF_out;
    Image<double>   image;
    int             zdim, ydim, xdim;
    double          factor=-1;
    bool            linear;
    bool            fourier;
    int             nThreads;
    Matrix2D< double > A(3, 3), B(4, 4);
    A.initIdentity();
    B.initIdentity();

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = getParameter(argc, argv, "-i", NULL, 1, "Scale: Input file not found");
        fn_out   = getParameter(argc, argv, "-o", "");
        fn_oext  = getParameter(argc, argv, "-oext", "");
        if (fn_input.isMetaData())
            SF.read(fn_input);
        if (checkParameter(argc,argv,"-factor"))
        {
            factor=textToFloat(getParameter(argc, argv, "-factor"));
            if (factor<=0)
                REPORT_ERROR(1,"Factor must be a positive number");
        }
        else
        {
            zdim = textToInteger(getParameter(argc, argv, "-zdim", "0"));
            ydim = textToInteger(getParameter(argc, argv, "-ydim", "0"));
            xdim = textToInteger(getParameter(argc, argv, "-xdim"));
        }
        linear = checkParameter(argc, argv, "-linear");
        fourier = checkParameter(argc, argv, "-fourier");

        if (ydim == 0)
            ydim = xdim;
        if (zdim == 0)
            zdim = xdim;
        nThreads=textToInteger(getParameter(argc, argv, "-thr", "1"));
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }

    try
    {
        // Scale a single image -------------------------------------------------
        if (!fn_input.isMetaData())
        {
            image.read(fn_input);
            image().setXmippOrigin();
            if (image().getDim()==2)
            {
                if (factor>0)
                {
                    ydim=YSIZE(image())*factor;
                    xdim=XSIZE(image())*factor;
                }

                if (fourier)
                    selfScaleToSizeFourier(ydim,xdim,image(),nThreads);
                else if (linear)
                    selfScaleToSize(LINEAR,image(),xdim, ydim);
                else
                    selfScaleToSize(BSPLINE3,image(),xdim, ydim);
            }
            else
            {
                if (factor>0)
                {
                    zdim=ZSIZE(image())*factor;
                    ydim=YSIZE(image())*factor;
                    xdim=XSIZE(image())*factor;
                }
                if (linear)
                    selfScaleToSize(LINEAR,image(),xdim, ydim, zdim);
                else
                    selfScaleToSize(BSPLINE3,image(),xdim, ydim, zdim);
            }
            if (fn_out == "")
                image.write(fn_input);
            else
                image.write(fn_out);
        }
        // Scale a selection file ------------------------------------------------
        else
        {
            SF.read(fn_input);

            // Initialize progress bar
            time_config();
            int i = 0;
            //Image<double> image;

            init_progress_bar(SF.size());
            FOR_ALL_OBJECTS_IN_METADATA(SF)
            {
            	 // FIXME: This should disappear, now provoking segfault
                SF.getValue(MDL_IMAGE,fn_in);
                if (fn_oext == "")
                    fn_out = fn_in;
                else
                    fn_out = fn_in.without_extension() + "." + fn_oext;

                image.read(fn_in);
                image().setXmippOrigin();
                if (image().getDim()==2)
                {
                    if (factor>0)
                    {
                        ydim=YSIZE(image())*factor;
                        xdim=XSIZE(image())*factor;
                    }

                    if (fourier)
                        selfScaleToSizeFourier(ydim,xdim,image(),nThreads);
                    else if (linear)
                        selfScaleToSize(LINEAR,image(),xdim, ydim);
                    else
                        selfScaleToSize(BSPLINE3,image(),xdim, ydim);
                }
                else
                {
                    if (factor>0)
                    {
                        zdim=ZSIZE(image())*factor;
                        ydim=YSIZE(image())*factor;
                        xdim=XSIZE(image())*factor;
                    }
                    if (linear)
                        selfScaleToSize(LINEAR,image(),xdim, ydim, zdim);
                    else
                        selfScaleToSize(BSPLINE3,image(),xdim, ydim, zdim);
                }
                image.write(fn_out);
                SF_out.addObject();
                SF_out.setValue(MDL_IMAGE,fn_out);

                if (i++ % 25 == 0)
                    progress_bar(i);
            }
            progress_bar(SF.size());

            SF_out.write((SF.getFilename()).insert_before_extension(fn_oext));
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cerr << "Purpose:\n";
    std::cerr << "    Scale images/volumes to a given size\n";

    std::cerr << "Usage: scale <parameters>\n"
    << "   -i <image or volume> [-o <image_out or volume_out]\n"
    << "   -i <selfile> [-oext <output extension>]\n"
    << "  [-xdim <new x dimension>]\n"
    << "  [-ydim <new y dimension=new x dimension>]\n"
    << "  [-zdim <new z dimension=new x dimension>]\n"
    << "  [-factor <scale factor>]\n"
    << "  [-linear]         : Use bilinear/trilinear interpolation\n"
    << "  [-fourier]        : Use padding/windowing in Fourier Space (only for 2D)\n"
    << "  [-thr n]          : Use n threads, only implemented for fourier interpolation\n";

}
