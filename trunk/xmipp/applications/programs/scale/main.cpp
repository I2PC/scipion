/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2000)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/volume.h>
#include <data/image.h>
#include <data/args.h>
#include <data/selfile.h>
#include <data/gridding.h>
#include <data/fftw.h>

void Usage();

int main(int argc, char **argv)
{
    FileName        fn_input, fn_output, fn_oext, fn_in, fn_out;
    SelFile         SF, SF_out;
    ImageXmipp      image;
    VolumeXmipp     volume;
    int             zdim, ydim, xdim;
    double          factor=-1;
    bool            gridding;
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
        if (!Is_ImageXmipp(fn_input) && !Is_VolumeXmipp(fn_input))
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
        //gridding = checkParameter(argc, argv, "-gridding");
        linear = checkParameter(argc, argv, "-linear");
        fourier =   checkParameter(argc, argv, "-fourier");

        if (ydim == 0) ydim = xdim;
        if (zdim == 0) zdim = xdim;
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
        if (Is_ImageXmipp(fn_input))
        {
            image.read(fn_input);
            if (factor>0)
            {
                ydim=YSIZE(image())*factor;
                xdim=XSIZE(image())*factor;
            }
            
	    //if (gridding)
	    //{
		//KaiserBessel kb;
		//Matrix2D<double> Maux;
		//produceReverseGriddingMatrix2D(image(),Maux,kb);
		//DIRECT_MAT_ELEM(A, 0, 0) = (double) xdim / (double) XSIZE(image());
		//DIRECT_MAT_ELEM(A, 1, 1) = (double) ydim / (double) YSIZE(image());
		//applyGeometryReverseGridding(image(), A, Maux, kb, IS_NOT_INV, WRAP, xdim, ydim);
	    //}
            if (fourier)
            {
                selfScaleToSizeFourier(ydim,xdim,image(),nThreads);
            }
            else if (linear)
            {
		image().selfScaleToSize(ydim, xdim);
            }
	    else
	    {
		image().selfScaleToSizeBSpline(3, ydim, xdim);
	    }
            if (fn_out == "") image.write(fn_input);
            else            image.write(fn_out);

            // Scale a single volume ------------------------------------------------
        }
        else if (Is_VolumeXmipp(fn_input))
        {
            volume.read(fn_input);
            if (factor>0)
            {
                zdim=ZSIZE(volume())*factor;
                ydim=YSIZE(volume())*factor;
                xdim=XSIZE(volume())*factor;
            }
	    if (gridding)
	    {
		KaiserBessel kb;
		Matrix3D<double> Maux;
		produceReverseGriddingMatrix3D(volume(),Maux,kb);
		DIRECT_MAT_ELEM(B, 0, 0) = (double) xdim / (double) XSIZE(volume());
		DIRECT_MAT_ELEM(B, 1, 1) = (double) ydim / (double) YSIZE(volume());
		DIRECT_MAT_ELEM(B, 2, 2) = (double) zdim / (double) ZSIZE(volume());
		applyGeometryReverseGridding(volume(), B, Maux, kb, IS_NOT_INV, WRAP, xdim, ydim, zdim);
	    }
            else if (linear)
            {
		volume().selfScaleToSize(zdim, ydim, xdim);
            }
	    else
	    {
		volume().selfScaleToSizeBSpline(3, zdim, ydim, xdim);
	    }
            if (fn_out == "") volume.write(fn_input);
            else            volume.write(fn_out);

        }
	// Scale a selection file ------------------------------------------------
        else
        {
            SF.read(fn_input);

            SF_out.clear();
            // Initialise progress bar
            time_config();
            int i = 0;
            init_progress_bar(SF.ImgNo());
            while (!SF.eof())
            {
                fn_in = SF.NextImg();
                if (fn_in=="") break;
                if (fn_oext == "") fn_out = fn_in;
                else             fn_out = fn_in.without_extension() + "." + fn_oext;
                // Process an image ...............................................
                if (Is_ImageXmipp(fn_in))
                {
                    image.read(fn_in);
                    if (factor>0)
                    {
                        ydim=YSIZE(image())*factor;
                        xdim=XSIZE(image())*factor;
                    }
        	    if (fourier)
        	    {
                	selfScaleToSizeFourier(ydim,xdim,image(),nThreads);
        	    }
		    /*if (gridding)
		    {
			KaiserBessel kb;
			Matrix2D<double> Maux;
			produceReverseGriddingMatrix2D(image(),Maux,kb);
			DIRECT_MAT_ELEM(A, 0, 0) = (double) xdim / (double) XSIZE(image());
			DIRECT_MAT_ELEM(A, 1, 1) = (double) ydim / (double) YSIZE(image());
			applyGeometryReverseGridding(image(), A, Maux, kb, IS_NOT_INV, WRAP, xdim, ydim);
		    }*/
                    else if (linear)
                    {
			image().selfScaleToSize(ydim, xdim);
                    }
		    else
		    {
			image().selfScaleToSizeBSpline(3, ydim, xdim);
		    }
                    image.write(fn_out);
                    // Process a volume ...............................................
                }
                else if (Is_VolumeXmipp(fn_in))
                {
                    volume.read(fn_in);
                    if (factor>0)
                    {
                        zdim=ZSIZE(volume())*factor;
                        ydim=YSIZE(volume())*factor;
                        xdim=XSIZE(volume())*factor;
                    }
		    if (gridding)
		    {
			KaiserBessel kb;
			Matrix3D<double> Maux;
			produceReverseGriddingMatrix3D(volume(),Maux,kb);
			DIRECT_MAT_ELEM(B, 0, 0) = (double) xdim / (double) XSIZE(volume());
			DIRECT_MAT_ELEM(B, 1, 1) = (double) ydim / (double) YSIZE(volume());
			DIRECT_MAT_ELEM(B, 2, 2) = (double) zdim / (double) ZSIZE(volume());
			applyGeometryReverseGridding(volume(), B, Maux, kb, IS_NOT_INV, WRAP, xdim, ydim, zdim);
		    }
                    else if (linear)
                    {
			volume().selfScaleToSize(zdim, ydim, xdim);
                    }
		    else
		    {
			volume().selfScaleToSizeBSpline(3, zdim, ydim, xdim);
		    }
                    volume.write(fn_out);
                    // Not a Spider file ..............................................
                }
                else
                    std::cout << fn_in << " is not a SPIDER file\n";
                SF_out.insert(fn_out);

                if (i++ % 25 == 0) progress_bar(i);
            }
            progress_bar(SF.ImgNo());

            SF_out.write((SF.name()).insert_before_extension(fn_oext));
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
    //<< "  [-gridding]       : Use reverse gridding for interpolation\n"
    << "  [-linear]         : Use bilinear/trilinear interpolation\n"
    << "  [-fourier]        : Use padding/windowing in Fourier Space (only for 2D)\n"
    << "  [-thr n]          : Use n threads, only implemented for fourier interpolation\n";
    
}
