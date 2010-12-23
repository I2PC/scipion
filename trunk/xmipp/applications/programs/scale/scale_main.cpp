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
#include <data/program.h>

class ProgScale: public XmippMetadataProgram
{
protected:
    int             zdim, ydim, xdim;
    double          factor;
    bool            linear;
    bool            fourier;
    int             nThreads;
    Matrix2D< double > A, B;//(3, 3), B(4, 4);

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("Scale images/volumes to a given size");

        addParamsLine(" -dim <x> <y=x> <z=x>      :New x,y and z dimensions");
        addParamsLine(" or -factor <scale_factor> :Scale by a factor");
        addParamsLine(" [-interp <interpolation_type=spline>] : Interpolation type to be used. ");
        addParamsLine("      where <interpolation_type>");
        addParamsLine("        spline          : Use spline interpolation");
        addParamsLine("        linear          : Use bilinear/trilinear interpolation");
        addParamsLine("        fourier <thr=1> : Use padding/windowing in Fourier Space (only for 2D)");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        if (checkParam("-factor"))
        {
            factor = getDoubleParam("-factor");
            //Some extra validations for factor
            if (factor <= 0)
                REPORT_ERROR(ERR_VALUE_INCORRECT,"Factor must be a positive number");
        }
        else
        {
            xdim = getDoubleParam("-dim", 0);
            ydim = (String(getParam("-dim", 1)) == "x") ? xdim : getDoubleParam("-dim", 1);
            zdim = (String(getParam("-dim", 2)) == "x") ? xdim : getDoubleParam("-dim", 2);
        }
        linear = (getParam("-interp") == "linear");
        fourier = (getParam("-interp") == "fourier");
        if (fourier)
            nThreads = getIntParam("-interp", "fourier");
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, long int objId)
    {
        Image<double> img; img.read(fnImg);
        img().setXmippOrigin();

        if (img().getDim()==2)
        {
            if (factor>0)
            {
                ydim=YSIZE(img())*factor;
                xdim=XSIZE(img())*factor;
            }

            if (fourier)
                selfScaleToSizeFourier(ydim,xdim,img(),nThreads);
            else if (linear)
                selfScaleToSize(LINEAR,img(),xdim, ydim);
            else
                selfScaleToSize(BSPLINE3,img(),xdim, ydim);
        }
        else
        {
            if (factor>0)
            {
                zdim=ZSIZE(img())*factor;
                ydim=YSIZE(img())*factor;
                xdim=XSIZE(img())*factor;
            }
            if (linear)
                selfScaleToSize(LINEAR,img(),xdim, ydim, zdim);
            else
                selfScaleToSize(BSPLINE3,img(),xdim, ydim, zdim);
        }

        img.write(fnImgOut);
    }

public:
    /** Constructor */
    ProgScale()
    {
        factor = -1;
        A.initIdentity(3);
        B.initIdentity(4);
    }
}
; //end of class ProgScale


int main(int argc, char **argv)
{
    ProgScale program;
    program.read(argc, argv);
    program.tryRun();
} //main

