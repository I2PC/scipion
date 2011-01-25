/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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

#include "make_spectra.h"
#include <data/args.h>

// Empty constructor -------------------------------------------------------
ProgMakeSpectra::ProgMakeSpectra()
{
	produces_an_output=true;
	each_image_produces_an_output=false;
}

// Usage -------------------------------------------------------------------
void ProgMakeSpectra::defineParams()
{
	addUsageLine("Computes the rotational spectrum of a set of images");
    XmippMetadataProgram::defineParams();
	addParamsLine("   --r1 <lowRadius>            : Lowest Integration radius");
    addParamsLine("   --r2 <highRadius>           : Highest Integration radius");
    addParamsLine("  [--x0 <xcenter=-1>]          : In physical units.");
    addParamsLine("  [--y0 <ycenter=-1>]          : By default, the Xmipp origin");
    addParamsLine("  [--low  <lowerHarmonic=  1>] : Lower harmonic to compute");
    addParamsLine("  [--high <higherHarmonic=15>] : Higher harmonic to compute");
}

// Read from command line --------------------------------------------------
void ProgMakeSpectra::readParams()
{
	XmippMetadataProgram::readParams();
    fn_out = getParam("-o");
    rot_spt.rl = getIntParam("--r1");
    rot_spt.rh = getIntParam("--r2");
    rot_spt.dr = 1;
    rot_spt.x0 = getDoubleParam("--x0");
    rot_spt.y0 = getDoubleParam("--y0");
    rot_spt.numin = getIntParam("--low");
    rot_spt.numax = getIntParam("--high");
}

// Show --------------------------------------------------------------------
void ProgMakeSpectra::show()
{
	XmippMetadataProgram::show();
    std::cout << rot_spt << std::endl;
}

// Process an image --------------------------------------------------------
void ProgMakeSpectra::processImage(const FileName &fnImg, const FileName &fnImgOut,
	long int objId)
{
	Image<double> I;
	I.readApplyGeo(fnImg,mdIn,objId);
    rot_spt.compute_rotational_spectrum(I(), rot_spt.rl, rot_spt.rh,
                                        rot_spt.dr, rot_spt.rh - rot_spt.rl);
    Harmonics.push_back(rot_spt.rot_spectrum);
    Img_name.push_back(fnImg);
}

// Finish processing -------------------------------------------------------
void ProgMakeSpectra::postProcess()
{
    std::ofstream fh_out;
    fh_out.open(fn_out.c_str());

    if (!fh_out)
        REPORT_ERROR(ERR_IO_NOTOPEN, fn_out);
    if (Harmonics.size() != 0)
    {
        fh_out << XSIZE(Harmonics[0]) << " " << Harmonics.size() << std::endl;
        int imax = Harmonics.size();
        for (int i = 0; i < imax; i++)
        {
            double norm = Harmonics[i].sum() / 100.0;
            for (int j = 0; j < XSIZE(Harmonics[i]); j++)
                fh_out << floatToString(Harmonics[i](j) / norm, 6, 4) << " ";
            fh_out << Img_name[i] << std::endl;
        }
    }
    fh_out.close();
}
