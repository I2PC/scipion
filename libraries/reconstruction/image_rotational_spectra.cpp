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

#include "image_rotational_spectra.h"
#include <data/metadata_extension.h>

// Usage -------------------------------------------------------------------
void ProgMakeSpectra::defineParams()
{
	addUsageLine("Computes the rotational spectrum of a set of images.");
	addUsageLine("+This program generates a Fourier-Bessel decomposition of each image listed ");
	addUsageLine("+in an image collection. For each one, it creates a vector with the armonic ");
	addUsageLine("+energy percentage as a function of the radius. ");
	addUsageLine("+ ");
	addUsageLine("+The rotational spectrum is computed around the point (x0,y0), by default the center of the image.");
	addUsageLine("+The length of the rotational spectrum is defined by the lowest and highest harmonic.");
	addUsageLine("+Each harmonic is calculated by the integration of the image over a ring defined by r1 and r2.");
	addUsageLine("+ ");
	addUsageLine("+The rotational spectrum was defined in Crowther, R.A., and Amos, L.A. ");
	addUsageLine("+Harmonic analysis of electron microscopy images with rotational symmetry.");
	addUsageLine("+J. Mol. Biol.: 60, 123-130 (1971).");
	addSeeAlsoLine("classify_kerdensom, image_vectorize");
	addParamsLine("   -i <file>                   : Input image, selfile or stack");
	addParamsLine("   -o <metadata>               : Output vector metadata");
	addParamsLine("  [--r1 <lowRadius=15>]        : Lowest Integration radius (as a percentage of the total radius)");
    addParamsLine("  [--r2 <highRadius=80>]       : Highest Integration radius (as a percentage of the total radius)");
    addParamsLine("  [--x0 <xcenter=-1>]          : In physical units.");
    addParamsLine("  [--y0 <ycenter=-1>]          : By default, the Xmipp origin");
    addParamsLine("  [--low  <lowerHarmonic=  1>] : Lower harmonic to compute");
    addParamsLine("  [--high <higherHarmonic=15>] : Higher harmonic to compute");
    addExampleLine("xmipp_image_rotational_spectra -i images.stk -o spectra.xmd");
}

// Read from command line --------------------------------------------------
void ProgMakeSpectra::readParams()
{
	fn_in = getParam("-i");
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
	if (verbose==0)
		return;
	std::cout
	<< "Input file: " << fn_in << std::endl
	<< "Output file: " << fn_out << std::endl
	;
    std::cout << rot_spt << std::endl;
}

// Run -----------------------------------------------------------------
void ProgMakeSpectra::run()
{
    MetaData vectorContent, vectorHeader;
    vectorHeader.setColumnFormat(false);
    bool first=true;
    size_t order=0;
    FileName fnImg;
    Image<double> I;
    MetaData MD;
    MD.read(fn_in);
    std::ofstream fhOutRaw;
    MultidimArray<float> spectrum;
    FileName fnOutRaw=fn_out.withoutExtension()+".vec";
    fnOutRaw.deleteFile();

    // Convert the radii from percentages to actual pixels
    size_t Xdim, Ydim, Zdim, Ndim;
    getImageSize(MD,Xdim,Ydim,Zdim,Ndim);
    rot_spt.rl=(int)((rot_spt.rl/100.0)*Xdim/2);
    rot_spt.rh=(int)((rot_spt.rh/100.0)*Xdim/2);

    MDRow row;
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        MD.getValue(MDL_IMAGE,fnImg,__iter.objId);
    	I.readApplyGeo(fnImg,MD,__iter.objId);
        rot_spt.compute_rotational_spectrum(I(), rot_spt.rl, rot_spt.rh,
                                            rot_spt.dr, rot_spt.rh - rot_spt.rl);
        fhOutRaw.open(fnOutRaw.c_str(),std::ios::app | std::ios::binary);
        if (!fhOutRaw)
        	REPORT_ERROR(ERR_IO_NOWRITE,fnOutRaw);
        typeCast(rot_spt.rot_spectrum,spectrum);
        fhOutRaw.write((char*)MULTIDIM_ARRAY(spectrum),XSIZE(spectrum)*sizeof(float));
        fhOutRaw.close();

        // Create header
        if (first)
        {
            size_t vectorSize=MULTIDIM_SIZE(rot_spt.rot_spectrum);
            size_t outId=vectorHeader.addObject();
            vectorHeader.setValue(MDL_XSIZE,XSIZE(I()),outId);
            vectorHeader.setValue(MDL_YSIZE,YSIZE(I()),outId);
            vectorHeader.setValue(MDL_ZSIZE,ZSIZE(I()),outId);
            vectorHeader.setValue(MDL_COUNT,MD.size(),outId);
            vectorHeader.setValue(MDL_CLASSIFICATION_DATA_SIZE,vectorSize,outId);
            vectorHeader.write(formatString("vectorHeader@%s",fn_out.c_str()));
            first=false;
        }

        // Save this image in the output metadata
        MD.getRow(row,__iter.objId);
        row.setValue(MDL_ORDER,order++);
        vectorContent.addRow(row);
    }
    vectorContent.write(formatString("vectorContent@%s",fn_out.c_str()),MD_APPEND);
}
