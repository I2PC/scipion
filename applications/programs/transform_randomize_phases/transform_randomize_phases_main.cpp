/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <data/xmipp_program.h>
#include <data/xmipp_fftw.h>

class ProgRandomizePhases: public XmippMetadataProgram
{
public:
	String freqMode;
	double w;
	double Ts;
    void defineParams()
    {
        each_image_produces_an_output = true;
        allow_apply_geo = false;
        XmippMetadataProgram::defineParams();
        addUsageLine("This program randomizes all phases beyond a certain frequency");
        addParamsLine("  --freq <mode=digital>                  : Frequency beyond which phases will be randomized");
        addParamsLine("        where <mode>");
        addParamsLine("              discrete <w=0.25>          : Discrete frequencies are normalized to 0.5");
        addParamsLine("              continuous <w=15> <s=1>    : Continuous frequencies are in Angstroms. s=Sampling rate (Angstroms/pixel), w=Freq. in Angstroms");
        addExampleLine("Randomize all phases beyond 15 A in a volume whose sampling rate is 2 A/pixel",false);
        addExampleLine("xmipp_transform_randomize_phases -i volume.vol --freq continuous 15 2");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        freqMode=getParam("--freq");
        w=getDoubleParam("--freq",1);
        if (freqMode == "continuous")
        {
        	Ts=getDoubleParam("--freq",2);
        	w=Ts/w;
        }
        if (w>0.5 || w<0)
        	REPORT_ERROR(ERR_ARG_BADCMDLINE,"Discrete frequency must be between 0 and 0.5");
    }

    void show()
    {
    	XmippMetadataProgram::show();
    	std::cout << "Discrete freq: " << w << std::endl;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
    	Image<double> I;
    	I.read(fnImg);
    	randomizePhases(I(),w);
        I.write(fnImgOut);
    }
};

RUN_XMIPP_PROGRAM(ProgRandomizePhases)
