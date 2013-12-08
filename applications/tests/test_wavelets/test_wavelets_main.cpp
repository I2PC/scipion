/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/multidim_array.h>
#include <reconstruction/fringe_processing.h>
#include <data/xmipp_image.h>
#include <data/wavelet.h>
#include <iostream>
#include "../../../external/gtest-1.6.0/fused-src/gtest/gtest.h"

// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
class WaveletTests : public ::testing::Test
{
protected:
    //init metadatas
    virtual void SetUp()
    {
#define len 128
        //get example down1_42_Periodogramavg.psd
        if (chdir(((String)(getXmippPath() + (String)"/resources/test/filters")).c_str())==-1)
        	REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot change directory");
        Image<double> img;
        img.read("KLH.tif");
        im = img();
    }

    //Image to be processed:
    MultidimArray<double> im;

};

TEST_F(WaveletTests, phaseCongMono)
{
    MultidimArray< double > Or,Ph,Energy,lowPass,Radius;
    MultidimArray< std::complex <double> > H;
    int nScale = 2;
    double minWaveLength=80;
    double mult = 1.25;
    double sigmaOnf = 2;

    phaseCongMono(im,Or,Ph,Energy,lowPass,Radius,H,nScale,minWaveLength,mult,sigmaOnf);
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


