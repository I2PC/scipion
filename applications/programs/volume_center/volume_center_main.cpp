/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
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

#include <data/mask.h>
#include <data/xmipp_program.h>

class ProgVolumeCenter: public XmippMetadataProgram
{
public:
    Mask mask_prm;

    // Define parameters
    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("Finds the 3D center of mass and center the volume around this point");
        mask_prm.defineParams(this,INT_MASK,NULL,"Restrict the center of mass to the mask area.");
        addExampleLine("xmipp_volume_center -i volume.vol -o volumeCentered.vol");
    }

    // Read parameters
    void readParams()
    {
        XmippMetadataProgram::readParams();
        mask_prm.allowed_data_types = INT_MASK;
        if (checkParam("--mask"))
        	mask_prm.readParams(this);
    }

    // Show
    void show()
    {
    	if (verbose==0)
    		return;
        XmippMetadataProgram::show();
        mask_prm.show();
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        // Read input volume
        Image<double> volume;
        volume.read(fnImg);
        volume().setXmippOrigin();
        mask_prm.generate_mask(volume());

        // Compute center of mass
        Matrix1D<double> centerOfMass;
        volume().centerOfMass(centerOfMass, &mask_prm.get_binary_mask());

        // Move origin to that center of mass
        selfTranslate(BSPLINE3,volume(),-centerOfMass, DONT_WRAP);
        volume.write(fnImgOut);
    }
};

int main(int argc, char **argv)
{
	ProgVolumeCenter prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
