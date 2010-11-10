/***************************************************************************
 *
 * Authors:    Sjors Scheres
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

#include <data/args.h>
#include <data/image.h>
#include <data/metadata.h>
#include <data/progs.h>

/* PROGRAM ----------------------------------------------------------------- */

class ProgHeaderExtract: public XmippMetadataProgram
{
private:
    bool round_shifts;
    double           xx, yy;

protected:
    void defineParams()
    {
        produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("Extract the geometric transformation (angles & shifts) in the header of 2D-images");
        addUsageLine("and write values to a metadata file.");
        addParamsLine("   [--round_shifts]    :Round shifts to integers.");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        round_shifts = checkParam("--round_shifts");
    }

    void postProcess()
    {
        mdIn.write(fn_out);
    }

    void processImage()
    {
        img.read(fnImg, false);
        xx = img.Xoff();
        yy = img.Yoff();
        if (round_shifts)
        {
            xx = (double)ROUND(xx);
            yy = (double)ROUND(yy);
        }
        mdIn.setValue(MDL_ANGLEROT,  img.rot());
        mdIn.setValue(MDL_ANGLETILT, img.tilt());
        mdIn.setValue(MDL_ANGLEPSI,  img.psi());
        mdIn.setValue(MDL_SHIFTX,    xx );
        mdIn.setValue(MDL_SHIFTY,    yy );
        mdIn.setValue(MDL_WEIGHT,    img.weight());
        mdIn.setValue(MDL_SCALE,     img.scale());
        mdIn.setValue(MDL_FLIP,      img.flip());
    }
}
;// end of class ProgHeaderExtract

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    try
    {
        ProgHeaderExtract program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }
}
