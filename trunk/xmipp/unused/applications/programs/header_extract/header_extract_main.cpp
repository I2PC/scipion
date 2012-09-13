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

        addUsageLine("Extract the geometric transformation from header");
        addUsageLine("and write values to a metadata file.");
        addExampleLine("Example of use: Extract headers from in_file and overwrite out_file",false);
        addExampleLine("   xmipp_header_extract -i in_file.sel -o out_file.doc");
        addExampleLine("Example of use: Extract headers from in_file and append blockname named ONE to out_file",false);
        addExampleLine("   xmipp_header_extract -i in_file.sel --bn ONE --mode append -o out_file.doc");

        XmippMetadataProgram::defineParams();
        addParamsLine("   [--round_shifts]    :Round shifts to integers.");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        round_shifts = checkParam("--round_shifts");
    }

    void postProcess()
    {
        mdIn.write(fn_out, mode);
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.read(fnImg, _HEADER_ALL);
        //ApplyGeo(fnImg,mdIn,objId);
        xx = img.Xoff();
        yy = img.Yoff();
        if (round_shifts)
        {
            xx = (double)ROUND(xx);
            yy = (double)ROUND(yy);
        }
        mdIn.setValue(MDL_ANGLE_ROT,  img.rot(),objId);
        mdIn.setValue(MDL_ANGLE_TILT, img.tilt(),objId);
        mdIn.setValue(MDL_ANGLE_PSI,  img.psi(),objId);
        mdIn.setValue(MDL_SHIFT_X,    xx ,objId);
        mdIn.setValue(MDL_SHIFT_Y,    yy ,objId);
        mdIn.setValue(MDL_WEIGHT,    img.weight(),objId);
        mdIn.setValue(MDL_SCALE,     img.scale(),objId);
        mdIn.setValue(MDL_FLIP,      img.flip(),objId);
    }
}
;// end of class ProgHeaderExtract

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgHeaderExtract program;
    program.read(argc, argv);
    program.tryRun();
}
