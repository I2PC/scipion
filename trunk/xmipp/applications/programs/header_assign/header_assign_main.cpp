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

class ProgHeaderAssign: public XmippMetadataProgram
{
private:
    double          rot, tilt, psi, xshift, yshift, weight,scale;
    bool            mirror, round_shifts;
    int             levels, labelsnumber;
    std::vector<MDLabel> activeLabels;

protected:
    void defineParams()
    {
        each_image_produces_an_output = true;
        apply_geo = false;
        addUsageLine("Assign rotation angles, origin offsets (and optionally weights and mirror flags)");
        addUsageLine("read from input file to the headers of images. ");

        addExampleLine("Example of use: Assign to headers metadata contained in file in_file.doc",false);
        addExampleLine("   xmipp_header_assign -i in_file.doc");
        addExampleLine("Example of use: Assign to headers metadata contained in file in_file.doc with blockname named ONE",false);
        addExampleLine("The metadata assign is append to the file out_file.doc",false);
        addExampleLine("   xmipp_header_extract -i in_file.doc  --mode append -o out_file.doc:ONE");

        XmippMetadataProgram::defineParams();
        addParamsLine("   [--round_shifts]    :Round shifts to integers");
        addParamsLine("   [--levels <n=0>]    :Levels of pyramidal reduction, n=1, 2, ...");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        round_shifts = checkParam("--round_shifts");
        levels = getIntParam("--levels");
    }

    void preProcess()
    {
        activeLabels = mdIn.getActiveLabels();
        labelsnumber = activeLabels.size();
    }

    void postProcess()
    {
        if(fn_in != fn_out)
            mdIn.write(fn_out,mode);
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {

        Image<double> img;
        img.read(fnImg);
        for (int iter = 0; iter < labelsnumber; iter++)
        {
            switch (activeLabels[iter])
            {
            case MDL_ANGLEROT:
                mdIn.getValue( MDL_ANGLEROT, rot, objId);
                img.setRot(rot);
                break;
            case MDL_ANGLETILT:
                mdIn.getValue( MDL_ANGLETILT, tilt, objId);
                img.setTilt(tilt);
                break;
            case MDL_ANGLEPSI:
                mdIn.getValue( MDL_ANGLEPSI, psi, objId);
                img.setPsi(psi);
                break;
            case MDL_SHIFTX:
                mdIn.getValue( MDL_SHIFTX, xshift, objId);
                if (levels != 0)
                    xshift /= pow(2.0, levels);
                if (round_shifts)
                    xshift = (float)ROUND(xshift);
                img.setXoff(xshift);
                break;
            case MDL_SHIFTY:
                mdIn.getValue(MDL_SHIFTY, yshift, objId);
                if (levels != 0)
                    yshift /= pow(2.0, levels);
                if (round_shifts)
                    yshift = (float)ROUND(yshift);
                img.setYoff(yshift);
                break;
            case MDL_WEIGHT:
                mdIn.getValue( MDL_WEIGHT, weight, objId);
                img.setWeight(weight);
                break;
            case MDL_SCALE:
                mdIn.getValue( MDL_SCALE, scale, objId);
                img.setScale(scale);
                break;
            case MDL_FLIP:
                mdIn.getValue( MDL_FLIP, mirror, objId);
                img.setFlip(mirror);
                break;
            default:
                break;
            }
        }
        img.write(fnImg);
    }
}
;// end of class ProgHeaderAssign

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgHeaderAssign program;
    program.read(argc, argv);
    program.tryRun();
}



