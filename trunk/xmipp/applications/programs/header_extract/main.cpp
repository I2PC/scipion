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

class ProgHeaderExtract: public ProgHeader
{
private:
    bool round_shifts;
    double           xx, yy;

protected:
    void defineParams()
    {
        addUsageLine("Extracts the geometric transformation (angles & shifts) in the header of 2D-images.");
        addParamsLine("   -i <metadata>      :metadata file with input.");
        addParamsLine("   alias --input;");
        addParamsLine("   [-o <metadata> ]    :metadata file, by default data is stored in input metaData file.");
        addParamsLine("   alias --output;");
        addParamsLine("   [-round_shifts]    :Round shifts to integers.");
    }

    void readParams()
    {
        ProgHeader::readParams();
        fn_out = checkParam("-o") ? getParam("-o") : fn_in;
        round_shifts = checkParam("-round_shifts");
    }

    void preprocess()
{}

    void postprocess()
    {
        md_input.write(fn_out);
    }

    void headerProcess(const FileName &fn_img)
    {

        img.read(fn_img, false);
        xx = (double) img.Xoff();
        yy = (double) img.Yoff();
        if (round_shifts)
        {
            xx = (double)ROUND(xx);
            yy = (double)ROUND(yy);
        }
        md_input.setValue(MDL_ANGLEROT,  (double) img.rot());
        md_input.setValue(MDL_ANGLETILT, (double) img.tilt());
        md_input.setValue(MDL_ANGLEPSI,  (double) img.psi());
        md_input.setValue(MDL_SHIFTX,    xx );
        md_input.setValue(MDL_SHIFTY,    yy );
        md_input.setValue(MDL_WEIGHT,    (double) img.weight());
        md_input.setValue(MDL_FLIP,      (bool)   img.flip());

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
