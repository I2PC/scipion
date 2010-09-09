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

class ProgHeaderAssign: public ProgHeader
{
private:
    double          rot, tilt, psi, xshift, yshift, weight;
    bool            mirror, round_shifts;
    int             levels, labelsnumber;
    std::vector<MDLabel> activeLabels;

protected:
    void defineParams()
    {
        addUsageLine("Set the geometric transformation (angles & shifts) in the header of 2D-images.");
        addParamsLine("   -i <metadata>      :Metadata file with input\n");
        addParamsLine("   alias --input;");
        addParamsLine("   [-round_shifts]    :Round shifts to integers");
        addParamsLine("   [-levels <n=0>]    :Levels of pyramidal reduction, n=1, 2, ...");
    }

    void readParams()
    {
        ProgHeader::readParams();
        round_shifts = checkParam("-round_shifts");
        levels = getIntParam("-levels");
    }

    void preprocess()
    {
        std::vector<MDLabel> activeLabels = md_input.getActiveLabels();
        labelsnumber = activeLabels.size();
    }

    void postprocess()
    {
        md_input.write(fn_out);
    }

    void headerProcess(FileName &fn_img)
    {

        img.read(fn_img);
        for (int iter = 0; iter < labelsnumber; iter++)
        {
            switch (activeLabels[iter])
            {
            case MDL_ANGLEROT:
                md_input.getValue( MDL_ANGLEROT, rot);
                img.setRot(rot);
                break;
            case MDL_ANGLETILT:
                md_input.getValue( MDL_ANGLETILT, tilt);
                img.setTilt(tilt);
                break;
            case MDL_ANGLEPSI:
                md_input.getValue( MDL_ANGLEPSI, psi);
                img.setPsi(psi);
                break;
            case MDL_SHIFTX:
                md_input.getValue( MDL_SHIFTX, xshift);
                if (levels != 0)
                    xshift /= pow(2.0, levels);
                if (round_shifts)
                    xshift = (float)ROUND(xshift);
                img.setXoff(xshift);
                break;
            case MDL_SHIFTY:
                md_input.getValue(MDL_SHIFTY, yshift);
                if (levels != 0)
                    yshift /= pow(2.0, levels);
                if (round_shifts)
                    yshift = (float)ROUND(yshift);
                img.setYoff(yshift);
                break;
            case MDL_WEIGHT:
                md_input.getValue( MDL_WEIGHT, weight);
                img.setWeight(weight);
                break;
            case MDL_FLIP:
                md_input.getValue( MDL_FLIP, mirror);
                img.setFlip(mirror);
                break;
            default:
                break;
            }
        }
        img.write(fn_img);

    }
}
;// end of class ProgHeaderAssign

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  try
  {
      ProgHeaderAssign program;
      program.read(argc, argv);
      program.run();
  }
  catch (XmippError xe)
  {
      std::cerr << xe;
  }
}



