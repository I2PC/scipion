/***************************************************************************
 *
 * Authors:    Sjors Scheres
 *             Slavica Jonic
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
 * MERCHANTABILITY or FITNESS FO A PARTICULAR PURPOSE.  See the
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

class ProgHeaderPrint: public ProgHeader
{
protected:
    void defineParams()
    {
        addUsageLine("Print information from the header of 2D-images.");

        addParamsLine(" -i <metadata>   :MetaData file with images or an individual image.");
        addParamsLine(" alias --input;");
    }

    void preprocess()
    {
        std::cout << " Printing the header ... " << std::endl;
    }

    void postprocess()
    {}

    void headerProcess(FileName &fn_img)
    {
        img.read(fn_img, false, -1, false);
        std::cout << img;
        std::cout << std::endl;
    }
}
;// end of class ProgHeaderPrint

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    try
    {
        ProgHeaderPrint program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }
}

