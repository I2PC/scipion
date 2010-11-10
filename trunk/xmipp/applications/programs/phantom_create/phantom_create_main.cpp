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

#include <data/phantom.h>
#include <data/program.h>

class ProgPhantomCreate: public XmippProgram
{
protected:
    FileName          fn_phantom;
    FileName          fn_vol;
    Phantom           phantom;
    Image<double>     vol;

    void defineParams()
    {
        addUsageLine("Create phantom XMIPP volumes");
        addUsageLine("from a phantom feature description file.");

        addParamsLine("-i <description_file>");
        addParamsLine("-o <output_file>");
    }

    void readParams()
    {
        fn_phantom = getParam("-i");
        fn_vol = getParam("-o");
    }

public:
    void run()
    {
        // Read description file .............................................
        phantom.read(fn_phantom);
        // Generate volume and write .........................................
        phantom.draw_in(vol());
        vol.write(fn_vol);
    }
}
;//end of class ProgPhantomCreate

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    try
    {
        ProgPhantomCreate program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }
}
