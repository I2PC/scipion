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
#include <data/xmipp_program.h>

class ProgPhantomCreate: public XmippProgram
{
protected:
    FileName          fn_phantom;
    FileName          fn_vol;
    Phantom           phantom;
    Image<double>     vol;

    void defineParams()
    {
    	addParamsLine("-i <description_file> : Input file with the mathematical features");
    	addParamsLine("-o <output_file>      : Output volume in voxels");
        addUsageLine("Create phantom volume from a feature description file with two Metadatas.");
        addUsageLine("+You may define a mathematical phantom from its geometrical features");
        addUsageLine("+(cubes, cylinders, ...) and translate it into a voxel volume with this program.");
        addUsageLine("+ File format can be found in this [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/FileFormats][link]].");
        addExampleLine("xmipp_phantom_create -i phantom.descr -o phantom.vol");
        addExampleLine("++In the following link you can find an example of phantom description file:",false);
        addExampleLine("++",false);
        addExampleLine("++http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/phantom_descr.descr",false);
    }

    void readParams()
    {
        fn_phantom = getParam("-i");
        fn_vol = getParam("-o");
    }

public:
    void run()
    {
        phantom.read(fn_phantom);
        phantom.draw_in(vol());
        vol.write(fn_vol);
    }
}
;//end of class ProgPhantomCreate

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgPhantomCreate program;
    program.read(argc, argv);
    return program.tryRun();
}
