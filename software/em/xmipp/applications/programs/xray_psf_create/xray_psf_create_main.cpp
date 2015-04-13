/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#include <data/psf_xr.h>
#include <data/xmipp_program.h>

class ProgPSFXrCreate: public XmippProgram
{

protected:
    /// Filename with the Microscope Parameters.
    FileName fnParam, fnPSF;
    /// Number of threads;
    int nThr;
    /// PSF
    XRayPSF psf;

    void defineParams()
    {
        //Usage
        addUsageLine("Create a volume with the 3D PSF of an X-ray microscope.");
        addUsageLine("A param file can be passed or directly setting the microscope parameters.");
        addUsageLine("The program generates a PSF volume file and its associated info file.");
        //See Also
        addSeeAlsoLine("xray_project");
        addParamsLine("[-i <psf_param_file>] : XRay-Microscope parameters file.");
        addParamsLine(" alias --input;");
        addParamsLine("[-o <output_name_file>]  : Name for output files. It creates a PSF volume file and a PSF parameters file.");
        addParamsLine(" alias --output;");
        psf.defineParams(this);
        // Examples
        addExampleLine("The parameters are in a file",false);
        addExampleLine("xmipp_xray_psf_create -i psf560.xmd -o psf560.vol");
        addExampleLine("The parameters are given in the command line",false);
        addExampleLine("xmipp_xray_psf_create -o psf900.vol -lambda 2.5 -zones 900");
        addExampleLine("In the following link you can find an example of X-ray microscope parameters file:",false);
        addExampleLine(" ",false);
        addExampleLine("http://sourceforge.net/p/testxmipp/code/ci/3.0/tree/input/xray_psf.xmd",false);
    }

    void readParams()
    {
        psf.verbose = verbose;

        if (checkParam("-i"))
        {
            fnParam = getParam("-i");
            /* This forces always the creation of the PSF Volume file, even if the input already includes the
             name of a volume in the "image" label
             */
            psf.read(fnParam, false);
        }
        else
            psf.readParams(this);

        if (checkParam("-o"))
            fnPSF = getParam("-o");
        else
            fnPSF = fnParam;
    }

public:

    void run()
    {
        psf.generatePSF();
        psf.write(fnPSF);
    }
};

int main(int argc, char *argv[])
{
    ProgPSFXrCreate program;
    program.read(argc, argv);
    return program.tryRun();
}
