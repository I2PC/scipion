/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#include <data/metadata.h>
#include <data/ctf.h>

#include <data/xmipp_program.h>

class ProgCtfCreateCtfdat : public XmippProgram
{

public:

    FileName fn_sel, fn_doc, fn_ctf, fn_param, oroot;
    bool do_mode2;

public:

    void readParams()
    {
        fn_sel = getParam("-i");
        oroot = getParam("--oroot");

        if ((do_mode2 = checkParam("--ctfs")))
        {
            fn_ctf = getParam("--ctfs");
        }
        else
        {
            fn_doc = getParam("--defocus", 0);
            fn_param = getParam("--defocus", 1);
        }
    }

    void defineParams()
    {
        addUsageLine("Create CTFdat files from a selfile that contains image selfiles for each micrograph.");
        addUsageLine("+This program is useful when interacting with other packages outputs. There are two ");
        addUsageLine("+main usage modes. In the first mode you should provide a ctfparam file with common");
        addUsageLine("+options of each micrograph and a metadata with defocus values for each one. ");
        addUsageLine("+In the second mode you should provide a metadata with filenames of the ctfparams for");
        addUsageLine("+each  micrograph.");

        addParamsLine("  -i <selfile>           : Input selfile of selfiles for each micrograph");
        addParamsLine("  --oroot <rootname=out>  : Root name for output files ");
        addParamsLine("  --defocus <metadata> <common_ctf> : CTFparam file with common parameters ");
        addParamsLine("                         : Metadata with defocus values for each micrograph ");
        addParamsLine("                         : this file may have either a single column (defocus)");
        addParamsLine("                         : or three columns (defocusU, defocusV, azimuth_angle)");
        addParamsLine("  or --ctfs <selfile>    : Selfile of CTF param files for each micrograph ");

        addExampleLine("Example of use: Sample at MODE 1", false);
        addExampleLine("   ctf_create_ctfdat -i input.sel --defocus defocus.doc input.ctfparam");
        addExampleLine("Example of use: Sample at MODE 2", false);
        addExampleLine("   ctf_create_ctfdat -i input.sel --ctfs ctfparams.sel");
    }


    void run()
    {
        MetaData mdIn, SFind, mdCtf, ctfdat;
        FileName fnsel, fnimg, fnctf;

        mdIn.read(fn_sel);

        if (do_mode2)//Mode 2: read directly the ctf file
        {
            mdCtf.read(fn_ctf);
        }
        else //Mode 1: create ctf file from defocus file and a common ctfparam file
        {
            // Write param files for each micrograph to disc and make an internal SFctf
            CTFDescription ctf;
            MetaData DFdef;
            double defU, defV, azi;
            int ii = 0;
            DFdef.read(fn_doc);
            ctf.read(fn_param);
            ctf.enable_CTF = true;
            ctf.enable_CTFnoise = false;

            if (mdIn.size() != DFdef.size())
                REPORT_ERROR(ERR_IO_SIZE, "Sizes between input images (-i) and docfile(-doc) should be the same!!! ");

            FOR_ALL_OBJECTS_IN_METADATA(DFdef)
            {
                ii++;
                DFdef.getValue(MDL_CTF_DEFOCUSU, defU, __iter.objId);
                if (!DFdef.getValue(MDL_CTF_DEFOCUSU,defV,__iter.objId))
                    defV=defU;
                if (!DFdef.getValue(MDL_CTF_DEFOCUS_ANGLE,azi,__iter.objId))
                    azi=0;

                ctf.DeltafU = defU;
                ctf.DeltafV = defV;
                ctf.azimuthal_angle = azi;
                fnctf.compose(oroot, ii, "ctfparam");
                ctf.write(fnctf);
                mdCtf.setValue(MDL_CTF_MODEL, fnctf, mdCtf.addObject());
                if (verbose)
                  std::cout << formatString(" Saved CTF parameter file %s for micrograph number %d", fnctf.c_str(), ii) << std::endl;
            }
        }

        // For both modes
        if (mdIn.size() != mdCtf.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Selfiles of options -i and -ctfs have unequal number of entries! ");

        size_t id;

        FOR_ALL_OBJECTS_IN_METADATA2(mdIn, mdCtf)
        {
            mdIn.getValue(MDL_SELFILE,fnsel,__iter.objId);
            mdCtf.getValue(MDL_CTF_MODEL,fnctf,__iter2.objId);
            SFind.read(fnsel);
            FOR_ALL_OBJECTS_IN_METADATA(SFind)
            {
                SFind.getValue(MDL_IMAGE,fnimg,__iter.objId);
                id = ctfdat.addObject();
                ctfdat.setValue(MDL_IMAGE,fnimg, id);
                ctfdat.setValue(MDL_CTF_MODEL,fnctf, id);
            }
        }

        FileName fn = oroot + ".ctfdat";
        ctfdat.write(fn);
        if (verbose)
          std::cout << " Saved CTFdat file as " << fn << std::endl << " Done! "<< std::endl;
    }
};


int main(int argc, char **argv)
{
    ProgCtfCreateCtfdat program;
    program.read(argc, argv);
    return program.tryRun();
}

