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

#include <data/argsparser.h>
#include <data/program.h>
#include <data/micrograph.h>
#include <data/args.h>

class ProgMicrographScissor: public XmippProgram
{
private:

protected:
    void defineParams()
    {
        addUsageLine ("Cut the images marked with xmipp_mark.");

        addParamsLine("  -i <input_untilted_micrograph>     : From which the untilted images will be cutted");
        addParamsLine("     alias --untilted;");
        addParamsLine("  [--orig <original_micrograph>]      : unless this parameter is specified");
        addParamsLine("  [-t <input_tilted_micrograph>]      : From which the   tilted images will be cutted");
        addParamsLine("     alias --tilted;");

        addParamsLine("  -o <output_stack>                   : Name for cutted images (stack fileName)");
        addParamsLine("     alias --untiltfn;");

        addParamsLine("  [--tiltfn <output_stack>]           : Name for tilted images (Stack FileName)");
        addParamsLine("     requires --untiltfn;                                                         ");

        addParamsLine("  [--pos <position_file>]             : file with particle coordinates (NOT for pairs)");
        addParamsLine("  [--down_transform <float=1.>]       : The transformation matrix was determined with this downsampling rate");

        addParamsLine("  --Xdim <window_X_dim>               : In pixels");
        addParamsLine("  [--Ydim <window_Y_dim>]             : If not given Ydim=Xdim");
        addParamsLine("  [--start <N=1>]                     : Number of the first image");
        addParamsLine("  [--invert]                          : Invert contrast");
        addParamsLine("  [--log]                             : Take logarithm (compute transmitance)");
        addParamsLine("  [--rmStack]                         : By default files are added to stack");

        addUsageLine ("Examples:");
        addUsageLine ("   xmipp_micrograph_scissor -i g7107.raw --pos g7107.raw.Common.pos -o kk.mrcs --Xdim 64");
        addUsageLine ("   xmipp_micrograph_scissor -i g7107.raw -t g7106.raw  --Xdim 64 -o u_.img --tiltfn t_.img");
    }
    FileName fn_micrograph,fn_root;
    FileName fn_tilted, fn_root_tilted;
    FileName fn_orig, fn_pos;
    int      startN;
    bool     pair_mode;
    int      Ydim, Xdim;
    bool     reverse_endian;
    bool     compute_transmitance;
    bool     compute_inverse ;
    double   down_transform;
    bool     rmStack;

    void readParams()
    {
        fn_micrograph = getParam  ("-i");
        pair_mode     = checkParam("--tilted");
        fn_root       = getParam  ("--untiltfn");
        Xdim          = getIntParam("--Xdim");
        if (checkParam("--Ydim"))
            Ydim      = getIntParam("--Ydim");
        else
            Ydim = Xdim;

        startN               = getIntParam("--start");
        compute_inverse      = checkParam("--invert");
        compute_transmitance = checkParam("--log");
        rmStack              = checkParam("--rmStack");

        if (!pair_mode)
        {
            fn_pos        = getParam("--pos");
            if (checkParam("--orig"))
            	fn_orig       = getParam("--orig");
        }
        else
        {
            fn_root_tilted = getParam("--tiltfn");
            fn_tilted      = getParam("--tilted");
        }
        down_transform = getDoubleParam("--down_transform");
    }
public:
    void run()
    {

        if (!pair_mode)
        {
            Micrograph m;
            m.open_micrograph(fn_micrograph);
            m.set_window_size(Xdim, Ydim);
            m.read_coordinates(0, fn_pos);
            if (down_transform != 1.)
                m.scale_coordinates(1./down_transform);
            m.add_label("");
            m.set_transmitance_flag(compute_transmitance);
            m.set_inverse_flag(compute_inverse);
            m.produce_all_images(0, fn_root, startN, fn_orig, 0.,0.,0., rmStack);
        }
        else
        {
            MetaData auxMd;
            // Read angles
            FileName fn_ang = fn_micrograph.withoutExtension();
            fn_ang = fn_ang.addExtension("ang");
            double alpha_u, alpha_t, tilt_angle;
            auxMd.read(fn_ang);
            auxMd.getValue(MDL_ANGLEPSI,alpha_u);
            auxMd.getValue(MDL_ANGLEPSI2,alpha_t);
            auxMd.getValue(MDL_ANGLETILT,tilt_angle);

            // Generate the images for the untilted image
            Micrograph m;
            m.open_micrograph(fn_micrograph);
            m.set_window_size(Xdim, Ydim);
            m.read_coordinates(0, fn_micrograph.addExtension("Common.pos"));
            m.add_label("");
            m.set_transmitance_flag(compute_transmitance);
            m.set_inverse_flag(compute_inverse);
            m.produce_all_images(0, fn_root, startN, "", alpha_u,0.,0.,rmStack);
            m.close_micrograph();

            // Generate the images for the tilted image
            Micrograph mt;
            mt.open_micrograph(fn_tilted);
            mt.set_window_size(Xdim, Ydim);
            mt.read_coordinates(0, fn_tilted.addExtension("Common.pos"));
            mt.add_label("");
            mt.set_transmitance_flag(compute_transmitance);
            mt.set_inverse_flag(compute_inverse);
            mt.produce_all_images(0, fn_root_tilted, startN, "", 0., tilt_angle, alpha_t, rmStack);
            mt.close_micrograph();
        }
    }
};

//#endif
void Usage();

int main(int argc, char **argv)
{
    try
    {
        ProgMicrographScissor program;
        program.read(argc, argv);
        program.run();

    }
    catch (XmippError e)
    {
        std::cerr << e.msg <<std::endl;
    }

}




