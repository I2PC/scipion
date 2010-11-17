/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
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

#include <data/mask.h>
#include <data/args.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/image.h>

#include <cstdio>

#include <data/program.h>


class ProgMask: public Mask_Params, public XmippProgram
{
public:

    Mask_Params     mask_prm;
    FileName        oext, fn_in, fn_out, fn_mask;
    int             save_mask;
    int             create_mask;
    int             count_above;
    double          th_above;
    int             count_below;
    double          th_below;
    bool            apply_geo;
    double          subs_val;
    std::string     str_subs_val;
    int             count;

    /// Metadata Blockname (several metadata may go in the same file)
    std::string blockName;

    //bool            _produces_an_output;
    //bool            _each_image_produces_an_output;

    //    ProgMask()
    //    {
    //        _produces_an_output=false;
    //        _each_image_produces_an_output=false;
    //
    //    }
    /* Read Parameters ----------------------------------------------------------------- */
    /* Define Parameters ----------------------------------------------------------------- */
    void defineParams()
    {
        addUsageLine("Apply a mask");
        addUsageLine("Example of use: Sample at circular mask inside radius 72");
        addUsageLine("   xmipp_mask  -i reference.vol -o output_volume.vol --mask circular -72");
        addUsageLine("As above but save mask");
        addUsageLine("   xmipp_mask  -i reference.vol --create_mask  output_volume.vol --mask circular -25");
        addUsageLine("mask a selection file");
        addUsageLine("   xmipp_mask  -i t7_10.sel -oext msk --mask circular -72");

        addParamsLine("   -i <input_volume>           : Input Volume");
        addParamsLine(" alias --input;");
        addParamsLine(" [--bn <blockName=\"\">]   : Block name for metadata file");
        addParamsLine(" alias --blockname;");

        addParamsLine("  [-o <root_file_name>]         : Root for output files");
        addParamsLine("  [-oext <extension=\"\">] :  Output file format extension.");
        addParamsLine("   [--save_mask]     : Apply and save mask");
        addParamsLine("   [--dont_apply_geo]    : Dont apply (opposite) geometric transformation as stored in header of 2D-image");
        addParamsLine("   [--create_mask <output_mask_file>]  : Don't apply and save mask");
        addParamsLine("   [--count_above <th>]                : Voxels within mask >= th");
        addParamsLine("   [--count_below <th>]                : Voxels within mask <= th");
        addParamsLine("   [--substitute <val=\"0\">]  : Value outside the mask: userProvidedValue|min|max|avg");
        commonCommandLine
    }
    void readParams()
    {
        fn_in = getParam("-i");
        blockName = getParam("--blockname");
        fn_out = checkParam("-o") ? getParam("-o") : fn_in;
        oext = getParam("-oext");

        save_mask    = checkParam("--save_mask");
        count_above  = checkParam("--count_above");
        apply_geo    = !checkParam("--dont_apply_geo");
        if (count_above)
            th_above  = getDoubleParam("-count_above");
        count_below  = checkParam("--count_below");
        if (count_below)
            th_below  = getDoubleParam("--count_below");
        create_mask  = checkParam("--create_mask");
        if (create_mask)
            fn_mask  = getParam("--create_mask");
        mask_prm.read(argc, argv);
        str_subs_val = getParam("--substitute");
        count = count_below || count_above;

        processCommonCommandLine
    }

    void readParamsMask()
{}

    /* RuMAIN -------------------------------------------------------------------- */
    void run()
    {
        MetaData        SF_in, SF_out;
        Image<double>   image;
        int             max_length = 0;

        // Read list of images --------------------------------------------------
        //This should be converted to the new xmipp 3.0 standard
        std::cerr << "fn_in " << fn_in << std::endl;
        if (!fn_in.isMetaData())
        {
            SF_in.addObject();
            SF_in.setValue(MDL_IMAGE,fn_in);
        }
        else
            SF_in.read(fn_in);
        int Nimg=SF_in.size();

        // Mask a selection file ------------------------------------------------
        if (create_mask && Nimg>1)
            REPORT_ERROR(ERR_MD_NOOBJ, "Mask: Cannot create a mask for a selection file\n");

        // Initialize progress bar
        time_config();
        int i = 0;
        if (!count && Nimg>1)
            init_progress_bar(Nimg);
        else
            max_length = MaxFileNameLength(SF_in);

        // Process all selfile
        FOR_ALL_OBJECTS_IN_METADATA(SF_in)
        {
            // In and out filenames ...........................................
            FileName fn_in;
            SF_in.getValue(MDL_IMAGE,fn_in);
            if (Nimg>1)
            {
                if (oext == "")
                    fn_out = fn_in;
                else
                {
                    fn_out = fn_in.withoutExtension() + "." + oext;
                    SF_out.addObject();
                    SF_out.setValue(MDL_IMAGE,fn_out);
                }
            }
            else
            {
                if (fn_out=="")
                    fn_out=fn_in;
            }

            // Read image
            image.read(fn_in);
            image().setXmippOrigin();

            // Generate mask
            if (ZSIZE(image())>1)
                apply_geo=false;
            if (apply_geo)
            {
                if (mask_prm.x0 + mask_prm.y0 != 0.)
                    REPORT_ERROR(ERR_ARG_INCORRECT, "Mask: -center option cannot be combined with apply_geo; use -dont_apply_geo");
                else
                    // Read geometric transformation from the image and store for mask
                    mask_prm.mask_geo = image.getTransformationMatrix();
            }
            mask_prm.generate_mask(image());

            // Apply mask
            if (!create_mask)
            {
                if      (str_subs_val=="min")
                    subs_val=image().computeMin();
                else if (str_subs_val=="max")
                    subs_val=image().computeMax();
                else if (str_subs_val=="avg")
                    subs_val=image().computeAvg();
                else
                    subs_val=textToFloat(str_subs_val);
                mask_prm.apply_mask(image(), image(), subs_val, apply_geo);
                if (!count)
                    if (fn_out == "")
                        image.write(fn_in);
                    else
                        image.write(fn_out);
            }
            else
                mask_prm.write_mask(fn_mask);

            // Count
            if (count)
            {
                if (mask_prm.datatype() == INT_MASK)
                {
                    int count;
                    std::string elem_type="pixels";
                    if (ZSIZE(image())>1)
                        elem_type="voxels";
                    if      (count_above && !count_below)
                    {
                        std::cout << stringToString(fn_in,max_length)
                        << " number of " << elem_type << " above " << th_above;
                        count=count_with_mask_above(mask_prm.get_binary_mask(),
                                                    image(),th_above);
                    }
                    else if (count_below && !count_above)
                    {
                        std::cout << stringToString(fn_in,max_length)
                        << " number of " << elem_type << " below " << th_below;
                        count=count_with_mask_below(mask_prm.get_binary_mask(),
                                                    image(),th_below);
                    }
                    else if (count_below && count_above)
                    {
                        std::cout << stringToString(fn_in,max_length)
                        << " number of " << elem_type << " above " << th_above
                        << " and below " << th_below << " = ";
                        count=count_with_mask_between(mask_prm.get_binary_mask(),
                                                      image(),th_above,th_below);
                    }
                    std::cout << " = " << count << std::endl;
                }
                else
                    std::cerr << "Cannot count pixels with a continuous mask\n";
            }

            if (i++ % 25 == 0 && !count)
                progress_bar(i);
        }
        if (!count)
            progress_bar(Nimg);
        if (oext != "" && Nimg>1)
        {
            fn_out = fn_in.insertBeforeExtension(oext);
            SF_out.write(fn_out);
        }

        // Save masks -----------------------------------------------------------
        if (save_mask)
            mask_prm.write_mask("mask.spi");

        exit(0);
    }


    void processImage()
{}

}
;


/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    try
    {
        ProgMask program;
        /*
         if (checkParameter(argc,argv,"--create_mask"))
            program._produces_an_output = true;
        else
            program._each_image_produces_an_output = true;
         */
        program.XmippProgram::read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
        return 1;
    }
    return 0;
}

