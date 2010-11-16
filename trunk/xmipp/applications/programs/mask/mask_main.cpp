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


class ProgMask: public Mask_Params, public XmippMetadataProgram
{
public:

    Mask_Params     mask_prm;
    FileName        oext, fn_output, fn_input, fn_mask;
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

    /* Read Parameters ----------------------------------------------------------------- */
    /* Define Parameters ----------------------------------------------------------------- */
    void defineParams()
    {

        addUsageLine("Apply a mask");
        addUsageLine("Example of use: Sample at circular mask inside radius 72");
        addUsageLine("   xmipp_mask  -i reference.vol -o output_volume.vol -mask circular -72");

        //each_image_produces_an_output = true;
        produces_an_output = true;
        XmippMetadataProgram::defineParams();
        defineParamsMask();

        addParamsLine("   [--save_mask]     : Apply and save mask");
        addParamsLine("   [--dont_apply_geo]    : Dont apply (opposite) geometric transformation as stored in header of 2D-image");
        addParamsLine("   [--create_mask <output_mask_file>]  : Don't apply and save mask");
        addParamsLine("   [--count_above <th>]                : Voxels within mask >= th");
        addParamsLine("   [--count_below <th>]                : Voxels within mask <= th");
        //!a
        //addParamsLine("   [--substitute <val=0>]  : Value outside the mask: min|max|avg");
        addParamsLine("   [--substitute <val=min>]  : Value outside the mask: min|max|avg");
        addParamsLine("        where <val> min max avg");
    }

    //DEFINE PARAMETERS
    void defineParamsMask()
    {

        addParamsLine("   [--center <x0=0> <y0=0> <z0=0>]: Center of the mask");

        addParamsLine("   [--mask <mask_type=cirular>]");
        addParamsLine("        where <mask_type> ");
        addUsageLine("INT MASK");
        addParamsLine("         circular <R>        : circle/sphere mask");
        addParamsLine("                :if R>0 => outside R");
        addParamsLine("                :if R<0 => inside  R");
        addParamsLine("         DWT_circular <R> <smin> <smax>: circle/sphere mask");
        addParamsLine("                : smin and smax define the scales to be kept");
        addParamsLine("         rectangular <Xrect> <Yrect> <Zrect=0>: 2D or 3D rectangle");
        addParamsLine("                                         :if X,Y,Z > 0 => outside rectangle");
        addParamsLine("                                          :if X,Y,Z < 0 => inside rectangle");
        addParamsLine("         crown <R1> <R2>    : 2D or 3D crown");
        addParamsLine("                                          :if R1,R2 > 0 => outside crown");
        addParamsLine("                                          :if R1,R2 < 0 => inside crown");
        addParamsLine("         cylinder <R> <H>   : 2D circle or 3D cylinder");
        addParamsLine("                                         :if R,H > 0 => outside cylinder");
        addParamsLine("                                          :if R,H < 0 => inside cylinder");
        addParamsLine("         cone <theta>       : 3D cone (parallel to Z) ");
        addParamsLine("                                          :if theta > 0 => outside cone");
        addParamsLine("                                         :if theta < 0 => inside cone");
        addParamsLine("         wedge <th0> <thF>  : 3D missing-wedge mask for data ");
        addParamsLine("                                          :collected between tilting angles ");
        addParamsLine("                                          :th0 and thF (around the Y-axis) ");
        //addParamsLine("              |-mask <binary file>      : Read from file");

        addUsageLine("DOUBLE MASK");
        addParamsLine("         gaussian <sigma>   : 2D or 3D gaussian");
        addParamsLine("                              :if sigma > 0 => outside gaussian");
        addParamsLine("                              : if sigma < 0 => inside gaussian");
        addParamsLine("         raised_cosine <R1> <R2>: 2D or 3D raised_cosine");
        addParamsLine("                              : if R1,R2 > 0 => outside sphere");
        addParamsLine("                              : if R1,R2 < 0 => inside sphere");
        addParamsLine("         raised_crown <R1> <R2> <pixwidth>: 2D or 3D raised_crown");
        addParamsLine("                              : if R1,R2 > 0 => outside sphere");
        addParamsLine("                              : if R1,R2 < 0 => inside sphere");
        addParamsLine("         blob_circular <R1> <blob_radius>: 2D or 3D blob circular");
        addParamsLine("                              : if blob_radius > 0 => outside sphere");
        addParamsLine("                              : if blob_radius < 0 => inside sphere");
        addParamsLine("         blob_crown <R1> <R2> <blob_radius>: 2D or 3D blob_crown");
        addParamsLine("                              : if blob_radius > 0 => outside sphere");
        addParamsLine("                              : if blob_radius < 0 => inside sphere");
        //
        addParamsLine("         blackman           : 2D or 3D Blackman mask");
        addParamsLine("                             :  always inside blackman");
        addParamsLine("         sinc <w>          : 2D or 3D sincs");
        addParamsLine("                             :  if w > 0 => outside sinc");
        addParamsLine("                             :  if w < 0 => inside sinc");

        addParamsLine("   [ -m <blob_order=2>]       : Order of blob");
        //addParamsLine("     requires --xx;");
        //addParamsLine("     requires --blob_circular;");
        //addParamsLine("     requires --blob_crown;");

        addParamsLine("   [ -a <blob_alpha=10.4>]    : Alpha of blob");
        //addParamsLine("     requires blob_circular blob_crown;");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();

        fn_input = fn_in;
        fn_output = fn_out;

        readParamsMask();

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

    }

    void readParamsMask()
    {
        x0 = y0 = z0 = 0;
        x0 = getDoubleParam("--center",0);
        y0 = getDoubleParam("--center",1);
        z0 = getDoubleParam("--center",2);

        mask_type = getParam("--mask");

        // Circular mask ........................................................
        if (mask_type == "circular")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            R1 = getDoubleParam("--mask","circular");
            if (R1 < 0)
            {
                Mask_Params::mode = INNER_MASK;
                R1 = ABS(R1);
            }
            else if (R1 > 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: circular mask with radius 0");
            type = BINARY_CIRCULAR_MASK;
        }
        // Circular DWT mask ....................................................
        else if (mask_type == "DWT_circular")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            R1 = getDoubleParam("--mask","DWT_circular",0);
            smin = getIntParam("--mask","DWT_circular",1);
            smax = getIntParam("--mask","DWT_circular",2);
            quadrant = getParam("--mask","DWT_circular",3);
            type = BINARY_DWT_CIRCULAR_MASK;
        }
        // Rectangular mask .....................................................
        else if (mask_type == "rectangular")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            Xrect = getIntParam("--mask","rectangular",0);
            Yrect = getIntParam("--mask","rectangular",1);
            Zrect = getIntParam("--mask","rectangular",2);

            if (Xrect < 0 && Yrect < 0 && Zrect <= 0)
            {
                Mask_Params::mode = INNER_MASK;
                Xrect = ABS(Xrect);
                Yrect = ABS(Yrect);
                Zrect = ABS(Zrect);
            }
            else if (Xrect > 0 && Yrect > 0 && Zrect >= 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: cannot determine mode for rectangle");
            type = BINARY_FRAME_MASK;
        }
        // Cone mask ............................................................
        else if (mask_type == "cone")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            R1 = getDoubleParam("--mask","cone");
            if (R1 < 0)
            {
                Mask_Params::mode = INNER_MASK;
                R1 = ABS(R1);
            }
            else
                Mask_Params::mode = OUTSIDE_MASK;
            type = BINARY_CONE_MASK;
        }
        // Wedge mask ............................................................
        else if (mask_type == "wedge")
        {
            if (!(allowed_data_types & DOUBLE_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            R1 = getDoubleParam("--mask","wedge",0);
            R2 = getDoubleParam("--mask","wedge",1);
            type = BINARY_WEDGE_MASK;
        }
        // Crown mask ...........................................................
        else if (mask_type == "crown")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            R1 = getDoubleParam("--mask","crown",0);
            R2 = getDoubleParam("--mask","crown",1);
            if (R1 < 0 && R2 < 0)
            {
                Mask_Params::mode = INNER_MASK;
                R1 = ABS(R1);
                R2 = ABS(R2);
            }
            else if (R1 > 0 && R2 > 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: cannot determine mode for crown");
            type = BINARY_CROWN_MASK;
        }
        // Cylinder mask ........................................................
        else if (mask_type == "cylinder")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            R1 = getDoubleParam("--mask","cylinder",0);
            H = getDoubleParam("--mask","cylinder",1);
            if (R1 < 0 && H < 0)
            {
                Mask_Params::mode = INNER_MASK;
                R1 = ABS(R1);
                H = ABS(H);
            }
            else if (R1 > 0 && H > 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: cannot determine mode for cylinder");
            type = BINARY_CYLINDER_MASK;
        }
        // Gaussian mask ........................................................
        else if (mask_type == "gaussian")
        {
            if (!(allowed_data_types & DOUBLE_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: continuous masks are not allowed");
            sigma = getDoubleParam("--mask","gaussian");
            if (sigma < 0)
            {
                Mask_Params::mode = INNER_MASK;
                sigma = ABS(sigma);
            }
            else
                Mask_Params::mode = OUTSIDE_MASK;
            type = GAUSSIAN_MASK;
        }
        // Raised cosine mask ...................................................
        else if (mask_type == "raised_cosine")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: continuous masks are not allowed");
            R1 = getDoubleParam("--mask","raised_cosine",0);
            R2 = getDoubleParam("--mask","raised_cosine",1);
            if (R1 < 0 && R2 < 0)
            {
                Mask_Params::mode = INNER_MASK;
                R1 = ABS(R1);
                R2 = ABS(R2);
            }
            else if (R1 > 0 && R2 > 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: cannot determine mode for raised_cosine");
            type = RAISED_COSINE_MASK;
        }
        // Raised crown mask ....................................................
        else if (mask_type == "raised_crown")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: continuous masks are not allowed");
            R1 = getDoubleParam("--mask","raised_crown",0);
            R2 = getDoubleParam("--mask","raised_crown",1);
            pix_width = getDoubleParam("--mask","raised_crown",2);
            if (R1 < 0 && R2 < 0)
            {
                Mask_Params::mode = INNER_MASK;
                R1 = ABS(R1);
                R2 = ABS(R2);
            }
            else if (R1 > 0 && R2 > 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: cannot determine mode for raised_cosine");
            type = RAISED_CROWN_MASK;
        }
        // Blob circular mask ....................................................
        else if (mask_type == "blob_circular")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: continuous masks are not allowed");
            R1 = getDoubleParam("--mask","blob_circular",0);
            double aux = getDoubleParam("--mask","blob_circular",1);
            blob_radius= ABS(aux);
            if (aux < 0)
                Mask_Params::mode = INNER_MASK;
            else if (aux > 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: cannot determine mode for blob_circular");
            type = BLOB_CIRCULAR_MASK;

            blob_order= getDoubleParam("-m");
            blob_alpha= getDoubleParam("-a");

            // Raised crown mask ....................................................
        }
        else if (mask_type == "blob_crown")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: continuous masks are not allowed");
            R1 = getDoubleParam("--mask","blob_crown",0);
            R2 = getDoubleParam("--mask","blob_crown",1);
            double aux = getDoubleParam("--mask","blob_crown",2);
            blob_radius= ABS(aux);
            if (aux < 0)
                Mask_Params::mode = INNER_MASK;
            else if (aux > 0)
                Mask_Params::mode = OUTSIDE_MASK;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: cannot determine mode for blob_crown");
            type = BLOB_CROWN_MASK;
            blob_order= getDoubleParam("-m");
            blob_alpha= getDoubleParam("-a");

        }
        // Blackman mask ........................................................
        else if (mask_type == "blackman")
        {
            Mask_Params::mode = INNER_MASK;
            type = BLACKMAN_MASK;
        }
        // Sinc mask ............................................................
        else if (mask_type == "sinc")
        {
            if (!(allowed_data_types & INT_MASK))
                REPORT_ERROR(ERR_ARG_INCORRECT, "Mask_Params: binary masks are not allowed");
            omega = getDoubleParam("--mask","sinc");
            if (omega < 0)
            {
                Mask_Params::mode = INNER_MASK;
                omega = ABS(omega);
            }
            else
                Mask_Params::mode = OUTSIDE_MASK;
            type = SINC_MASK;
        }
        else
        {
            //fn_mask = argv[i+1];
            type = READ_MASK;
        }

    }

    /* RuMAIN -------------------------------------------------------------------- */
    void run()
    {
        MetaData        SF_in, SF_out;
        Image<double>   image;
        int             max_length = 0;

        // Read list of images --------------------------------------------------
        if (image.isImage(fn_input))
        {
            SF_in.addObject();
            SF_in.setValue(MDL_IMAGE,fn_input);
        }
        else
            SF_in.read(fn_input);
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
                        image.write(fn_input);
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
            fn_out = fn_input.insertBeforeExtension(oext);
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
        program.XmippMetadataProgram::read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
        return 1;
    }
    return 0;
}

