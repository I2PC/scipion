/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#include "progs.h"
#include "selfile.h"
#include "args.h"

/* Common functions -------------------------------------------------------- */
void Prog_parameters::read(int argc, char **argv)
{
    fn_in  = getParameter(argc, argv, "-i");
    fn_out = getParameter(argc, argv, "-o", "");
    oext   = getParameter(argc, argv, "-oext", "");
    oroot  = getParameter(argc, argv, "-oroot", "");
    // For each_image_produces_an_output there exists no possibility to apply_geo
    // This because it would require a back-transformation, which deteriorates the images
    if (!each_image_produces_an_output)
    {
        apply_geo = !checkParameter(argc, argv, "-dont_apply_geo");
    }
}

void Prog_parameters::show()
{
    cout << "Input File: " << fn_in << endl;
    if (apply_geo && !Is_VolumeXmipp(fn_in))
        cout << "Applying transformation stored in header of 2D-image" << endl;
    if (each_image_produces_an_output)
    {
        if (fn_out != "")
            cout << "Output File: " << fn_out << endl;
        if (oext != "")
            cout << "Output Extension: " << oext << endl;
        if (oroot != "")
            cout << "Output Root: " << oroot << endl;
    }
}

void Prog_parameters::usage()
{
    cerr << "   -i <input file>          : either an image/volume or a selection file\n";
    if (each_image_produces_an_output)
    {
        cerr << "  [-o <output file>]        : if wanted in case of a single image\n"
        << "  [-oext <extension>]       : if wanted in case of a selection file\n"
        << "  [-oroot <root>]           : if wanted in case of a selection file\n";
    }
    else
    {
        cerr  << "  [-dont_apply_geo]         : for 2D-images: do not apply transformation stored in the header\n";
    }
}

void Prog_parameters::get_input_size(int &Zdim, int &Ydim, int &Xdim)
{
    if (Is_ImageXmipp(fn_in))
    {
        ImageXmipp I(fn_in);
        Zdim = 1;
        Ydim = YSIZE(I());
        Xdim = XSIZE(I());
    }
    else if (Is_VolumeXmipp(fn_in))
    {
        VolumeXmipp V(fn_in);
        Zdim = ZSIZE(V());
        Ydim = YSIZE(V());
        Xdim = XSIZE(V());
    }
    else
    {
        SelFile SF;
        SF.read(fn_in);
        Zdim = 1;
        SF.ImgSize(Ydim, Xdim);
    }
}

int Prog_parameters::get_images_to_process()
{
    if (Is_ImageXmipp(fn_in))
        return 1;
    else if (Is_VolumeXmipp(fn_in))
        return 1;
    else
    {
        SelFile SF;
        SF.read(fn_in);
        return SF.ImgNo();
    }
}

/* With arguments ---------------------------------------------------------- */
void SF_main(int argc, char **argv,
             Prog_parameters *prm,
             void *process_img,
             void *process_vol, int operation_mode)
{
    // Some variables .......................................................
    ImageXmipp img, orig_img;
    VolumeXmipp vol;
    FourierImageXmipp IMG;
    FourierVolumeXmipp VOL;
    SelFile SF_in, SF_out;

    bool(*fi2i)(ImageXmipp &, const Prog_parameters *) = NULL;
    bool(*fi2I)(ImageXmipp &, FourierImageXmipp &, const Prog_parameters *) = NULL;
    bool(*fI2i)(FourierImageXmipp &, ImageXmipp &, const Prog_parameters *) = NULL;
    bool(*fI2I)(FourierImageXmipp &, const Prog_parameters *) = NULL;
    bool(*fi2F)(ImageXmipp &, const FileName &, const Prog_parameters *) = NULL;
    bool(*fv2v)(VolumeXmipp &, const Prog_parameters *) = NULL;
    bool(*fv2V)(VolumeXmipp &, FourierVolumeXmipp &, const Prog_parameters *) = NULL;
    bool(*fV2v)(FourierVolumeXmipp &, VolumeXmipp &, const Prog_parameters *) = NULL;
    bool(*fV2V)(FourierVolumeXmipp &, const Prog_parameters *) = NULL;
    bool(*fv2F)(VolumeXmipp &, const FileName &, const Prog_parameters *) = NULL;
    bool(*fF2F)(const FileName &, const FileName &, const Prog_parameters *) = NULL;

    // Read command line ....................................................
    try
    {
        (*(Prog_parameters *)prm).read(argc, argv);
    }
    catch (Xmipp_error XE)
    {
        cerr << XE;
        cerr << "Usage: \n"
        << "   " << argv[0] << " [options]\n";
        prm->usage();
        exit(1);
    }
    prm->show();

    try
    {
        if (!exists(prm->fn_in))
            EXIT_ERROR(1, (string)argv[0] + ": " + prm->fn_in + " doesn't exist");
        bool success;

        // For single image .....................................................
        FileName fn_out;
        if (prm->fn_out == "")
            fn_out = prm->fn_in;
        else
            fn_out = prm->fn_out;
        if (Is_ImageXmipp(prm->fn_in))
        {
            img.read(prm->fn_in, false, false, prm->apply_geo);
            img().setXmippOrigin();
            switch (operation_mode)
            {
            case IMAGE2IMAGE:
                fi2i = (bool(*)(ImageXmipp &, const Prog_parameters *)) process_img;
                success = fi2i(img, prm);
                if (prm->each_image_produces_an_output)
                {
                    /*       if (prm->apply_geo) {
                            cerr << "BUG: apply_geo and each_image_produces_an_output should not co-exist"
                                 << " the only exception is the apply_geo program";
                           }
                    */     img.write(fn_out);
                }
                break;
            case IMAGE2FOURIER:
                fi2I = (bool(*)(ImageXmipp &, FourierImageXmipp &, const Prog_parameters *)) process_img;
                success = fi2I(img, IMG, prm);
                if (prm->each_image_produces_an_output)
                    IMG.write(fn_out);
                break;
            case IMAGE2FILE:
                fi2F = (bool(*)(ImageXmipp &, const FileName &, const Prog_parameters *)) process_img;
                success = fi2F(img, fn_out, prm);
                break;
            }
            // For single volume ....................................................
        }
        else if (Is_VolumeXmipp(prm->fn_in))
        {
            vol.read(prm->fn_in);
            vol().setXmippOrigin();
            switch (operation_mode)
            {
            case IMAGE2IMAGE:
                fv2v = (bool(*)(VolumeXmipp &, const Prog_parameters *)) process_vol;
                success = fv2v(vol, prm);
                if (prm->each_image_produces_an_output)
                    vol.write(fn_out);
                break;
            case IMAGE2FOURIER:
                fv2V = (bool(*)(VolumeXmipp &, FourierVolumeXmipp &, const Prog_parameters *)) process_vol;
                success = fv2V(vol, VOL, prm);
                if (prm->each_image_produces_an_output)
                    VOL.write(fn_out);
                break;
            case IMAGE2FILE:
                fv2F = (bool(*)(VolumeXmipp &, const FileName &, const Prog_parameters *)) process_vol;
                success = fv2F(vol, fn_out, prm);
                break;
            }
            // For single Fourier image .............................................
        }
        else if (Is_FourierImageXmipp(prm->fn_in))
        {
            IMG.read(prm->fn_in);
            IMG().setXmippOrigin();
            switch (operation_mode)
            {
            case FOURIER2FOURIER:
                fI2I = (bool(*)(FourierImageXmipp &, const Prog_parameters *)) process_img;
                success = fI2I(IMG, prm);
                if (prm->each_image_produces_an_output)
                    IMG.write(fn_out);
                break;
            case FOURIER2IMAGE:
                fI2i = (bool(*)(FourierImageXmipp &, ImageXmipp &, const Prog_parameters *)) process_img;
                success = fI2i(IMG, img, prm);
                if (prm->each_image_produces_an_output)
                    img.write(fn_out);
                break;
            }
            // For single Fourier volume ............................................
        }
        else if (Is_FourierVolumeXmipp(prm->fn_in))
        {
            VOL.read(prm->fn_in);
            VOL().setXmippOrigin();
            switch (operation_mode)
            {
            case FOURIER2FOURIER:
                fV2V = (bool(*)(FourierVolumeXmipp &, const Prog_parameters *)) process_img;
                success = fV2V(VOL, prm);
                if (prm->each_image_produces_an_output)
                    VOL.write(fn_out);
                break;
            case FOURIER2IMAGE:
                fV2v = (bool(*)(FourierVolumeXmipp &, VolumeXmipp &, const Prog_parameters *)) process_img;
                success = fV2v(VOL, vol, prm);
                if (prm->each_image_produces_an_output)
                    vol.write(fn_out);
                break;
            }
        }
        else
        {
            // For a selection file .................................................
            SF_in.read(prm->fn_in);
            SF_out.clear();

            // Initialise progress bar
            time_config();
            int i = 0;
            if (prm->allow_time_bar)
                init_progress_bar(SF_in.ImgNo());
            int istep = CEIL((double)SF_in.ImgNo() / 60.0);

            // Process all selfile
            while (!SF_in.eof())
            {
                FileName fn_read;
                fn_read = SF_in.NextImg();
                if (prm->each_image_produces_an_output)
                    if (prm->oext == "" && prm->oroot == "")
                        prm->fn_out = fn_read;
                    else if (prm->oroot != "")
                    {
                        prm->fn_out = prm->oroot + fn_read.without_root();
                        if (operation_mode == IMAGE2FOURIER)
                            if (fn_read.get_extension() == "xmp")
                                prm->fn_out = prm->fn_out.without_extension() + ".fft";
                            else
                                prm->fn_out = prm->fn_out.without_extension() + ".fft3";
                        if (operation_mode == FOURIER2IMAGE)
                            if (fn_read.get_extension() == "fft")
                                prm->fn_out = prm->fn_out.without_extension() + ".xmp";
                            else
                                prm->fn_out = prm->fn_out.without_extension() + ".vol";
                    }
                    else
                    {
                        prm->fn_out = fn_read.without_extension() + "." + prm->oext;
                    }

                if (Is_ImageXmipp(fn_read))
                {
                    img.read(fn_read, false, false, prm->apply_geo);
                    img().setXmippOrigin();
                    switch (operation_mode)
                    {
                    case IMAGE2IMAGE:
                        fi2i = (bool(*)(ImageXmipp &, const Prog_parameters *)) process_img;
                        success = fi2i(img, prm);
                        if (prm->each_image_produces_an_output)
                        {
                            /*
                            if (prm->apply_geo) {
                              cerr << "BUG: apply_geo and each_image_produces_an_output should not co-exist"
                                   << " the only exception is the apply_geo program";
                            }
                            */
                            img.write(prm->fn_out);
                        }
                        break;
                    case IMAGE2FOURIER:
                        fi2I = (bool(*)(ImageXmipp &, FourierImageXmipp &, const Prog_parameters *)) process_img;
                        success = fi2I(img, IMG, prm);
                        if (prm->each_image_produces_an_output)
                            IMG.write(prm->fn_out);
                        break;
                    case IMAGE2FILE:
                        fi2F = (bool(*)(ImageXmipp &, const FileName &, const Prog_parameters *)) process_img;
                        success = fi2F(img, prm->fn_out, prm);
                        break;
                    }
                }
                else if (Is_VolumeXmipp(fn_read))
                {
                    vol.read(fn_read);
                    vol().setXmippOrigin();
                    switch (operation_mode)
                    {
                    case IMAGE2IMAGE:
                        fv2v = (bool(*)(VolumeXmipp &, const Prog_parameters *)) process_vol;
                        fv2v(vol, prm);
                        if (prm->each_image_produces_an_output)
                            vol.write(prm->fn_out);
                        break;
                    case IMAGE2FOURIER:
                        fv2V = (bool(*)(VolumeXmipp &, FourierVolumeXmipp &, const Prog_parameters *)) process_vol;
                        success = fv2V(vol, VOL, prm);
                        if (prm->each_image_produces_an_output)
                            VOL.write(prm->fn_out);
                        break;
                    case IMAGE2FILE:
                        fv2F = (bool(*)(VolumeXmipp &, const FileName &, const Prog_parameters *)) process_vol;
                        success = fv2F(vol, prm->fn_out, prm);
                        break;
                    }
                }
                else if (Is_FourierImageXmipp(fn_read))
                {
                    IMG.read(fn_read);
                    IMG().setXmippOrigin();
                    switch (operation_mode)
                    {
                    case FOURIER2FOURIER:
                        fI2I = (bool(*)(FourierImageXmipp &, const Prog_parameters *)) process_img;
                        success = fI2I(IMG, prm);
                        if (prm->each_image_produces_an_output)
                            IMG.write(prm->fn_out);
                        break;
                    case FOURIER2IMAGE:
                        fI2i = (bool(*)(FourierImageXmipp &, ImageXmipp &, const Prog_parameters *)) process_img;
                        success = fI2i(IMG, img, prm);
                        if (prm->each_image_produces_an_output)
                            img.write(prm->fn_out);
                        break;
                    }
                }
                else if (Is_FourierVolumeXmipp(fn_read))
                {
                    VOL.read(fn_read);
                    VOL().setXmippOrigin();
                    switch (operation_mode)
                    {
                    case FOURIER2FOURIER:
                        fV2V = (bool(*)(FourierVolumeXmipp &, const Prog_parameters *)) process_img;
                        success = fV2V(VOL, prm);
                        if (prm->each_image_produces_an_output)
                            VOL.write(prm->fn_out);
                        break;
                    case FOURIER2IMAGE:
                        fV2v = (bool(*)(FourierVolumeXmipp &, VolumeXmipp &, const Prog_parameters *)) process_img;
                        success = fV2v(VOL, vol, prm);
                        if (prm->each_image_produces_an_output)
                            vol.write(prm->fn_out);
                        break;
                    }
                }
                else if (operation_mode == FILE2FILE)
                {
                    fF2F = (bool(*)(const FileName &, const FileName &, const Prog_parameters *)) process_img;
                    success = fF2F(fn_read, prm->fn_out, prm);
                }

                if (prm->each_image_produces_an_output)
                    if (success)
                        SF_out.insert(prm->fn_out);
                    else
                        SF_out.insert(prm->fn_out, SelLine::DISCARDED);

                if (i++ % istep == 0 && prm->allow_time_bar)
                    progress_bar(i);
            }
            if (prm->allow_time_bar)
                progress_bar(SF_in.ImgNo());
            if (prm->each_image_produces_an_output)
                if (prm->oext != "")
                {
                    prm->fn_out = prm->fn_in.insert_before_extension(prm->oext);
                    SF_out.write(prm->fn_out);
                }
        }
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        exit(1);
    }
}
