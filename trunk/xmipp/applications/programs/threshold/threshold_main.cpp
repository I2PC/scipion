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

#include <data/image.h>
#include <data/mask.h>
#include <data/args.h>

void Usage();

#define SET_SUBS_VAL(I, subs_val, str_subs_val) \
    if      (str_subs_val=="min") subs_val=I.computeMin(); \
    else if (str_subs_val=="max") subs_val=I.computeMax(); \
    else                          subs_val=textToFloat(str_subs_val);

int main(int argc, char **argv)
{
    Image<double>   I, Vdist, Vlabel;
    FileName        fn_in, fn_out, fn_dist, fn_label;
    double          th_below, th_above, th;
    double          dmin, dmax;
    int             enable_th_below, enable_th_above, enable_th_abs_below;
    int             enable_dmin, enable_dmax;
    int             binarize;
    bool            enable_substitute;
    bool            enable_random_substitute;
    bool            apply_geo;
    double          new_val, old_val, avg_val, sig_val;
    std::string     str_new_val, str_old_val;
    double          accuracy;
    MaskProgram     mask_prm(INT_MASK);
    MultidimArray<int> * mask;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_in      = getParameter(argc, argv, "-i");
        fn_out     = getParameter(argc, argv, "-o", fn_in.c_str());
        fn_dist    = getParameter(argc, argv, "-dist", "");
        fn_label   = getParameter(argc, argv, "-label", "");
        if ((enable_th_abs_below = checkParameter(argc, argv, "-abs_below")))
            th_below = textToFloat(getParameter(argc, argv, "-abs_below"));
        if ((enable_th_below = checkParameter(argc, argv, "-below")))
            th_below = textToFloat(getParameter(argc, argv, "-below"));
        if ((enable_th_above = checkParameter(argc, argv, "-above")))
            th_above = textToFloat(getParameter(argc, argv, "-above"));
        if ((enable_dmin = checkParameter(argc, argv, "-dmin")))
            dmin = textToFloat(getParameter(argc, argv, "-dmin"));
        if ((enable_dmax = checkParameter(argc, argv, "-dmax")))
            dmax = textToFloat(getParameter(argc, argv, "-dmax"));
        binarize = checkParameter(argc, argv, "-binarize");
        if (binarize)
        {
            if (enable_th_above)
                th = th_above;
            else if (enable_th_below)
                th = th_below;
            else if (enable_th_abs_below)
                th = th_below;
            else
                th = 0;
        }
        int i;
        if ((i = paremeterPosition(argc, argv, "-substitute")) != -1)
        {
            enable_substitute = true;
            if (i + 2 >= argc)
                REPORT_ERROR(ERR_ARG_DEPENDENCE, "Threshold: Not enough parameters behind -substitute\n");
            str_old_val = argv[i+1];
            str_new_val = argv[i+2];
            accuracy = textToFloat(getParameter(argc, argv, "-accuracy", "0"));
        }
        else
            enable_substitute = false;
        if ((i = paremeterPosition(argc, argv, "-random_substitute")) != -1)
        {
            enable_random_substitute = true;
            if (i + 3 >= argc)
                REPORT_ERROR(ERR_ARG_DEPENDENCE, "Threshold: Not enough parameters behind -substitute\n");
            old_val = textToFloat(argv[i+1]);
            avg_val = textToFloat(argv[i+2]);
            sig_val = textToFloat(argv[i+3]);
            accuracy = textToFloat(getParameter(argc, argv, "-accuracy", "0"));
        }
        else
            enable_random_substitute = false;
        // Read mask stuff
        mask_prm.read(argc, argv);
        apply_geo = !checkParameter(argc, argv, "-dont_apply_geo");

    }

    catch (XmippError Xe)
    {
        std::cout << Xe;
        Usage();
        mask_prm.usage();
        exit(1);
    }

    try
    {
        I.read(fn_in);
        if (ZSIZE(I())==1)
        {
            mask_prm.mask_geo = I.getTransformationMatrix();
            mask_prm.generate_mask(I(),apply_geo);
        }
        else
            mask_prm.generate_mask(I());
        mask = & (mask_prm.get_binary_mask());

        // Apply density restrictions -------------------------------------------
        if (enable_th_below)
            I().threshold("below", th_below, th_below, mask);
        if (enable_th_above)
            I().threshold("above", th_above, th_above, mask);
        if (enable_th_abs_below)
            I().threshold("abs_below", th_below, 0, mask);

        // Apply substitution ---------------------------------------------------
        if (enable_substitute)
        {
            SET_SUBS_VAL(I(), new_val, str_new_val);
            SET_SUBS_VAL(I(), old_val, str_old_val);
            I().substitute(old_val, new_val, accuracy, mask);
        }

        // Apply random substitution --------------------------------------------
        if (enable_random_substitute)
            I().randomSubstitute(old_val, avg_val, sig_val, accuracy, mask);

        // Apply distance restrictions ------------------------------------------
        if (ZSIZE(I())>1 && fn_dist != "")
        {
            Vdist.read(fn_dist);
            if (fn_label == "")
                REPORT_ERROR(ERR_ARG_MISSING, "Threshold: You must supply a label volume\n");
            Vlabel.read(fn_label);
            FOR_ALL_ELEMENTS_IN_ARRAY3D(MULTIDIM_ARRAY(Vdist))
            {
                if (Vlabel(k, i, j) > 0)
                    Vdist(k, i, j) *= -1;
                if (enable_dmin && Vdist(k, i, j) < dmin)
                    I(k, i, j) = 0;
                if (enable_dmax && Vdist(k, i, j) > dmax)
                    I(k, i, j) = 0;
            }
        }

        // Binarize -------------------------------------------------------------
        if (binarize)
            I().binarize(th);

        // Write output volume --------------------------------------------------
        I.write(fn_out);
    }
    catch (XmippError Xe)
    {
        std::cout << Xe;
    }
    exit(0);
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cout << "threshold [Parameters]:\n"
    << "   -i <File_in>                   : Input Xmipp Volume or Image\n"
    << "  [-o <File_out>]                 : If not given, the same as input\n"
    << "  [-dist <Distance volume>        : of the same size as the input\n"
    << "   -label <Label volume>          : of the same size as the input\n"
    << "  [-dmin <dmin>]                  : remove voxels whose distance is smaller\n"
    << "  [-dmax <dmax>]]                 : remove voxels whose distance is greater\n"
    << "  [-below <th>]                   : remove voxels below this threshold\n"
    << "  [-above <th>]                   : remove voxels above this threshold\n"
    << "  [-abs_below <th>]               : remove voxels below this threshold\n"
    << "                                  : they will be set to 0\n"
    << "  [-binarize]]                    : binarize output\n"
    << "  [-substitute <old_val> <new_val>: where a value can be\n"
    << "                                    (<val>|min|max)\n"
    << "  [-random_substitute <old_val> <avg_val> <std_val>: where avg_val and sig_val \n"
    << "                                    are the mean and stddev of a Gaussian distribution \n"
    << "    [-accuracy <accuracy=0>]]     : when substituting this value\n"
    << "                                    determines if two values are the same\n"
    << "  [-dont_apply_geo]               : dont apply (opposite) header transformation to mask\n";

}
