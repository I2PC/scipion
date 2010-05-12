/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2006)
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

#include "angular_assign_for_tomogram.h"
#include <data/projection.h>

#include <data/args.h>
#include <data/filters.h>

// Empty constructor =======================================================
AlignmentTomography::AlignmentTomography()
{
    rot=tilt=psi=x=y=corr=0;
    fn_img="";
}

// Read arguments ==========================================================
void Prog_angular_predict_tomography_prm::read(int argc, char **argv)
{
    fn_ref = getParameter(argc, argv, "-ref");
    fn_sel = getParameter(argc, argv, "-sel");
    fn_out = getParameter(argc, argv, "-oroot");
    max_rot_change = textToFloat(getParameter(argc, argv, "-max_rot_change", "5"));
    max_tilt_change = textToFloat(getParameter(argc, argv, "-max_tilt_change", "2"));
    max_psi_change = textToFloat(getParameter(argc, argv, "-max_psi_change", "5"));
    rot_step = textToFloat(getParameter(argc, argv, "-rot_step", "1"));
    tilt_step = textToFloat(getParameter(argc, argv, "-tilt_step", "1"));
    psi_step = textToFloat(getParameter(argc, argv, "-psi_step", "1"));
    max_shift_change = textToFloat(getParameter(argc, argv, "-max_shift_change", "10"));
    shift_step = textToFloat(getParameter(argc, argv, "-shift_step", "2"));
    adjustGray = checkParameter(argc, argv, "-adjustGray");
    produce_side_info();
}

// Show ====================================================================
void Prog_angular_predict_tomography_prm::show()
{
    std::cout << "Reference images:   " << fn_ref           << std::endl
    << "Input images:       " << fn_sel           << std::endl
    << "Ouput rootname:     " << fn_out           << std::endl
    << "Max rot change:     " << max_rot_change   << " step: " << rot_step << std::endl
    << "Max tilt change:    " << max_tilt_change  << " step: " << tilt_step << std::endl
    << "Max psi change:     " << max_psi_change   << " step: " << psi_step << std::endl
    << "Max shift change:   " << max_shift_change << " step: " << shift_step << std::endl
    ;
}

// usage ===================================================================
void Prog_angular_predict_tomography_prm::usage()
{
    std::cerr
    << "   -ref <volume>             : Reference volume\n"
    << "   -sel <selfile>            : Images to align\n"
    << "   -oroot <rootname>         : rootname for the output\n"
    << "  [-max_rot_change <ang=5>]  : Maximum change allowed in rot\n"
    << "  [-max_tilt_change <ang=2>] : Maximum change allowed in tilt\n"
    << "  [-max_psi_change <ang=5>]  : Maximum change allowed in psi\n"
    << "  [-max_shift_change <r=10>] : Maximum change allowed in shift\n"
    << "  [-rot_step <ang=1>]        : Rot search step\n"
    << "  [-tilt_step <ang=1>]       : Tilt search step\n"
    << "  [-psi_step <ang=3>]        : Psi search step\n"
    << "  [-shift_step <r=2>]        : Step in shift in pixels\n"
    << "  [-adjustGray]              : Adjust also gray values\n"
    ;
}

// Produce side information ================================================
void Prog_angular_predict_tomography_prm::produce_side_info()
{
    V.read(fn_ref);
    V().setXmippOrigin();
    MetaData SF(fn_sel);
    FileName fnAux;
    bool masksPresent=SF.getValue(MDL_MASK,fnAux);
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        AlignmentTomography dummy;
        Image<double> I, Imask;
        MultidimArray<int> mask;
        FileName fnImg;
        SF.getValue(MDL_IMAGE,fnImg);
        I.read(fnImg);
        I().setXmippOrigin();
        dummy.rot=I.rot();
        dummy.tilt=I.tilt();
        dummy.psi=I.psi();
        dummy.x=I.Xoff();
        dummy.y=I.Yoff();
        dummy.fn_img=I.name();
        Projection theo;
        project_Volume(V(), theo, YSIZE(V()), XSIZE(V()),
                       dummy.rot, dummy.tilt, dummy.psi);
        selfTranslate(LINEAR,I(),vectorR2(dummy.x,dummy.y),DONT_WRAP);
        if (masksPresent)
        {
            FileName fnMask;
            SF.getValue(MDL_MASK,fnMask);
            Imask.read(fnMask);
            Imask().setXmippOrigin();
            dummy.fn_mask=Imask.name();
            selfTranslate(LINEAR,Imask(),vectorR2(dummy.x,dummy.y),DONT_WRAP);
            Imask().binarize(0.5);
            typeCast(Imask(),mask);
        }
        dummy.corr=correlation_index(theo(),I(),&mask);
        list_of_assigned.push_back(dummy);
    }
}

// Predict shift and psi -----------------------------------------------------
//#define DEBUG
void Prog_angular_predict_tomography_prm::predict_angles(int i)
{
    AlignmentTomography newAlignment=list_of_assigned[i];
    Image<double> I, Imask;
    I.read(list_of_assigned[i].fn_img);
    if (list_of_assigned[i].fn_mask!="")
        Imask.read(list_of_assigned[i].fn_mask);
    MultidimArray<int> mask;
    typeCast(Imask(),mask);

    double rot0=list_of_assigned[i].rot-max_rot_change;
    double rotF=list_of_assigned[i].rot+max_rot_change;
    double tilt0=list_of_assigned[i].tilt-max_tilt_change;
    double tiltF=list_of_assigned[i].tilt+max_tilt_change;
    if (i>0)
        tilt0=XMIPP_MAX(tilt0,
                        (list_of_assigned[i].tilt+list_of_assigned[i-1].tilt)/2);
    if (i<list_of_assigned.size()-1)
        tiltF=XMIPP_MIN(tiltF,
                        (list_of_assigned[i].tilt+list_of_assigned[i+1].tilt)/2);

    for (double rot = rot0; rot <= rotF; rot += rot_step)
        for (double tilt = tilt0; tilt <= tiltF; tilt += tilt_step)
        {
            // Take a projection from the given direction
            Projection theo;
            project_Volume(V(), theo, YSIZE(V()), XSIZE(V()), rot, tilt, 0);
            double theo_avg, theo_stddev, min_val, max_val;
            computeStats_within_binary_mask(mask, theo(),
                                            theo_avg, theo_stddev, min_val, max_val);
            theo() -= theo_avg;

            // Compare it to all possible rotations and shifts
            // of the experimental image
            Image<double> Ip, Imaskp;
            Matrix1D<double> shift(2);
            double xC=list_of_assigned[i].x;
            double yC=list_of_assigned[i].y;
            double x0=xC-max_shift_change;
            double xF=xC+max_shift_change;
            double y0=yC-max_shift_change;
            double yF=yC+max_shift_change;
            for (double x = x0; x <= xF; x += shift_step)
                for (double y = y0; y <= yF; y += shift_step)
                {
                    if ((x-xC)*(x-xC) + (y-yC)*(y-yC) >
                        max_shift_change*max_shift_change)
                        continue;
                    double psi0=list_of_assigned[i].psi-max_psi_change;
                    double psiF=list_of_assigned[i].psi+max_psi_change;
                    for (double psi = psi0; psi <= psiF; psi += psi_step)
                    {
                        // Shift image if necessary
                        if (x == 0 && y == 0)
                        {
                            Ip() = I();
                            Imaskp() = Imask();
                        }
                        else
                        {
                            VECTOR_R2(shift, x, y);
                            translate(LINEAR, Ip(), I(), shift, DONT_WRAP);
                            translate(LINEAR, Imaskp(), Imask(), shift, DONT_WRAP);
                        }

                        // Rotate image if necessary
                        // Adding 2 is a trick to avoid that the 0, 90, 180 and 270
                        // are treated in a different way
                        if (psi != 0)
                        {
                            selfRotate(LINEAR,Ip(),psi + 2, DONT_WRAP);
                            selfRotate(LINEAR,Ip(),-2, DONT_WRAP);
                            selfRotate(LINEAR,Imaskp(),psi, DONT_WRAP);
                        }
                        Imaskp().binarize(0.5);
                        typeCast(Imaskp(),mask);

                        // Compute the correlation index
                        double read_avg, read_stddev;
                        computeStats_within_binary_mask(mask, Ip(),
                                                        read_avg, read_stddev, min_val, max_val);
                        double correlation_index = 0;
                        double N=0;
                        FOR_ALL_ELEMENTS_IN_ARRAY2D(Ip())
                        if (Imaskp(i,j))
                        {
                            correlation_index += (Ip(i, j) - read_avg) * theo(i, j);
                            N++;
                        }
                        correlation_index /= N;
                        correlation_index /= read_stddev * theo_stddev;

                        // Keep the value
                        if (correlation_index>newAlignment.corr)
                        {
                            newAlignment.rot = rot;
                            newAlignment.tilt = tilt;
                            newAlignment.psi = psi;
                            newAlignment.x = x;
                            newAlignment.y = y;
                            newAlignment.corr = correlation_index;
#ifdef DEBUG

                            ImageXmipp save;
                            save() = theo() - theo_avg;
                            save.write("PPPtheo.xmp");
                            save() = Ip() - read_avg;
                            save.write("PPPexp.xmp");
                            multiplyElements((theo() - theo_avg), (Ip() - read_avg), save());
                            save.write("PPPprod.xmp");
                            char c;
                            std::cin >> c;
#endif

                        }
                    }
                }
        }

    // Select best alignment
    list_of_assigned[i]=newAlignment;
}
#undef DEBUG

// Finish processing ---------------------------------------------------------
void Prog_angular_predict_tomography_prm::run()
{
    MetaData DF;
    for (int i=0; i<list_of_assigned.size(); i++)
    {
        // Read input image
        Image<double> I;
        I.read(list_of_assigned[i].fn_img);
        I().setXmippOrigin();

        // Get angles and shifts
        predict_angles(i);

        // Store them in the output docfile
        DF.addObject();
        DF.setValue(MDL_IMAGE,fn_out+"_img_"+integerToString(i,4)+".xmp");
        DF.setValue(MDL_IMAGE_ORIGINAL,list_of_assigned[i].fn_img);
        DF.setValue(MDL_MASK,fn_out+"_mask_"+integerToString(i,4)+".xmp");
        DF.setValue(MDL_ANGLEROT,list_of_assigned[i].rot);
        DF.setValue(MDL_ANGLETILT,list_of_assigned[i].tilt);
        DF.setValue(MDL_ANGLEPSI,list_of_assigned[i].psi);
        DF.setValue(MDL_SHIFTX,list_of_assigned[i].x);
        DF.setValue(MDL_SHIFTY,list_of_assigned[i].y);
        DF.setValue(MDL_MAXCC,list_of_assigned[i].corr);

        // Shift the image and mask
        selfTranslate(LINEAR,I(),
                      vectorR2(list_of_assigned[i].x,list_of_assigned[i].y),DONT_WRAP);
        Image<double> Imask;
        MultidimArray<int> mask;
        const MultidimArray<int> *maskPtr=NULL;
        if (list_of_assigned[i].fn_mask!="")
        {
            Imask.read(list_of_assigned[i].fn_mask);
            Imask().setXmippOrigin();
            selfTranslate(LINEAR,Imask(),
                          vectorR2(list_of_assigned[i].x,list_of_assigned[i].y),DONT_WRAP);
            Imask().binarize(0.5);
            typeCast(Imask(),mask);
            maskPtr=&mask;
        }

        // Adjust Gray values
        if (adjustGray)
        {
            Projection theo;
            project_Volume(V(), theo, YSIZE(V()), XSIZE(V()),
                           list_of_assigned[i].rot,
                           list_of_assigned[i].tilt,
                           list_of_assigned[i].psi);
            I().rangeAdjust(theo(),maskPtr);
        }

        // Set angles and offsets
        I.setEulerAngles(list_of_assigned[i].rot,
                         list_of_assigned[i].tilt,
                         list_of_assigned[i].psi);
        I.setShifts(0,0);
        Imask.setEulerAngles(list_of_assigned[i].rot,
                             list_of_assigned[i].tilt,
                             list_of_assigned[i].psi);
        Imask.setShifts(0,0);

        I.write(fn_out+"_img_"+integerToString(i,4)+".xmp");
        I.write(fn_out+"_mask_"+integerToString(i,4)+".xmp");
    }
    DF.write(fn_out+".doc");
}
