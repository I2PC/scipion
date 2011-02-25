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
void Prog_angular_predict_tomography_prm::readParams()
{
    fn_ref = getParam("--ref");
    fn_sel = getParam("--sel");
    fn_out = getParam("--oroot");
    max_rot_change = getDoubleParam("--max_rot_change");
    max_tilt_change = getDoubleParam("--max_tilt_change");
    max_psi_change = getDoubleParam("--max_psi_change");
    rot_step = getDoubleParam("--rot_step");
    tilt_step = getDoubleParam("--tilt_step");
    psi_step = getDoubleParam("--psi_step");
    max_shift_change = getDoubleParam("--max_shift_change");
    shift_step = getDoubleParam("--shift_step");
    adjustGray = checkParam("--adjustGray");
    produce_side_info();
}

// Show ====================================================================
void Prog_angular_predict_tomography_prm::show()
{
    std::cout
    << "Reference images:   " << fn_ref           << std::endl
    << "Input images:       " << fn_sel           << std::endl
    << "Ouput rootname:     " << fn_out           << std::endl
    << "Max rot change:     " << max_rot_change   << " step: " << rot_step << std::endl
    << "Max tilt change:    " << max_tilt_change  << " step: " << tilt_step << std::endl
    << "Max psi change:     " << max_psi_change   << " step: " << psi_step << std::endl
    << "Max shift change:   " << max_shift_change << " step: " << shift_step << std::endl
    ;
}

// usage ===================================================================
void Prog_angular_predict_tomography_prm::defineParams()
{
    addUsageLine("Realign a tilt series");
    addParamsLine("   --ref <volume>             : Reference volume");
    addParamsLine("   --sel <selfile>            : Images to align");
    addParamsLine("   --oroot <rootname>         : rootname for the output");
    addParamsLine("  [--max_rot_change <ang=5>]  : Maximum change allowed in rot");
    addParamsLine("  [--max_tilt_change <ang=2>] : Maximum change allowed in tilt");
    addParamsLine("  [--max_psi_change <ang=5>]  : Maximum change allowed in psi");
    addParamsLine("  [--max_shift_change <r=10>] : Maximum change allowed in shift");
    addParamsLine("  [--rot_step <ang=1>]        : Rot search step");
    addParamsLine("  [--tilt_step <ang=1>]       : Tilt search step");
    addParamsLine("  [--psi_step <ang=3>]        : Psi search step");
    addParamsLine("  [--shift_step <r=2>]        : Step in shift in pixels");
    addParamsLine("  [--adjustGray]              : Adjust also gray values");
}

// Produce side information ================================================
void Prog_angular_predict_tomography_prm::produce_side_info()
{
    V.read(fn_ref);
    V().setXmippOrigin();
    MetaData SF(fn_sel);
    FileName fnAux;
    bool masksPresent = SF.getValue(MDL_MASK, fnAux, SF.firstObject());
    Image<double> I, Imask;
    Projection theo;
    MultidimArray<int> mask;
    FileName fnImg;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        AlignmentTomography dummy;
        SF.getValue(MDL_IMAGE,fnImg,__iter.objId);
        I.readApplyGeo(fnImg,SF,__iter.objId);
        I().setXmippOrigin();
        dummy.rot=I.rot();
        dummy.tilt=I.tilt();
        dummy.psi=I.psi();
        dummy.x=I.Xoff();
        dummy.y=I.Yoff();
        dummy.fn_img=fnImg;
        projectVolume(V(), theo, YSIZE(V()), XSIZE(V()),
                      dummy.rot, dummy.tilt, dummy.psi);
        selfTranslate(LINEAR,I(),vectorR2(dummy.x,dummy.y),DONT_WRAP);
        if (masksPresent)
        {
            FileName fnMask;
            SF.getValue(MDL_MASK,fnMask,__iter.objId);
            Imask.read(fnMask);
            Imask().setXmippOrigin();
            dummy.fn_mask=fnMask;
            selfTranslate(LINEAR,Imask(),vectorR2(dummy.x,dummy.y),DONT_WRAP);
            Imask().binarize(0.5);
            typeCast(Imask(),mask);
            dummy.corr=correlation_index(theo(),I(),&mask);
        }
        else
            dummy.corr=correlation_index(theo(),I());
        list_of_assigned.push_back(dummy);
    }
}

// Lookf for best angles -------------------------------------------------
//#define DEBUG
void Prog_angular_predict_tomography_prm::predict_angles(int i)
{
    AlignmentTomography newAlignment=list_of_assigned[i];
    Image<double> I, Imask;
    I.read(list_of_assigned[i].fn_img);
    bool masksPresent=list_of_assigned[i].fn_mask!="";
    if (masksPresent)
        Imask.read(list_of_assigned[i].fn_mask);

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

    Projection theo;
    Matrix2D<double> M;
    MultidimArray<int> mask;
    for (double rot = rot0; rot <= rotF; rot += rot_step)
        for (double tilt = tilt0; tilt <= tiltF; tilt += tilt_step)
        {
            // Take a projection from the given direction
            projectVolume(V(), theo, YSIZE(V()), XSIZE(V()), rot, tilt, 0);
            alignImages(theo(),I(),M);

            // Compute the correlation index
            double newCorr;
            if (masksPresent)
            {
                Imask.read(list_of_assigned[i].fn_mask);
                Imask().setXmippOrigin();
                selfApplyGeometry(LINEAR,Imask(),M,IS_NOT_INV,DONT_WRAP);
                Imask().binarize(0.5);
                typeCast(Imask(),mask);
                newCorr=correlation_index(theo(),I(),&mask);
            }
            else
                newCorr=correlation_index(theo(),I());

            // Keep the value
            if (newCorr>newAlignment.corr)
            {
                newAlignment.rot = rot;
                newAlignment.tilt = tilt;
                bool flip;
                double scale;
                newAlignment.M=M;
                transformationMatrix2Parameters(M, flip, scale, newAlignment.x, newAlignment.y,
                                                newAlignment.psi);
                newAlignment.corr = newCorr;
#ifdef DEBUG
                ImageXmipp save;
                save() = theo();
                save.write("PPPtheo.xmp");
                save() = Ip();
                save.write("PPPexp.xmp");
                save()=theo()-Ip();
                save.write("PPPdiff.xmp");
                char c;
                std::cin >> c;
#endif
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
    FileName fnMaskOut, fnImgOut;
    Projection theo;
    Image<double> Imask, I;

    if (adjustGray && exists(fn_out+".stk"))
    	unlink((fn_out+".stk").c_str());

    for (int i=0; i<list_of_assigned.size(); i++)
    {
        // Get angles and shifts
        predict_angles(i);

        // Store them in the output docfile
        size_t id = DF.addObject();
        if (adjustGray)
        	fnImgOut.compose(i+1,fn_out+".stk");
        else
        	fnImgOut=list_of_assigned[i].fn_img;
        DF.setValue(MDL_IMAGE,fnImgOut,id);
        DF.setValue(MDL_IMAGE_ORIGINAL,list_of_assigned[i].fn_img,id);
        bool masksPresent=list_of_assigned[i].fn_mask!="";
        if (masksPresent)
        {
            fnMaskOut.compose(i+1,fn_out+"_mask.stk");
            DF.setValue(MDL_MASK,fnMaskOut,id);
        }
        DF.setValue(MDL_ANGLEROT,list_of_assigned[i].rot,id);
        DF.setValue(MDL_ANGLETILT,list_of_assigned[i].tilt,id);
        DF.setValue(MDL_ANGLEPSI,list_of_assigned[i].psi,id);
        DF.setValue(MDL_SHIFTX,list_of_assigned[i].x,id);
        DF.setValue(MDL_SHIFTY,list_of_assigned[i].y,id);
        DF.setValue(MDL_MAXCC,list_of_assigned[i].corr,id);

        // Adjust Gray values
        if (adjustGray)
        {
        	I.read(list_of_assigned[i].fn_img);
        	I().setXmippOrigin();
            projectVolume(V(), theo, YSIZE(V()), XSIZE(V()),
                          list_of_assigned[i].rot,
                          list_of_assigned[i].tilt,
                          list_of_assigned[i].psi);
            MultidimArray<int> mask;
            if (masksPresent)
            {
                Imask.read(list_of_assigned[i].fn_mask);
                Imask().setXmippOrigin();
                selfApplyGeometry(LINEAR,Imask(),list_of_assigned[i].M,IS_NOT_INV,DONT_WRAP);
                Imask().binarize(0.5);
                typeCast(Imask(),mask);
            }
            I().rangeAdjust(theo(),&mask);
            I.write(fnImgOut);
        }
    }
    DF.write(fn_out+".doc");
}
