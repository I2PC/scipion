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

#include "tomo_align_refinement.h"
#include <data/projection.h>

#include <data/args.h>
#include <data/filters.h>
#include <data/xmipp_image_base.h>

// Empty constructor =======================================================
AlignmentTomography::AlignmentTomography()
{
    rot=tilt=psi=x=y=corr=0;
    fn_img="";
}

// Read arguments ==========================================================
void ProgTomoAlignRefinement::readParams()
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
    generateAligned = checkParam("--generateAligned");
}

// Show ====================================================================
void ProgTomoAlignRefinement::show()
{
    std::cout
    << "Reference images:   " << fn_ref           << std::endl
    << "Input images:       " << fn_sel           << std::endl
    << "Ouput rootname:     " << fn_out           << std::endl
    << "Max rot change:     " << max_rot_change   << " step: " << rot_step << std::endl
    << "Max tilt change:    " << max_tilt_change  << " step: " << tilt_step << std::endl
    << "Max psi change:     " << max_psi_change   << " step: " << psi_step << std::endl
    << "Max shift change:   " << max_shift_change << " step: " << shift_step << std::endl
    << "Adjust gray:        " << adjustGray       << std::endl
    << "Generate aligned:   " << generateAligned  << std::endl
    ;
}

// usage ===================================================================
void ProgTomoAlignRefinement::defineParams()
{
    addUsageLine("Realign a tilt series with respect to a volume following a single ");
    addUsageLine("particles approach, i.e., the volume is reprojected and the alignment");
    addUsageLine("parameters are reoptimized.");
    addParamsLine("   --ref <volume>             : Reference volume");
    addParamsLine("   --sel <selfile>            : Images to align");
    addParamsLine("   --oroot <rootname>         : rootname for the output");
    addParamsLine("                              : rootname.doc contains a selfile with the");
    addParamsLine("                              : images realigned");
    addParamsLine("                              : rootname.stk contains the aligned and");
    addParamsLine("                              : adjusted images in case --adjustGray or");
    addParamsLine("                              : --generateAligned are given");
    addParamsLine("  [--max_rot_change <ang=2>]  : Maximum change allowed in rot");
    addParamsLine("  [--max_tilt_change <ang=2>] : Maximum change allowed in tilt");
    addParamsLine("  [--max_psi_change <ang=2>]  : Maximum change allowed in psi");
    addParamsLine("  [--max_shift_change <r=10>] : Maximum change allowed in shift");
    addParamsLine("  [--rot_step <ang=0.25>]     : Rot search step");
    addParamsLine("  [--tilt_step <ang=0.25>]    : Tilt search step");
    addParamsLine("  [--psi_step <ang=0.25>]     : Psi search step");
    addParamsLine("  [--shift_step <r=2>]        : Step in shift in pixels");
    addParamsLine("  [--adjustGray]              : Adjust also gray values");
    addParamsLine("  [--generateAligned]         : Generate aligned images");
}

// Produce side information ================================================
//#define DEBUG
void ProgTomoAlignRefinement::produce_side_info()
{
    V.read(fn_ref);
    V().setXmippOrigin();
    MetaData SF(fn_sel);
    Image<double> I;
    FileName fnImg;
    ApplyGeoParams params;
    params.datamode = HEADER;

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        AlignmentTomography dummy;
        SF.getValue(MDL_IMAGE,fnImg,__iter.objId);
        I.readApplyGeo(fnImg,SF,__iter.objId, params);
        I().setXmippOrigin();
        dummy.rot=I.rot();
        dummy.tilt=I.tilt();
        dummy.psi=I.psi();
        dummy.x=I.Xoff();
        dummy.y=I.Yoff();
        dummy.fn_img=fnImg;
        dummy.corr=-1;
        list_of_assigned.push_back(dummy);
    }
}
#undef DEBUG

// Look for best angles -------------------------------------------------
//#define DEBUG
void ProgTomoAlignRefinement::predict_angles(size_t idx,
        const FileName &fnImgOut)
{
    AlignmentTomography newAlignment=list_of_assigned[idx];
    Image<double> I, Ip;
    MultidimArray<int> mask;
    I.read(list_of_assigned[idx].fn_img);
    I().setXmippOrigin();
    mask.resizeNoCopy(I());

    double rot0=list_of_assigned[idx].rot-max_rot_change;
    double rotF=list_of_assigned[idx].rot+max_rot_change;
    double tilt0=list_of_assigned[idx].tilt-max_tilt_change;
    double tiltF=list_of_assigned[idx].tilt+max_tilt_change;
    if (idx>0)
        tilt0=XMIPP_MAX(tilt0,
                        (list_of_assigned[idx].tilt+list_of_assigned[idx-1].tilt)/2);
    if (idx<list_of_assigned.size()-1)
        tiltF=XMIPP_MIN(tiltF,
                        (list_of_assigned[idx].tilt+list_of_assigned[idx+1].tilt)/2);

    Projection theo;
    Matrix2D<double> M;
#ifdef DEBUG

    std::cout << "original idx,rot,tilt,psi,x,y=" << idx << " "
    << newAlignment.rot << " " << newAlignment.tilt
    << " " << newAlignment.psi << " "
    << newAlignment.x << " " << newAlignment.y << std::endl;
#endif

    for (double rot = rot0; rot <= rotF; rot += rot_step)
        for (double tilt = tilt0; tilt <= tiltF; tilt += tilt_step)
        {
            Ip()=I();
            // Take a projection from the given direction
            projectVolume(V(), theo, YSIZE(V()), XSIZE(V()), rot, tilt, 0);
            mask.initConstant(1);
            FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
            if (IMGPIXEL(theo,i,j)==0 || IMGPIXEL(Ip,i,j)==0)
            {
                A2D_ELEM(mask,i,j)=0;
                IMGPIXEL(Ip,i,j)=0;
            }
            double newCorr=correlationIndex(theo(),Ip(),&mask);
            if (newCorr>newAlignment.corr)
            {
                newAlignment.rot = rot;
                newAlignment.tilt = tilt;
                newAlignment.corr=newCorr;
                if (adjustGray)
                    Ip().rangeAdjust(theo(),&mask);
                if (adjustGray || generateAligned)
                    Ip.write(fnImgOut);
#ifdef DEBUG

                std::cout << "    Improved " << idx << " rot=" << rot << " tilt=" << tilt
                << " x,y="
                << newAlignment.x << " " << newAlignment.y << " psi="
                << newAlignment.psi << " corr="
                << newCorr << std::endl;
                Image<double> save;
                save()=Ip();
                save.write("PPPexp.xmp");
                save()=theo();
                save.write("PPPtheo.xmp");
                typeCast(mask,save());
                save.write("PPPmask.xmp");
                save()=theo()-Ip();
                save.write("PPPdiff.xmp");
#endif

            }

            // Look for better alignment
            alignImages(theo(),Ip(),M, DONT_WRAP);

            // Measure the new correlation
            newCorr=correlationIndex(theo(),Ip(),&mask);

            // Keep the value
            if (newCorr>newAlignment.corr)
            {
                newAlignment.rot = rot;
                newAlignment.tilt = tilt;
                bool flip;
                double scale;
                newAlignment.M=M;
                transformationMatrix2Parameters2D(M, flip, scale, newAlignment.x, newAlignment.y,
                                                newAlignment.psi);
                if (adjustGray)
                    Ip().rangeAdjust(theo(),&mask);
                if (adjustGray || generateAligned)
                    Ip.write(fnImgOut);
#ifdef DEBUG

                std::cout << "    Improved " << idx << " rot=" << rot << " tilt=" << tilt
                << " x,y="
                << newAlignment.x << " " << newAlignment.y << " psi="
                << newAlignment.psi << " corr="
                << newCorr << std::endl;
                newAlignment.corr = newCorr;
                Image<double> save;
                save() = Ip();
                save.write("PPPexpRealigned.xmp");
                save()=theo()-Ip();
                save.write("PPPdiffRealigned.xmp");
                std::cout << "Press any key\n";
                char c;
                std::cin >> c;
#endif

            }
        }

    // Select best alignment
    list_of_assigned[idx]=newAlignment;
}
#undef DEBUG

// Finish processing ---------------------------------------------------------
void ProgTomoAlignRefinement::run()
{
    produce_side_info();
    MetaData DF;
    FileName fnMaskOut, fnImgOut;
    Projection theo;
    Image<double> Imask, I;
    FileName fn_tmp(fn_out + ".stk");
    if (adjustGray)
      fn_tmp.deleteFile();

    for (size_t i=0; i<list_of_assigned.size(); i++)
    {
        if (adjustGray || generateAligned)
            fnImgOut.compose(i+1,fn_out+".stk");
        else
            fnImgOut=list_of_assigned[i].fn_img;

        // Get angles and shifts
        predict_angles(i,fnImgOut);

        // Store them in the output docfile
        size_t id = DF.addObject();
        DF.setValue(MDL_IMAGE,fnImgOut,id);
        DF.setValue(MDL_IMAGE_ORIGINAL,list_of_assigned[i].fn_img,id);
        DF.setValue(MDL_ANGLE_ROT,list_of_assigned[i].rot,id);
        DF.setValue(MDL_ANGLE_TILT,list_of_assigned[i].tilt,id);
        DF.setValue(MDL_ANGLE_PSI,list_of_assigned[i].psi,id);
        DF.setValue(MDL_SHIFT_X,list_of_assigned[i].x,id);
        DF.setValue(MDL_SHIFT_Y,list_of_assigned[i].y,id);
        DF.setValue(MDL_MAXCC,list_of_assigned[i].corr,id);
    }
    DF.write(fn_out+".doc");
}
