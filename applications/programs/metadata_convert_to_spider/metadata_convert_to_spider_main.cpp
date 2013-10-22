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
#include <data/args.h>
#include <interface/spider.h>
#include <data/xmipp_program.h>

class ProgSpiderTranslate: public XmippMetadataProgram
{
public:
    DocFile DF_out;
    String action,ang1,ang2,ang3;
    bool new_style;
    int currentImage;
    ApplyGeoParams params;
    bool ReadImg;

    void defineParams()
    {
        produces_an_output = true;
        addUsageLine("Extracts information from a Xmipp selfile into a Spider docfile");
        XmippMetadataProgram::defineParams();
        addParamsLine("--action <action>                       : Choose one of the following choices");
        addParamsLine("   where <action>");
        addParamsLine("         extract_selfile <style=new>    : Valid styles are new or old");
        addParamsLine("                                        : The old style is 0=not selected, 1=selected");
        addParamsLine("                                        : The new style contains the numbers");
        addParamsLine("                                        : corresponding to the selected images");
        addParamsLine("         extract_angles <ang1=psi> <ang2=rot> <ang3=tilt> : Specify the order");
        addParamsLine("                                        : in which the angles should be written");
        addParamsLine("         generate_count                 : Generate a count file with as many");
        addParamsLine("                                        : numbers as enabled files in the metadata");
        addParamsLine("[--disregard_disabled]                  : Disregard disabled images from the metadata");
        addParamsLine("[--do_not_read_img]                     : Ignore image label");
    }

    void readParams()
    {
        remove_disabled=checkParam("--disregard_disabled");
        XmippMetadataProgram::readParams();
        action=getParam("--action");
        ReadImg=!checkParam("--do_not_read_img");
        if (action=="extract_selfile")
        {
            String style=getParam("--action",1);
            new_style = style=="new";
        }
        else if (action=="extract_angles")
        {
            ang1=getParam("--action",1);
            ang2=getParam("--action",2);
            ang3=getParam("--action",3);
        }
        params.datamode = HEADER;
    }

    void show()
    {
        if (!verbose)
            return;

        std::cout
        << "Input Selfile  = " << fn_in << std::endl
        << "Output Selfile = " << fn_out << std::endl
        << "Action         = " << action << std::endl;
        if (action=="extract_selfile")
            std::cout << "New style      = " << new_style << std::endl;
        else if (action=="extract_angles")
            std::cout << "Angle order    = " << ang1 << " " << ang2 << " " << ang3 << std::endl;
    }

    void preProcess()
    {
        DF_out.clear();
        DF_out.append_comment((std::string)"Translation for Spider of " + fn_in);
        currentImage=1;
        if (action=="extract_angles")
        {
            checkAngle(ang1);
            checkAngle(ang2);
            checkAngle(ang3);
            DF_out.append_comment((std::string)"Angle order: " + ang1 + " " + ang2 + " " + ang3);
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        if (action=="extract_selfile")
        {
            bool store = true;
            int enabled;
            Matrix1D<double> aux(1); // Auxiliary vector to be added to the docfile
            if (!rowIn.getValue( MDL_ENABLED, enabled))
                enabled=1;
            if (enabled == 1)
                aux(0) = (new_style) ? currentImage : 1;
            else
            {
                if (new_style)
                    store = false;
                else
                    aux(0) = 0;
            }
            if (store)
                DF_out.append_data_line(aux);
        }
        else if (action=="extract_angles")
        {
        	double rot,tilt,psi;
            if (ReadImg)
            {
                Image<double> img;
                img.readApplyGeo(fnImg,rowIn, params);
                rot = img.rot();
                tilt =img.tilt();
                psi =img.psi();
            }
            else
            {
            	rowIn.getValue(MDL_ANGLE_ROT, rot);
            	rowIn.getValue(MDL_ANGLE_TILT, tilt);
            	rowIn.getValue(MDL_ANGLE_PSI, psi);
            }

            DF_out.append_angles(rot,tilt,psi,
                                 ang1, ang2, ang3);

        }
        else if (action=="generate_count")
        {
            Matrix1D<double> aux(1); // Auxiliary vector to be added to the docfile
            aux(0) = currentImage;
            DF_out.append_data_line(aux);
        }
        currentImage++;
    }

    void postProcess()
    {
        DF_out.write(fn_out);
    }

}
;//end of class ProgSpiderTranslate

RUN_XMIPP_PROGRAM(ProgSpiderTranslate)
