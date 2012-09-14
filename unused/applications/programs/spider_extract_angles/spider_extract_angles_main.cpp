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
#include <data/program.h>

class ProgExtractAngles: public XmippMetadataProgram
{
protected:
    bool from_metadata;
    std::string       ang1,ang2,ang3;
    DocFile DF_out;

    void defineParams()
    {
        produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addParamsLine("[-order <ang1=psi> <ang2=rot> <ang3=tilt>]   : order of the angles");
        addParamsLine("[-from_metadata]         : Read angles from MetaData");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        from_metadata = checkParam("-from_metadata");
        ang1 = getParam("-order", 0);
        ang2 = getParam("-order", 1);
        ang3 = getParam("-order", 2);
    }

    void show()
    {
        if (!verbose)
            return;

        std::cout
        << "Input Selfile = " << fn_in << std::endl
        << "Angle 1 = " << ang1 << std::endl
        << "Angle 2 = " << ang2 << std::endl
        << "Angle 3 = " << ang3 << std::endl
        << "Output Selfile = " << fn_out << std::endl;
    }

    void preProcess()
    {
        checkAngle(ang1);
        checkAngle(ang2);
        checkAngle(ang3);

        DF_out.append_comment((std::string)"Angles for " + fn_in +
                              ".   Angle order: " + ang1 + " " + ang2 + " " + ang3);

        if (verbose)
            std::cerr << "Extracting angles ...\n";
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        if (from_metadata)
        {
            double rot, tilt, psi;
            mdIn.getValue(MDL_ANGLE_ROT,rot, objId);
            mdIn.getValue(MDL_ANGLE_TILT,tilt, objId);
            mdIn.getValue(MDL_ANGLE_PSI,psi, objId);
            DF_out.append_angles(rot, tilt, psi, ang1, ang2, ang3);
        }
        else
        {
            Image<double> img;
            img.readApplyGeo(fnImg,mdIn,objId, HEADER);
            DF_out.append_angles(img.rot(), img.tilt(), img.psi(),
                                 ang1, ang2, ang3);
        }
    }

    void postProcess()
    {
        DF_out.write(fn_out);
    }

}; ///end of class ProgSpiderExtract

int main(int argc, char *argv[])
{

    try
    {
        ProgExtractAngles program;
        program.read(argc, argv);
        program.run();

    }
    catch(XmippError xe)
    {
        std::cerr << xe;
    }
}
