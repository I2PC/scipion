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

class ProgSpiderRename: public XmippMetadataProgram
{
protected:

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        if (!checkParam("-o"))
            REPORT_ERROR(ERR_ARG_MISSING, "-o is required");
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        static int counter = 0;
        ++counter;
        FileName fnTemp(fnImgOut);
        fnTemp.compose(oroot, counter, oext);
        if (verbose)
            std::cout << "Renaming " << fn_in << " as " << fn_out << std::endl;
        std::string command = (std::string)"cp " + fn_in + " " + fn_out;
        system(command.c_str());
    }

}
; //end of class ProgSpiderRename

int main(int argc, char *argv[])
{
    ProgSpiderRename program;
    program.read(argc, argv);
    program.tryRun();
}
