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

class ProgSpiderTranslate: public XmippMetadataProgram
{
protected:
    DocFile DF_out;
    bool new_style;

    void defineParams()
    {
        produces_an_output = true;
        addUsageLine("Translates a Xmipp selfile into a Spider docfile")
        XmippMetadataProgram::defineParams();
        addParamsLine("[--new_style]         : Use Spider new style for docfiles");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        new_style = checkParam("--new_style");
    }

    void show()
    {
        if (!verbose)
            return;

        std::cout
        << "Input Selfile = " << fn_in << std::endl
        << "Spider Selfile style = " << new_style << std::endl
        << "Output Selfile = " << fn_out << std::endl;
    }

    void preProcess()
    {
        DF_out.clear();
        DF_out.append_comment((std::string)"Translation for Spider of " + fn_in);
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        static double aux;
        static int i=1;
        bool store = true;
        int enabled;
        mdIn.getValue( MDL_ENABLED, enabled, objId);

        if (enabled == 1)
        {
            aux = (new_style) ? i++ : 1;
        }
        else
        {
            if (!new_style)
                aux = 0;
            else
            {
                store = false;
                i++;
            }
        }
        if (store)
            DF_out.append_data_line(aux);
    }

    void postProcess()
    {
        DF_out.write(fn_out);
    }

}
;//end of class ProgSpiderTranslate

int main(int argc, char *argv[])
{
    ProgSpiderTranslate program;
    program.read(argc, argv);
    program.tryRun();
}
