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

class ProgGenerateCount: public XmippProgram
{
protected:
    int maxcount;
    FileName fn_out;

    void defineParams()
    {
        addParamsLine("-o <selfile>   : Output Spider SelFile");
        addParamsLine("-max <max_count>          : Generate corresponding Spider SelFile");
    }

    void readParams()
    {
        maxcount = getIntParam("-max");
        fn_out = getParam("-o");
    }

    void show()
    {
        if (!verbose)
            return;

        std::cout
        << "Max count = " << maxcount << std::endl
        << "Output Selfile = " << fn_out << std::endl;

    }


public:
    void run()
    {
        Matrix1D<double>   aux(1);

        DocFile DF_out;
        DF_out.append_comment((std::string)"Count for Spider up to " + integerToString(maxcount));

        for (aux(0) = 1; aux(0) <= maxcount; ++aux(0))
            DF_out.append_data_line(aux);

        DF_out.write(fn_out);
    }
};///end of class ProgGenerateCount

int main(int argc, char *argv[])
{

    try
    {
        ProgGenerateCount program;
        program.read(argc, argv);
        program.run();

    }
    catch(XmippError xe)
    {
        std::cerr << xe;
    }
}

