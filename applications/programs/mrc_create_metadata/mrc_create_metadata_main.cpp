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

class ProgMrcCreateMetaData: public XmippProgram
{
public:
    FileName fnStack, fnAngles, fnOut;

    void defineParams()
    {
    	addUsageLine("Produces a selfile valid for tomography with an MRC stack and a list of angles");
        addParamsLine(" --stack <stack>       : Stack file (not necessarily MRC)");
        addParamsLine(" --tiltAngles <rawtlt> : Tilt angle file");
        addParamsLine(" -o <metadata>         : Output metadata");
    }

    void readParams()
    {
    	fnStack = getParam("--stack");
    	fnAngles = getParam("--tiltAngles");
    	fnOut = getParam("-o");
    }

    void show()
    {
        std::cout
        << "Stack   = " << fnStack  << std::endl
        << "Angles  = " << fnAngles << std::endl
        << "Output  = " << fnOut << std::endl
        ;
    }

    void run()
    {
    	std::string extension=fnStack.getExtension();
    	if (extension=="mrc")
    		fnStack+=":mrcs";
    	else if (extension=="ali")
    		fnStack+=":mrcs";

    	std::ifstream fhAngles;
    	fhAngles.open(fnAngles.c_str());
    	if (!fhAngles)
    		REPORT_ERROR(ERR_IO_NOTOPEN,fnAngles);

    	Image<double> stack;
    	stack.read(fnStack, HEADER);
    	size_t xdim, ydim, zdim, ndim;
    	stack.getDimensions(xdim, ydim, zdim, ndim);

    	MetaData MD;
    	FileName fnSlice;
    	for (unsigned long slice = FIRST_IMAGE; slice <= ndim; slice++)
    	{
    		double angle;
    		fhAngles >> angle;
    		fnSlice.compose((int)slice,fnStack);
    		size_t objId=MD.addObject();
    		MD.setValue(MDL_IMAGE,fnSlice,objId);
    		MD.setValue(MDL_ANGLE_TILT,angle,objId);
    		MD.setValue(MDL_ENABLED,1,objId);
    	}
    	MD.write(fnOut);
    }
};

int main(int argc, char *argv[])
{
    	ProgMrcCreateMetaData program;
        program.read(argc, argv);
        program.tryRun();
}
