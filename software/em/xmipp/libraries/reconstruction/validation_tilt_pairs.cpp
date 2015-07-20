/***************************************************************************
 *
 * Authors:    Jose Luis Vilas          (jlvilas@cnb.csic.es)
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
#include "validation_tilt_pairs.h"



//Define Program parameters
void ProgValidationTiltPairs::defineParams()
{
    //Usage
    addUsageLine("Takes two coordinates sets and defines the coordinate transformation between them");
	addUsageLine("First set defines the untilted coordinates, second set defines the tilted coordinates");
	addParamsLine(" --untilted <metadata> : Metadata with untilted images");
	addParamsLine(" --tilted <metadata> : Metadata with tilted images");
}

//Read params
void ProgValidationTiltPairs::readParams()
{
    fntilted_In = getParam("--tilted");  //Set of tilted coordinates
	fnuntilted_In = getParam("--untilted");  //Set of untilted coordinates

	fnOut = getParam("-o");  //Output file

	MetaData MD_tilted, MD_untilted;
	String errorMessage;
	double Xoff;
	FOR_ALL_OBJECTS_IN_METADATA(fntilted_In) //loop through all lines
	{
	    if (MD_tilted.getValue(MDL_X, Xoff,__iter.objId))//get value for attribute X
	        std::cerr << "The value is: " <<  Xoff;
	    else{
	        errorMessage = formatString("Cannot find label %s ",(MDL::label2Str(MDL_X)).c_str());
	        REPORT_ERROR(ERR_MD_MISSINGLABEL,errorMessage);
	    }
	    MD_tilted.setValue(MDL_X,Xoff*2.);//store the double of the sampling rate
	}
	std::cout << MD_tilted <<std::endl;

}


