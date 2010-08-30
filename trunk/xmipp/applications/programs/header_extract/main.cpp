/***************************************************************************
 *
 * Authors:    Sjors Scheres
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
#include <data/image.h>
#include <data/metadata.h>

void Usage();

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    bool            round_shifts = false;
    double           xx, yy;
    FileName        fn_img, fn_out,fn_in;
    Image<double>      img;
    FileName        inputFile;

// Check command line options ===========================================
    try
    {
        
        fn_in  = getParameter(argc, argv, "-i"); 
	    if (checkParameter(argc, argv, "-o"))
            fn_out = getParameter(argc, argv, "-o");
	    else
	        fn_out = fn_in;   
            round_shifts = checkParameter(argc, argv, "-round_shifts");
    }
    catch (XmippError XE)
    {
        std::cout << XE;
        Usage();
    }
    MetaData SF(fn_in);
    SF.removeObjects(MDValueEQ(MDL_ENABLED, -1));
    // Extracting information  ==================================================
    try
    {

	    if(SF.isEmpty())
	    {
	        std::cerr << "Empty inputFile File\n";
	        exit(1);
	    }

	    FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
	        SF.getValue(MDL_IMAGE,fn_img);
            if (fn_img=="") break;
	    img.read(fn_img,false);
	    xx = (double) img.Xoff();
	    yy = (double) img.Yoff();
            if (round_shifts)
            {
                xx = (double)ROUND(xx);
                yy = (double)ROUND(yy);
            }
	        SF.setValue(MDL_ANGLEROT,  (double) img.rot());
	        SF.setValue(MDL_ANGLETILT, (double) img.tilt());
    	    SF.setValue(MDL_ANGLEPSI,  (double) img.psi());
	        SF.setValue(MDL_SHIFTX,    xx );
	        SF.setValue(MDL_SHIFTY,    yy );
	        SF.setValue(MDL_WEIGHT,    (double) img.weight());
	        SF.setValue(MDL_FLIP,      (bool)   img.flip());
        }
	
        SF.write(fn_out);
    }
    catch (XmippError XE)
    {
        std::cout << XE;
    }

}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cout << " Purpose:\n";
    std::cout << " Extracts the geometric transformation (angles & shifts) in the header of 2D-images.\n";
    std::cout << " Usage:\n";
    std::cout << "    header_extract \n";
    std::cout <<              "        -i <metadata>       : metadata selfile\n";
    std::cout << (std::string)"       [-o <docfile> ]     : metaData file, by default data\n" +
                              "                             is stored in input metaData file\n";
    std::cout <<              "       [-round_shifts]     : Round shifts to integers \n";
    exit(1);
}
