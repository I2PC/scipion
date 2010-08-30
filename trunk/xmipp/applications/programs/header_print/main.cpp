/***************************************************************************
 *
 * Authors:    Sjors Scheres
 *             Slavica Jonic
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
 * MERCHANTABILITY or FITNESS FO A PARTICULAR PURPOSE.  See the
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
    Image<double>      img;
    FileName        fn_input;
    MetaData SF;

    try
    {
        fn_input = getParameter(argc, argv, "-i");
        if (!fn_input.isMetaData())
        {
            SF.addObject();
            SF.setValue( MDL_IMAGE, fn_input);
            SF.setValue( MDL_ENABLED, 1);
        }
        else {
            SF.read( fn_input ,NULL);
            SF.removeObjects(MDValueEQ(MDL_ENABLED, -1));
        }
    }
    catch (XmippError XE)
    {
        std::cout << XE;
        Usage();
    }

    try
    {
        std::cout << " Printing the header ... " << std::endl;

        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fn_img;
            SF.getValue( MDL_IMAGE, fn_img); 
            if (fn_img=="") break;
            std::cout << "FileName     : " << fn_img << std::endl;

            img.read(fn_img, false, -1, false);

            std::cout << img;
	    std::cout << std::endl;
        }

    }
    catch (XmippError XE)
    {
        std::cout << XE;
    }
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    printf("Purpose:\n");
    printf(" Print information from the header of 2D-images.\n");
    printf("Usage:\n");
    printf("   header_print \n");
    printf("    -i                                   : metaDataFile with images or individual image\n");
    exit(1);
}

