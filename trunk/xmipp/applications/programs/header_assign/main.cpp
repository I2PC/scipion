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
    double          rot, tilt, psi, xshift, yshift;
    double           weight;
    bool            mirror;
    FileName        fn_img, fn_out, fn_in;
    ImageXmipp      img;
    int             levels,round_shifts;

// Check command line options ===========================================
    try
    {
        try
        {
        	round_shifts = checkParameter(argc, argv, "-round_shifts");
        	levels = textToInteger(getParameter(argc, argv, "-levels", "0"));
        	fn_in = getParameter(argc, argv, "-i");
        }
        catch (Xmipp_error XE)
        {
            std::cout << XE;
            Usage();
        }

        MetaData SF(fn_in);
        SF.removeObjects( MDL_ENABLED, -1 );
        std::vector< MetaDataLabel >::iterator strIt;

		long int ret=SF.firstObject();
		if(ret==MetaData::NO_OBJECTS_STORED)
		{
			std::cerr << "Empty inputFile File\n";
			exit(1);
		}
		do
		{
			SF.getValue( MDL_IMAGE, fn_img);
			if (fn_img=="") break;
			img.read(fn_img);
			img.clear_fFlag_flag();
			for( strIt  = SF.activeLabels.begin();
				 strIt != SF.activeLabels.end();
				 strIt ++ )
			{
				//std::cout << MetaDataContainer::decodeLabel(*strIt) << std::endl;


				switch ((*strIt)) {
					case MDL_ANGLEROT:
						SF.getValue( MDL_ANGLEROT, rot);
						img.set_rot(rot);
						break;
					case MDL_ANGLETILT:
						SF.getValue( MDL_ANGLETILT, tilt);
						img.set_tilt(tilt);
						break;
					case MDL_ANGLEPSI:
						SF.getValue( MDL_ANGLEPSI, psi);
						img.set_psi(psi);
						break;
					case MDL_SHIFTX:
						SF.getValue( MDL_SHIFTX, xshift);
					    if (levels != 0)
							xshift /= pow(2.0, levels);
		                if (round_shifts)
		                    xshift = ROUND(xshift);
						img.set_Xoff(xshift);
						break;
					case MDL_SHIFTY:
						SF.getValue(MDL_SHIFTY, yshift);
					    if (levels != 0)
							yshift /= pow(2.0, levels);
		                if (round_shifts)
		                    yshift = (float)ROUND(yshift);
						img.set_Yoff(yshift);
						break;
					case MDL_WEIGHT:
						SF.getValue( MDL_WEIGHT, weight);
						img.set_weight(weight);
						break;
					case MDL_FLIP:
						SF.getValue( MDL_FLIP, mirror);
						img.set_flip(mirror);
						break;
				}
			}
			img.write(fn_img);
		}
		while (SF.nextObject()!= MetaData::NO_MORE_OBJECTS);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }

}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    printf("Purpose:\n");
    printf(" Set the geometric transformation (angles & shifts) in the header of 2D-images.\n");
    printf("Usage:\n");
    printf("   header_assign  \n");
    printf("        -i <docfile>       : input metaData file\n");
    printf("       [-round_shifts]     : Round shifts to integers \n");
    printf("       [-levels <n=0>]     : Levels of pyramidal reduction, n=1, 2, ...\n");
    exit(1);
}

