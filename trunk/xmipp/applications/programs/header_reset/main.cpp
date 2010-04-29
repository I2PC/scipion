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
#include <data/selfile.h>
#include <data/metadata.h>

void Usage();

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ImageXmipp      img;
    FileName        fn_input;
    bool            tiltSeries;
    double          firstAngle, angularStep;
    MetaData SF;

    try
    {
        fn_input = getParameter(argc, argv, "-i");
        if (Is_ImageXmipp(fn_input))
        {
            SF.addObject();
            SF.setValue( MDL_IMAGE, fn_input);
            SF.setValue( MDL_ENABLED, 1);
        }
        else
            SF.read( fn_input ,NULL);

        tiltSeries=checkParameter(argc,argv,"-tiltSeries");
        if (tiltSeries)
        {
            int i=paremeterPosition(argc,argv,"-tiltSeries");
            if (i+2>=argc)
                REPORT_ERROR(1,"Not enough parameters after -tiltSeries");
            firstAngle=textToFloat(argv[i+1]);
            angularStep=textToFloat(argv[i+2]);
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
    }

    try
    {
        std::cerr << " Resetting all angles, origin offsets, weights and mirror flags to zero ... " << std::endl;
        if (tiltSeries)
            std::cerr << "Setting the tilt angles to a tilt series\n"
                      << "First angle=" << firstAngle << std::endl
                      << "Angular step=" << angularStep << std::endl;
        long int ret=SF.firstObject();
        double angle=firstAngle;
        do
        {
            FileName fn_img;
            SF.getValue( MDL_IMAGE, fn_img);
            if (fn_img=="") break;
            img.read(fn_img);
            img.clear_header();
            if (tiltSeries)
            {
                img.set_tilt(angle);
                angle+=angularStep;
            }
            img.write(img.name());
        }
        while (SF.nextObject()!= MetaData::NO_MORE_OBJECTS);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
    }
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    printf("Purpose:\n");
    printf(" Reset the geometric transformation (angles & shifts) in the header of 2D-images.\n");
    printf("Usage:\n");
    printf("   header_reset \n");
    printf("    -i                                   : metaDataFile with images or individual image\n");
    printf("   [-tiltSeries <firstAngle> <angleStep>]: Assign a regularly spaced angular distribution\n");
    exit(1);
}

