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
#include <data/metadata.h>
#include <data/progs.h>

class ProgHeaderReset: public XmippMetadataProgram
{
private:
    bool            tiltSeries;
    double          firstAngle, angularStep, angle;

protected:
    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();

        addUsageLine("Reset the geometric transformation (angles & shifts) in the header of 2D-images.");
        addParamsLine("   [--tilt_series <firstAngle> <angleStep>]: Assign a regularly spaced angular distribution.");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        tiltSeries = checkParam("--tilt_series");
        if (tiltSeries)
        {
            firstAngle = getDoubleParam("--tilt_series", 0);
            angularStep = getDoubleParam("--tilt_series", 1);
        }
    }

    void show()
    {
        std::cout << " Resetting all angles, origin offsets, weights and mirror flags to zero ... " << std::endl;

        if (tiltSeries)
        {
            std::cout << "Setting the tilt angles to a tilt series\n"
            << "First angle = " << firstAngle << std::endl
            << "Angular step = " << angularStep << std::endl;
            angle = firstAngle;
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg, mdIn, objId); //read data and header
        img.clearHeader();

        if (tiltSeries)
        {
            img.setTilt(angle);
            angle += angularStep;
        }
        double daux = (double)1.;
        img.setWeight(daux);
        img.write(fnImgOut);
    }
}
;// end of class ProgHeaderReset

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgHeaderReset program;
    program.read(argc, argv);
    program.tryRun();
}


