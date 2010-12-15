/***************************************************************************
 *
 * Authors:     Roberto Marabini roberto@cnb.csic.es)
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

#include <data/metadata.h>
#include <data/args.h>
#include <data/program.h>

class ProgAngleMain: public XmippProgram
{
    FileName fn_in, fn_out;
    std::string expressionStr;
    MetaData _mdIn;

protected:

    void defineParams()
    {
        addUsageLine("Convert angles from freealign to xmipp.");
        addUsageLine("Example:  ");
        addUsageLine("  convert_angle  -i a.doc -o b.doc ");
        addParamsLine("  -i <inputMetadata>     :MetaData input file name       ");
        addParamsLine("     alias --input;");
        addParamsLine("  -o <outputMetadata=\"/dev/stdout\"> :MetaData output file name, by default print to screen");
        addParamsLine("     alias --output;");
    }

    void readParams()
    {
        fn_in  = getParam("-i");
        fn_out = getParam("-o");
    }

public:
    void run()
    {
        double rot1,  tilt1,  psi1;
        double rot2,  tilt2,  psi2;
        _mdIn.read(fn_in);
        //test substract 90-omega and use all the symmetries on that
        FOR_ALL_OBJECTS_IN_METADATA(_mdIn)
        {
            _mdIn.getValue(MDL_ANGLEROT,rot1);
            _mdIn.getValue(MDL_ANGLETILT,tilt1);
            _mdIn.getValue(MDL_ANGLEPSI,psi1);

            rot2  = rot1;
            tilt2 = tilt1;
            psi2  = psi1;

            //EULER MATRICES
            Matrix2D<double> in(3, 3), temp(3,3), out(3,3);
            Euler_angles2matrix(rot1, tilt1, psi1, in);
            Euler_angles2matrix(90, 0, 0, temp);
            out = temp * in;
            Euler_angles2matrix(0, 90, 0, temp);
            out = out * temp;
            Euler_matrix2angles(out, rot2, tilt2, psi2);

            _mdIn.setValue(MDL_ANGLEROT,rot2);
            _mdIn.setValue(MDL_ANGLETILT,tilt2);
            _mdIn.setValue(MDL_ANGLEPSI,psi2);

        }

        _mdIn.write(fn_out);
    }

}
;//end of class ProgMetadataOperate

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    try
    {
        ProgAngleMain program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }
}
