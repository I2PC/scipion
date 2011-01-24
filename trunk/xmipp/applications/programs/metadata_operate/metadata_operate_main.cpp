/***************************************************************************
 *
 * Authors:     Sjors Scheres (scheres@cnb.csic.es)
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

class ProgMetadataOperate: public XmippProgram
{
    FileName fn_in, fn_out;
    std::string expressionStr;
    MetaData mdIn;

protected:

    void defineParams()
    {
        addUsageLine("Perform operations on MetaData columns. See examples below.");
        addExampleLine("  xmipp_metadata_operate  -i a.doc -o b.doc -e  \"angleRot=(angleRot*3.1416/180.)\"  ");
        addExampleLine("  xmipp_metadata_operate  -i a.doc -o b.doc -e  \"image=replace(image, 'xmp','spi')\"  ");
        addParamsLine("  -i <inputMetadata>     :MetaData input file name       ");
        addParamsLine("     alias --input;");
        addParamsLine("  -o <outputMetadata=\"/dev/stdout\"> :MetaData output file name, by default print to screen");
        addParamsLine("     alias --output;");
        addParamsLine("  -e <expression>                     :any valid operation in SQLite (expression must be between quotes) ");
        addParamsLine("     alias --expression;");
    }

    void readParams()
    {
        fn_in = getParam("-i");
        fn_out = getParam("-o");
        expressionStr = getParam("-e");
    }

public:
    void run()
    {
        mdIn.read(fn_in);
        mdIn.operate(expressionStr);
        mdIn.write(fn_out);
    }

}
;//end of class ProgMetadataOperate

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
        ProgMetadataOperate program;
        program.read(argc, argv);
        program.tryRun();
}
