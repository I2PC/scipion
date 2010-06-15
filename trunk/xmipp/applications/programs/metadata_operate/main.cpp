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

void Usage();

int main(int argc, char **argv)
{
    FileName fn_in, fn_out;
    std::string expressionStr;
    MetaData mdIn;

    try
    {
        fn_in = getParameter(argc, argv, "-i");
        fn_out = getParameter(argc, argv, "-o", "/dev/stdout");
        expressionStr = getParameter(argc, argv, "-e");
    }
    catch (Xmipp_error)
    {
        Usage();
        exit(1);
    }

    try
    {
        mdIn.read(fn_in);
        mdIn.operate(expressionStr);
        mdIn.write(fn_out);
    }
    catch (Xmipp_error)
    {
        std::cerr << "ERROR, exiting..." << std::endl;
        exit(1);
    }

}

void Usage()
{
    std::cout << "Usage: metadata_operate [options]\n"
            << "  Example:  \n"
            << " operate_metadata  -i a.doc -o b.doc -e  \"angleRot=(angleRot*3.1416/180.)\"  \n"
            << " operate_metadata  -i a.doc -o b.doc -e  \"image=replace(image, 'xmp','spi')\"  \n"
            << " Options:\n"
            << "  -i inputMetadata     MetaData input file name       \n"
            << "  -o outputMetadata    MetaData output file name, by default print to screen\n"
            << "  -e expression,       any valid operation in SQLite (expression must be between quotes) \n"

    ;
}
