/***************************************************************************
 *
 * Authors:     J.M de la Rosa (jmdelarosa@cnb.csic.es)
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
#include <data/program.h>

class ProgMetadataImport: public XmippProgram
{
    FileName fn_in, fn_out;
    String sep;
    StringVector labels;
    MetaData md;

protected:

    void defineParams()
    {
        addUsageLine("Import metadata from plain text files.");
        addExampleLine("  xmipp_metadata_operate  -i a.doc -o b.doc -e  \"angleRot=(angleRot*3.1416/180.)\"  ");
        addExampleLine("  xmipp_metadata_operate  -i a.doc -o b.doc -e  \"image=replace(image, 'xmp','spi')\"  ");
        addParamsLine("  -i <text_file>     :Input file text file       ");
        addParamsLine("     alias --input;");
        addParamsLine(" [ -o <output_metadata>]     :If not provided, the resulting metadata will be printed on screen");
        addParamsLine("     alias --output;");
        addParamsLine("  -l <...>                   :labels to be imported");
        addParamsLine("     alias --labels;");
        addParamsLine(" [ -m <metadata> ]           :merge the imported metadata to an existing one");
        addParamsLine("     alias --merge;");
        addParamsLine(" [ -s <sep=\" \">]             :Separator to be used, default is space");
        addParamsLine("     alias --separator;");

    }

    void readParams()
    {
        fn_in = getParam("-i");
        sep = getParam("-s");
        getListParam("-l", labels);
        if (labels.empty())
            REPORT_ERROR(ERR_ARG_BADCMDLINE, "You should provide at least one label to import");
    }

public:
    void run()
    {
      if (checkParam("-m"))
      {
        md.read(getParam("-m"));
        md.addPlain(fn_in, labels, sep);
      }
      else
        md.readPlain(fn_in, labels, sep);

        if (checkParam("-o"))
          md.write(getParam("-o"));
        else
          md.write(std::cout);
    }

}
;//end of class ProgMetadataOperate

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgMetadataImport program;
    program.read(argc, argv);
    return program.tryRun();
}
