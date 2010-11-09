/***************************************************************************
 * Authors:     roberto marabini roberto@cnb.csic.es
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your param) any later version.
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

#include <data/argsparser.h>
#include <data/program.h>
#include <string.h>
#include <data/metadata.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>

class ProgMetadataUtilities: public XmippProgram
{
private:

protected:
    void defineParams()
    {
        addUsageLine ("Perform several operation over the metadata files");

        addParamsLine("  [--label <l1> ]                 : metadata label");
        addParamsLine("     alias -l;");

        addParamsLine("  [--expression <e1> ]                 : constrain applied in select");
        addParamsLine("     alias -e;");

        addParamsLine("   [-o  <w1>]                          : Name of output metadata file");

        addParamsLine("   --union  <md1> <md2>             : union of metadata files md1 and md2");
        addParamsLine("     alias -u;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --intersection <md1> <md2>        : Intersection of md1 and md2");
        addParamsLine("     alias -i;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --subtraction <md1> <md2>         : subtraction of md1 and md2");
        addParamsLine("     alias -s;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --join <md1> <md2>                : inner join of md1 and md2 using label l1");
        addParamsLine("     alias -j;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --sort <md1>                      : sort metadata md1 using label l1");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --convert2db <md1>                : convert metadata to sqlite database");

        addParamsLine("or --copy  <md1> <path>               : copy files in metadata md1 to directory path (file names at lable column)");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --move  <md1> <path>               : move files in metadata md1 to directory path (file names at lable column)");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --delete  <md1>                    : delete files in metadata md1 (file names at label column)");
        addParamsLine("           requires --label;                                                         ");

        addParamsLine("or --select  <md1> <exp>       : create new metadata with those entries that satisfy the expression 'exp'");
        addParamsLine("           requires -o;                                                         ");

        addUsageLine ("Examples:");
        addUsageLine ("   xmipp_metadata_utilities --union         mD1.doc mD2.doc  -o out.doc --label image");
        addUsageLine ("   xmipp_metadata_utilities --intersection  mD1.doc mD2.doc  -o out.doc --label image");
        addUsageLine ("   xmipp_metadata_utilities --subtraction   mD1.doc mD2.doc  -o out.doc --label image");
        addUsageLine ("   xmipp_metadata_utilities --join j1.doc   mD1.doc          -o out.doc --label image");
        addUsageLine ("   xmipp_metadata_utilities --sort          mD1.doc          -o out.doc --label image");
        addUsageLine ("   xmipp_metadata_utilities --convert2db    mD1.doc          -o out.db; xmipp_sqlite3 out.db");
        addUsageLine ("   xmipp_metadata_utilities --copy mD1.doc kk                -o out.doc --label image ");
        addUsageLine ("   xmipp_metadata_utilities --delete out.doc                            --label image");
        addUsageLine ("   xmipp_metadata_utilities --select mD1.doc \"anglePsi > 0 AND shiftX > -0.5\" -o out.doc");

    }
    typedef enum {
        _unknown=0,
        _union=1,
        _intersection=2,
        _subtraction=3,
        _join=4,
        _copy=5,
        _move=6,
        _delete=7,
        _select=8,
        _sort=9,
        _convert2db=10
    } OperationType;
    OperationType operationType;
    MetaData inMD1, inMD2, outMD;
    FileName inFileName1, inFileName2, outFileName, tmpFileName;
    std::string _label;
    std::string expression;
    double min,max;

    void encode(const char * s)
    {
        operationType = _unknown;

        if (strcmp(s,"union") == 0)
            operationType = _union;

        else if (strcmp(s,"intersection") == 0)
            operationType = _intersection;

        else if (strcmp(s,"subtraction") == 0)
            operationType = _subtraction;

        else if (strcmp(s,"join") == 0)
            operationType = _join;

        else if (strcmp(s,"copy") == 0)
            operationType = _copy;

        else if (strcmp(s,"move") == 0)
            operationType = _move;

        else if (strcmp(s,"delete") == 0)
            operationType = _delete;

        else if (strcmp(s,"select") == 0)
            operationType = _select;

        else if (strcmp(s,"sort") == 0)
            operationType = _sort;

        else if (strcmp(s,"convert2db") == 0)
            operationType = _convert2db;
    }

    void readParams()
    {
        if (checkParam("-o"))
            outFileName = getParam("-o");
        if (checkParam("--label"))
            _label = getParam("--label");

        if (checkParam("--union"))
        {
            encode("union");
            inFileName1 = getParam("--union",0);
            inFileName2 = getParam("--union",1);
        }

        else if (checkParam("--intersection"))
        {
            encode("intersection");
            inFileName1 = getParam("--intersection",0);
            inFileName2 = getParam("--intersection",1);
        }

        else if (checkParam("--subtraction"))
        {
            encode("subtraction");
            inFileName1 = getParam("--subtraction",0);
            inFileName2 = getParam("--subtraction",1);
        }

        else if (checkParam("--join"))
        {
            encode("join");
            inFileName1 = getParam("--join",0);
            inFileName2 = getParam("--join",1);
        }

        else if (checkParam("--sort"))
        {
            encode("sort");
            inFileName1  = getParam("--sort",0);
            outFileName  = getParam("-o");
        }

        else if (checkParam("--convert2db"))
        {
            encode("convert2db");
            inFileName1 = getParam("--convert2db",0);
            inMD1.read(inFileName1);
            MDSql::dumpToFile(outFileName);
        }

        else if (checkParam("--move"))
        {
            encode("move");
            inFileName1 = getParam("--move",0);
            tmpFileName = getParam("--move",1);
        }

        else if (checkParam("--copy"))
        {
            encode("copy");
            inFileName1 = getParam("--copy",0);
            tmpFileName = getParam("--copy",1);
        }

        else if (checkParam("--delete"))
        {
            encode("delete");
            inFileName1 = getParam("--delete",0);
        }

        else if (checkParam("--select"))
        {
            encode("select");
            inFileName1 = getParam("--select",0);
            expression = getParam("--select",1);
        }
    }
public:
    void run()
    {
        switch (operationType)
        {
        case _union:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            inMD1.unionDistinct(inMD2, MDL::str2Label(_label));
            inMD1.write(outFileName);
            break;
        case _intersection:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            inMD1.intersection(inMD2, MDL::str2Label(_label));
            inMD1.write(outFileName);
            break;
        case _subtraction:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            inMD1.subtraction(inMD2, MDL::str2Label(_label));
            inMD1.write(outFileName);
            break;
        case _sort:
            inMD1.read(inFileName1);
            outMD.sort(inMD1, MDL::str2Label(_label));
            outMD.write(outFileName);
            break;
        case _join:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            MDSql::dumpToFile("pp.db");
            outMD.join(inMD1,inMD2, MDL::str2Label(_label));
            outMD.write(outFileName);
            break;
        case _copy:
            inMD1.read(inFileName1);
            //create dir
            if (!exists(tmpFileName))
            	if (mkpath(tmpFileName, 0755) != 0)
            		REPORT_ERROR(ERR_IO_NOPERM, "Run: Cannot create directory "+ tmpFileName);
            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
            {
                FileName inFnImg,outFnImg;
                inMD1.getValue(MDL::str2Label(_label),inFnImg);
                outFnImg = inFnImg.removeDirectories();
                outMD.addObject();
                outMD.setValue(MDL::str2Label(_label),outFnImg);
                outFnImg = tmpFileName + "/" + outFnImg;
                inFnImg.copyFile(outFnImg);
            }
            outMD.write(tmpFileName+"/"+outFileName);
            break;
        case _move:
            inMD1.read(inFileName1);
            //create dir
            if (!exists(tmpFileName))
            	if (mkpath(tmpFileName, 0755) != 0)
            		REPORT_ERROR(ERR_IO_NOPERM, "Run: Cannot create directory "+ tmpFileName);
            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
            {
                FileName inFnImg,outFnImg;
                inMD1.getValue(MDL::str2Label(_label),inFnImg);
                outFnImg = inFnImg.removeDirectories();
                outMD.addObject();
                outMD.setValue(MDL::str2Label(_label),outFnImg);
                outFnImg = tmpFileName + "/" + outFnImg;
                rename(inFnImg.c_str(),outFnImg.c_str());
            }
            outMD.write(tmpFileName+"/"+outFileName);
            break;
        case _delete:
            inMD1.read(inFileName1);
            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
            {
                FileName inFnImg;
                inMD1.getValue(MDL::str2Label(_label),inFnImg);
                remove(inFnImg.c_str());
                std::cerr << "Remove file: " << inFnImg <<std::endl;
            }
            break;
        case _select:
            {
                inMD1.read(inFileName1);
                outMD.importObjects(inMD1, MDExpression(expression));
                outMD.write(outFileName);
            }
            break;

        default:
            REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown operation");
        }
    }

}
;

int main(int argc, char **argv)
{
    try
    {
        ProgMetadataUtilities program;
        program.read(argc, argv);
        program.run();

    }
    catch (XmippError e)
    {
        std::cerr << e.msg <<std::endl;
    }
}
