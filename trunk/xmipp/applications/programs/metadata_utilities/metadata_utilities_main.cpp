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

    WriteModeMetaData mode;
protected:
    void defineParams()
    {
        addUsageLine ("Perform several operation over the metadata files");

        addParamsLine("  [--label <l1> ]                 : metadata label");
        addParamsLine("     alias -l;");

        addParamsLine(" [--mode+ <mode=overwrite>]   : Metadata writing mode.");
        addParamsLine("    where <mode>");
        addParamsLine("     overwrite   : Replace the content of the file with the Metadata");
        addParamsLine("     append      : Write the Metadata as a new block, removing the old one");

        addParamsLine("  [--expression <e1> ]                 : constrain applied in select");
        addParamsLine("     alias -e;");

        addParamsLine("   [-o  <md>]                          : Name of output metadata file");

        addParamsLine("   --union  <md11> <md22>             : union of metadata files md1 and md2");
        addParamsLine("     alias -u;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --intersection <md11> <md22>        : Intersection of md1 and md2");
        addParamsLine("     alias -i;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --subtraction <md11> <md22>         : subtraction of md1 and md2");
        addParamsLine("     alias -s;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --join <md11> <md22>                : inner join of md1 and md2 using label l1");
        addParamsLine("     alias -j;");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --sort <md11>                      : sort metadata md1 using label l1");
        addParamsLine("                                     : for sorting according to a component of a vector label");
        addParamsLine("                                     : use label:col, e.g., NMADisplacements:0");
        addParamsLine("                                     : The first column is column number 0");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --convert2db <md11>                : convert metadata to sqlite database");

        addParamsLine("or --copy  <md11> <path>               : copy files in metadata md1 to directory path (file names at lable column)");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --move  <md11> <path>               : move files in metadata md1 to directory path (file names at lable column)");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --delete  <md11>                    : delete files in metadata md1 (file names at label column)");
        addParamsLine("           requires --label;                                                         ");

        addParamsLine("or --select  <md11> <exp>       : create new metadata with those entries that satisfy the expression 'exp'");
        addParamsLine("           requires -o;                                                         ");

        addParamsLine("or --count  <md11>              : for each value of a given label create new metadata with the number of times the value appears");
        addParamsLine("           requires --label, -o;                                                         ");

        addParamsLine("or --size  <md11>              : metadata size");

        addExampleLine(" Concatenate two metadatas.", false);
        addExampleLine ("   xmipp_metadata_utilities --union         mD1.doc mD2.doc  -o out.doc --label image");
        addExampleLine(" Intersect two metadatas.", false);
        addExampleLine ("   xmipp_metadata_utilities --intersection  mD1.doc mD2.doc  -o out.doc --label image");
        addExampleLine(" Substract two metadatas.", false);
        addExampleLine ("   xmipp_metadata_utilities --subtraction   mD1.doc mD2.doc  -o out.doc --label image");
        addExampleLine(" Combine columns from both metadatas.", false);
        addExampleLine ("   xmipp_metadata_utilities --join j1.doc   mD1.doc          -o out.doc --label image");
        addExampleLine(" Sort the elements in metadata.", false);
        addExampleLine ("   xmipp_metadata_utilities --sort          mD1.doc          -o out.doc --label image");
        addExampleLine(" Dump metadata content to Sqlite3 database. (use xmipp_sqlite3 to visualize results)", false);
        addExampleLine ("   xmipp_metadata_utilities --convert2db    mD1.doc          -o out.db; xmipp_sqlite3 out.db");
        addExampleLine(" Copy files in metadata to a location.", false);
        addExampleLine ("   xmipp_metadata_utilities --copy mD1.doc kk                -o out.doc --label image ");
        addExampleLine(" Delete files in metadata.", false);
        addExampleLine ("   xmipp_metadata_utilities --delete out.doc                            --label image");
        addExampleLine(" Select elements in metadata that satisfy a given constrain.", false);
        addExampleLine ("   xmipp_metadata_utilities --select mD1.doc \"anglePsi > 0 AND shiftX > -0.5\" -o out.doc");
        addExampleLine(" Count number of images per CTF", false);
        addExampleLine ("   xmipp_metadata_utilities --count mD1.doc  -o out.doc --label CTFModel");
        addExampleLine(" Metadata Size", false);
        addExampleLine ("   xmipp_metadata_utilities --size mD1.doc");

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
        _convert2db=10,
        _count=11,
        _size=12
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

        else if (strcmp(s,"count") == 0)
            operationType = _count;

        else if (strcmp(s,"size") == 0)
            operationType = _size;
    }

    void readParams()
    {
        if (checkParam("-o"))
            outFileName = getParam("-o");
        if (checkParam("--label"))
            _label = getParam("--label");
        mode = metadataModeConvert(getParam("--mode"));

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

        else if (checkParam("--count"))
        {
            encode("count");
            inFileName1 = getParam("--count",0);
            outFileName  = getParam("-o");
        }
        else if (checkParam("--size"))
        {
            encode("size");
            inFileName1 = getParam("--size",0);
        }
    }
public:
    void run()
    {
        FileName inFnImg,outFnImg;
        size_t id;

        switch (operationType)
        {
        case _union:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            inMD1.unionDistinct(inMD2, MDL::str2Label(_label));
            inMD1.write(outFileName,mode);
            break;
        case _intersection:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            inMD1.intersection(inMD2, MDL::str2Label(_label));
            inMD1.write(outFileName,mode);
            break;
        case _subtraction:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            inMD1.subtraction(inMD2, MDL::str2Label(_label));
            inMD1.write(outFileName,mode);
            break;
        case _sort:
            inMD1.read(inFileName1);
            outMD.sort(inMD1, _label);
            outMD.write(outFileName,mode);
            break;
        case _join:
            inMD1.read(inFileName1);
            inMD2.read(inFileName2);
            MDSql::dumpToFile("pp.db");
            outMD.join(inMD1,inMD2, MDL::str2Label(_label));
            outMD.write(outFileName,mode);
            break;
        case _copy:
            inMD1.read(inFileName1);
            //create dir
            if (!exists(tmpFileName))
                if (mkpath(tmpFileName, 0755) != 0)
                    REPORT_ERROR(ERR_IO_NOPERM, "Run: Cannot create directory "+ tmpFileName);
            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
            {
                inMD1.getValue(MDL::str2Label(_label),inFnImg, __iter.objId);
                outFnImg = inFnImg.removeDirectories();
                id = outMD.addObject();
                outMD.setValue(MDL::str2Label(_label),outFnImg, id);
                outFnImg = tmpFileName + "/" + outFnImg;
                inFnImg.copyFile(outFnImg);
            }
            outMD.write(tmpFileName+"/"+outFileName,mode);
            break;
        case _move:
            inMD1.read(inFileName1);
            //create dir
            if (!exists(tmpFileName))
                if (mkpath(tmpFileName, 0755) != 0)
                    REPORT_ERROR(ERR_IO_NOPERM, "Run: Cannot create directory "+ tmpFileName);
            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
            {
                inMD1.getValue(MDL::str2Label(_label),inFnImg, __iter.objId);
                outFnImg = inFnImg.removeDirectories();
                id = outMD.addObject();
                outMD.setValue(MDL::str2Label(_label),outFnImg, id);
                outFnImg = tmpFileName + "/" + outFnImg;
                rename(inFnImg.c_str(),outFnImg.c_str());
            }
            outMD.write(tmpFileName+"/"+outFileName,mode);
            break;
        case _delete:
            inMD1.read(inFileName1);
            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
            {
                inMD1.getValue(MDL::str2Label(_label),inFnImg, __iter.objId);
                remove(inFnImg.c_str());
                std::cerr << "Remove file: " << inFnImg <<std::endl;
            }
            break;
        case _select:
            {
                inMD1.read(inFileName1);
                outMD.importObjects(inMD1, MDExpression(expression));
                outMD.write(outFileName,mode);
            }
            break;
        case _count:
            {
                inMD1.read(inFileName1);
                outMD.aggregate(inMD1, AGGR_COUNT,MDL::str2Label(_label),MDL_CTFMODEL,MDL_COUNT);
                outMD.write(outFileName,mode);
            }
            break;
        case _size:
            {
                inMD1.read(inFileName1);
                std::cout << inFileName1 + " size is: " << inMD1.size() << std::endl;
            }
            break;
        case _convert2db:
            {
                inMD1.read(inFileName1);
                MDSql::dumpToFile(outFileName);
            }
            break;

        default:
            REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown operation.");
        }
    }
}
;

int main(int argc, char **argv)
{

    ProgMetadataUtilities program;
    program.read(argc, argv);
    program.tryRun();


}
