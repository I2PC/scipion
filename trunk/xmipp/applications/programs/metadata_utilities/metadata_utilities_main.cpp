/***************************************************************************
 * Authors:     Roberto Marabini roberto@cnb.csic.es
 *              J.M. de la Rosa  jmdelarosa@cnb.csic.es
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


typedef enum { UNIFORM, GAUSSIAN, STUDENT } RandMode;

/** MDGenerator to generate random values on columns */
class MDRandGenerator: public MDValueGenerator
{
protected:
    double op1, op2, op3;
    RandMode mode;

    inline double getRandValue()
    {
        switch (mode)
        {
        case UNIFORM:
            return rnd_unif(op1, op2);
        case GAUSSIAN:
            return rnd_gaus(op1, op2);
        case STUDENT:
            return rnd_student_t(op3, op1, op2);
        }
    }
public:
    MDRandGenerator(double op1, double op2, const String &mode, double op3=0.)
    {
        static bool randomized = false;

        if (!randomized)//initialize random seed just once
        {
            randomize_random_generator();
            randomized = true;
        }
        this->op1 = op1;
        this->op2 = op2;
        this->op3 = op3;
        if (mode == "uniform")
            this->mode = UNIFORM;
        else if (mode == "gaussian")
            this->mode = GAUSSIAN;
        else if (mode == "student")
            this->mode = STUDENT;
        else
            REPORT_ERROR(ERR_PARAM_INCORRECT, formatString("Unknown random type '%s'", mode.c_str()));

    }

    bool fillValue(MetaData &md, size_t objId)
    {
        double aux = getRandValue();
        md.setValue(label, aux, objId);
    }

}
;//end of class MDRandGenerator

/** Class to fill columns with constant values */
class MDConstGenerator: public MDValueGenerator
{
public:
    String value;

    MDConstGenerator(const String &value)
    {
        this->value = value;
    }
    bool fillValue(MetaData &md, size_t objId)
    {
        md.setValueFromStr(label, value, objId);
    }
}
;//end of class MDConstGenerator

/** Class to fill columns with another metadata in row format */
class MDExpandGenerator: public MDValueGenerator
{
public:
    MetaData expMd;
    FileName fn;
    MDRow row;

    bool fillValue(MetaData &md, size_t objId)
    {
        if (md.getValue(label, fn, objId))
        {
            std::cerr << "expanding " << fn << std::endl;
            expMd.read(fn);
            if (expMd.getColumnFormat() || expMd.isEmpty())
                REPORT_ERROR(ERR_VALUE_INCORRECT, "Only can expand non empty and row formated metadatas");
            expMd.getRow(row, expMd.firstObject());
            md.setRow(row, objId);
        }
        else
            REPORT_ERROR(ERR_MD_BADLABEL, formatString("Can't expand missing label '%s'", MDL::label2Str(label).c_str()));
    }
}
;//end of class MDExpandGenerator

class ProgMetadataUtilities: public XmippProgram
{
private:
    WriteModeMetaData mode;
    FileName fn_in, fn_out, fn_md2;
    MetaData mdIn, md2;
    MDLabel label;
    std::vector<MDLabel> labels;
    String operation;
    bool doWrite;

protected:
    void defineParams()
    {
        addUsageLine ("Perform several operation over the metadata files.");
        addUsageLine ("If the -o option is not used, the original metadata");
        addUsageLine ("will be modified.");

        addParamsLine(" -i <metadata>                          : Input metadata file");
        addParamsLine("   [-o  <metadata>]                    : Output metadata file, if not provided result will be printed on screen");

        addParamsLine("  --set <set_operation> <label=image>    : Set operations");
        addParamsLine("         where <set_operation>");
        addParamsLine("   union  <md2>             : Union with metadata md2");
        addParamsLine("   intersection <md2>       : Intersection with metadata md2");
        addParamsLine("   subtraction <md2>        : Subtraction with metadata md2");
        addParamsLine("   join <md2>               : Inner join with md2 using label l1");
        addParamsLine("   merge <md2>              : Merge columns with md2, label is ignored");
        addParamsLine("                            : Both metadatas should have same size, and elements should be in same order,");
        addParamsLine("                            : if not, you should use 'join' instead, but this constrain having a common label");
        addParamsLine("   sort                     : Sort metadata using label l1");
        addParamsLine("                            : for sorting according to a component of a vector label");
        addParamsLine("                            : use label:col, e.g., NMADisplacements:0");
        addParamsLine("                            : The first column is column number 0");
        addParamsLine("           alias -s;                                             ");

        addParamsLine("or --operate <operation>     : Operations on the metadata structure");
        addParamsLine("         where <operation>");
        addParamsLine("    add_column <labels>                 : Add some columns(label list) to metadata");
        addParamsLine("    drop_column <labels>                : Drop some columns(label list) from metadata");
        addParamsLine("    modify_values <expression>                 : Use an SQLite expression to modify the metadata");
        addParamsLine("                                        : This option requires knowledge of basic SQL syntax(more specific SQLite");

        addParamsLine("or  --file <file_operation>     : File operations");
        addParamsLine("         where <file_operation>");
        addParamsLine("   copy <directory> <label=image>  : Copy files in metadata md1 to directory path (file names at label column)");
        addParamsLine("   move <directory>  <label=image> : Move files in metadata md1 to directory path (file names at label column)");
        addParamsLine("   delete  <label=image>      : Delete files in metadata md1 (file names at label column)");
        addParamsLine("   convert2db                 : Convert metadata to sqlite database");

        addParamsLine("or --query <query_operation>   : Query operations");
        addParamsLine("         where <query_operation>");
        addParamsLine("   select <expression>        : Create new metadata with those entries that satisfy the expression");
        addParamsLine("   count  <label>             : for each value of a given label create new metadata with the number of times the value appears");
        addParamsLine("   size                       : print Metadata size");

        addParamsLine("or --fill <labels> <fill_mode>                  : Fill a column values(should be of same type)");
        addParamsLine("   where <fill_mode>");
        addParamsLine("     constant  <value>                        : Fill with a constant value");
        addParamsLine("     rand_uniform  <a=0.> <b=1.>              : Follow a uniform distribution between a and b");
        addParamsLine("     rand_gaussian <mean=0.> <stddev=1.>      : Follow a gaussian distribution with mean and stddev");
        addParamsLine("     rand_student  <mean=0.> <stddev=1.> <df=3.> : Follow a student distribution with mean, stddev and df degrees of freedom.");
        addParamsLine("     expand                                   : Treat the column as the filename of an row metadata and expand values");
        addParamsLine(" [--mode+ <mode=overwrite>]   : Metadata writing mode.");
        addParamsLine("    where <mode>");
        addParamsLine("     overwrite   : Replace the content of the file with the Metadata");
        addParamsLine("     append      : Write the Metadata as a new block, removing the old one");

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
        addExampleLine(" Add column to metadata.", false);
        addExampleLine ("   xmipp_metadata_utilities --addColumn     mD1.doc          -o out.doc --label image");
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
        addExampleLine(" Metadata randomize double values", false);
        addExampleLine ("   xmipp_metadata_utilities --randValues mD1.doc gaussian 0. 0.15 -o out.doc --label scale");

    }

    void readParams()
    {
        fn_in = getParam("-i");
        mdIn.read(fn_in);
        doWrite = checkParam("-o");
        fn_out = doWrite ? getParam("-o") : fn_in;
    }

    void doSet()
    {
        operation = getParam("--set", 0);
        int labelIndex = 1;
        if (operation != "sort")
        {
            labelIndex = 2;
            md2.read(getParam("--set", 1));
        }
        else
            md2 = mdIn;
        MDLabel label = MDL::str2Label(getParam("--set", labelIndex));
        if (operation == "union")
            mdIn.unionDistinct(md2, label);
        else if (operation == "intersection")
            mdIn.intersection(md2, label);
        else if (operation == "subtraction")
            mdIn.subtraction(md2, label);
        else if (operation == "join")
        {
            MetaData md;
            md.join(mdIn, md2, label);
            mdIn = md;
        }
        else if (operation == "merge")
            mdIn.merge(md2);
        else if (operation == "sort")
            mdIn.sort(md2, label);
    }//end of function doSet

    void doOperate()
    {
        operation = getParam("--operate", 0);

        if (operation != "modify_values")
        {
            MDL::str2LabelVector(getParam("--operate", 1), labels);
            for (int i = 0; i < labels.size(); ++i)
            {
                if (operation == "add_column")
                    mdIn.addLabel(labels[i]);
                else if (operation == "drop_column")
                    mdIn.removeLabel(labels[i]);
            }
        }
        else // modify_values
            mdIn.operate(getParam("--operate", 1));
    }//end of function doOperate

    void doFill()
    {

        MDL::str2LabelVector(getParam("--fill", 0), labels);

        if (labels.empty())
            REPORT_ERROR(ERR_PARAM_INCORRECT, "You should provide at least one label to fill out");

        operation = getParam("--fill", 1);
        MDValueGenerator * generator;

        // Select wich generator to use
        if (operation == "constant")
            generator = new MDConstGenerator(getParam("--fill", 2));
        else if (operation.find("rand_") == 0)
        {
            double op1 = getDoubleParam("--fill", 2);
            double op2 = getDoubleParam("--fill", 3);
            double op3 = 0.;
            String type = findAndReplace(operation, "rand_", "");
            if (type == "student")
                op3 = getDoubleParam("--fill", 4);
            generator = new MDRandGenerator(op1, op2, type, op3);
        }
        else if (operation == "expand")
            generator = new MDExpandGenerator();

        //Fill columns
        for (int i = 0; i < labels.size(); ++i)
        {
            generator->label = labels[i];
            generator->fill(mdIn);
        }

        delete generator;
    }//end of function doFill

    void doQuery()
    {
        operation = getParam("--query", 0);
        String expression;

        if (operation == "count")//note MDL_CTFMODEL is a dummy parameter
        {
            label = MDL::str2Label(getParam("--query", 1));
            md2 = mdIn;
            mdIn.aggregate(md2, AGGR_COUNT,label,MDL_CTFMODEL,MDL_COUNT);
        }
        else if (operation == "select")
        {
            expression = getParam("--query", 1);
            md2 = mdIn;
            mdIn.importObjects(md2, MDExpression(expression));
        }
        else if (operation == "size")
        {
            doWrite = false;
            std::cout << fn_in + " size is: " << mdIn.size() << std::endl;
        }
    }//end of function doQuery

    void doFile()
    {
        operation = getParam("--file", 0);

        if (operation == "convert2db")
        {
            doWrite = false;
            MDSql::dumpToFile(fn_out);
        }
        else
        {
            bool doDelete;
            FileName path, inFnImg, outFnImg;

            if (!(doDelete = operation == "delete"))//copy or move
            {
                path = getParam("--file", 1);
                if (!exists(path))
                    if (mkpath(path, 0755) != 0)
                        REPORT_ERROR(ERR_IO_NOPERM, (String)"Cannot create directory "+ path);
            }
            label = MDL::str2Label(getParam("--file", doDelete ? 1 : 2));
            doWrite = !doDelete;

            FOR_ALL_OBJECTS_IN_METADATA(mdIn)
            {

                mdIn.getValue(label, inFnImg, __iter.objId);

                if (doDelete)
                    remove(inFnImg.c_str());
                else
                {
                    outFnImg = inFnImg.removeDirectories();
                    mdIn.setValue(label, outFnImg, __iter.objId);
                    outFnImg = path + "/" + outFnImg;

                    if (operation == "copy")
                        inFnImg.copyFile(outFnImg);
                    else if (operation == "move")
                        rename(inFnImg.c_str(), outFnImg.c_str());
                }
            }
            fn_out = path + "/" + fn_out.removeDirectories();
        }
    }

public:
    void run()
    {
        if (checkParam("--set"))
            doSet();
        else if (checkParam("--operate"))
            doOperate();
        else if (checkParam("--file"))
            doFile();
        else if (checkParam("--query"))
            doQuery();
        else if (checkParam("--fill"))
            doFill();

        if (doWrite)
            mdIn.write(fn_out);
    }
    //        case _copy:
    //            inMD1.read(inFileName1);
    //            //create dir
    //            if (!exists(tmpFileName))
    //                if (mkpath(tmpFileName, 0755) != 0)
    //                    REPORT_ERROR(ERR_IO_NOPERM, "Run: Cannot create directory "+ tmpFileName);
    //            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
    //            {
    //                inMD1.getValue(MDL::str2Label(_label),inFnImg, __iter.objId);
    //                outFnImg = inFnImg.removeDirectories();
    //                id = outMD.addObject();
    //                outMD.setValue(MDL::str2Label(_label),outFnImg, id);
    //                outFnImg = tmpFileName + "/" + outFnImg;
    //                inFnImg.copyFile(outFnImg);
    //            }
    //            outMD.write(tmpFileName+"/"+outFileName,mode);
    //            break;
    //        case _move:
    //            inMD1.read(inFileName1);
    //            //create dir
    //            if (!exists(tmpFileName))
    //                if (mkpath(tmpFileName, 0755) != 0)
    //                    REPORT_ERROR(ERR_IO_NOPERM, "Run: Cannot create directory "+ tmpFileName);
    //            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
    //            {
    //                inMD1.getValue(MDL::str2Label(_label),inFnImg, __iter.objId);
    //                outFnImg = inFnImg.removeDirectories();
    //                id = outMD.addObject();
    //                outMD.setValue(MDL::str2Label(_label),outFnImg, id);
    //                outFnImg = tmpFileName + "/" + outFnImg;
    //                rename(inFnImg.c_str(),outFnImg.c_str());
    //            }
    //            outMD.write(tmpFileName+"/"+outFileName,mode);
    //            break;
    //        case _delete:
    //            inMD1.read(inFileName1);
    //            FOR_ALL_OBJECTS_IN_METADATA(inMD1)
    //            {
    //                inMD1.getValue(MDL::str2Label(_label),inFnImg, __iter.objId);
    //                remove(inFnImg.c_str());
    //                std::cerr << "Remove file: " << inFnImg <<std::endl;
    //            }
    //            break;
    //        case _select:
    //            {
    //                inMD1.read(inFileName1);
    //                outMD.importObjects(inMD1, MDExpression(expression));
    //                outMD.write(outFileName,mode);
    //            }
    //            break;
    //        case _count:
    //            {
    //                inMD1.read(inFileName1);
    //                outMD.aggregate(inMD1, AGGR_COUNT,MDL::str2Label(_label),MDL_CTFMODEL,MDL_COUNT);
    //                outMD.write(outFileName,mode);
    //            }
    //            break;
    //        case _size:
    //            {
    //                inMD1.read(inFileName1);
    //                std::cout << inFileName1 + " size is: " << inMD1.size() << std::endl;
    //            }
    //            break;
    //        case _convert2db:
    //            {
    //                inMD1.read(inFileName1);
    //                MDSql::dumpToFile(outFileName);
    //            }
    //            break;
    //
    //        case _randValues:
    //            {
    //                inMD1.read(inFileName1);
    //                randomize_random_generator();
    //                outMD.randomizeDoubleValues(inMD1,MDL::str2Label(_label), rand_op1, rand_op2, randMode, rand_op3);
    //                outMD.write(outFileName);
    //            }
    //            break;
    //
    //
    //        default:
    //            REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown operation.");
    //        }
    //    }
}
;

int main(int argc, char **argv)
{
    ProgMetadataUtilities program;
    program.read(argc, argv);
    return program.tryRun();
}
