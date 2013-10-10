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
#include <data/xmipp_program.h>
#include <string.h>
#include <data/metadata.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <data/xmipp_funcs.h>

class ProgMetadataUtilities: public XmippProgram
{
private:
    WriteModeMetaData mode;
    FileName fn_in, fn_out, fn_md2;
    MetaData mdIn, md2;
    MDLabel label;
    std::vector<MDLabel> labels;
    String operation, order;
    bool doWrite;

protected:
    void defineParams()
    {
        addUsageLine("Perform several operations on metadata files. ");
        addUsageLine("If the -o option is not used the original metadata will be modified.");
        addUsageLine("+ Also you can use the --print option just to print out the result metadata to screen.");
        addUsageLine("+ The combination of -i and -o without other operations can serve to extract data blocks");
        addUsageLine("+ inside a medata and write to an independent one.");
        addSeeAlsoLine("metadata_import");

        addParamsLine(" -i <metadata>         : Input metadata file");
        addParamsLine("   [-o  <metadata>]    : Output metadata file, if not provided result will overwrite input file");

        addParamsLine("  [--set <set_operation> <md2_file> <label=image>]   : Set operations");
        addParamsLine("         where <set_operation>");
        addParamsLine("   union               : Union with metadata md2, duplicated values only will appear once");
        addParamsLine("   union_all           : Union with metadata md2, will repeat duplicated values");
        addParamsLine("   intersection        : Intersection with metadata md2");
        addParamsLine("   subtraction         : Subtraction with metadata md2");
        addParamsLine("   join                : Inner join with md2 using label l1");
        addParamsLine("   natural_join        : Natural  join with md2 using all common labels");
        addParamsLine("   merge               : Merge columns with md2, label is ignored");
        addParamsLine("                       : Both metadatas should have same size, and elements should be in same order,");
        addParamsLine("                       : if not, you should use 'join' instead, but this constrain having a common label");
        addParamsLine("           alias -s;                                             ");

        addParamsLine("or --operate <operation>     : Operations on the metadata structure");
        addParamsLine("         where <operation>");
        addParamsLine("   sort <label=image> <order=asc>: Sort metadata using a label as identifier");
        addParamsLine("                                 : for sorting according to a component of a vector label");
        addParamsLine("                                 : use label:col, e.g., NMADisplacements:0");
        addParamsLine("                                 : The first column is column number 0.");
        addParamsLine("                                 : order can be asc (ascending) or desc (descending)");
        addParamsLine("    random_subset <size>                : Extract a random subset without replacement of this metadata");
        addParamsLine("    bootstrap                           : Extract a bootstrap subset (with replacement) of this metadata");
        addParamsLine("    randomize                           : Randomize elements of metadata");
        addParamsLine("    keep_column <labels>                : Keep some columns(label list) from metadata");
        addParamsLine("    drop_column <labels>                : Drop some columns(label list) from metadata");
        addParamsLine("    rename_column <labels>              : Rename a column");
        addParamsLine("    modify_values <expression>          : Use an SQLite expression to modify the metadata");
        addParamsLine("                                        : This option requires knowledge of basic SQL syntax(more specific SQLite");
        addParamsLine("    expand <factor>                     : Expand the metadata content by union with himself ");
        addParamsLine(":+ Some of the function allowed are,");
        addParamsLine(":+ Math: +, -, *, /, abs");
        addParamsLine(":+       acos, asin, atan, atn2, atan2, acosh, asinh, atanh,");
        addParamsLine(":+       difference,");
        addParamsLine(":+       degrees, radians, cos, sin, tan, cot, cosh, sinh, tanh,");
        addParamsLine(":+       coth, exp,");
        addParamsLine(":+       log, log10, power, sign, sqrt, square, ceil, floor, pi.");
        addParamsLine(":+ String: replicate, charindex, leftstr, rightstr, ltrim, rtrim, trim,");
        addParamsLine(":+         replace, reverse, proper, padl, padr, padc, strfilter.");

        addParamsLine(":+ Aggregate: max, min, avg, sum, mstdev, variance, mode, median,");
        addParamsLine(":+            lower_quartile, upper_quartile.");
        addParamsLine("           alias -e;                                             ");


        addParamsLine("or  --file <file_operation>     : File operations");
        addParamsLine("         where <file_operation>");
        addParamsLine("   copy <directory> <label=image>  : Copy files in metadata md1 to directory path (file names at label column)");
        addParamsLine("   move <directory>  <label=image> : Move files in metadata md1 to directory path (file names at label column)");
        addParamsLine("   delete  <label=image>      : Delete files in metadata md1 (file names at label column)");
        //        addParamsLine("   convert2db                 : Convert metadata to sqlite database");
        //        addParamsLine("   convert2xml                : Convert metadata to xml file");
        addParamsLine("   import_txt <labels>        : Import a text file specifying its columns");
        addParamsLine("           alias -f;                                             ");

        addParamsLine("or --query <query_operation>   : Query operations");
        addParamsLine("         where <query_operation>");
        addParamsLine("   select <expression>        : Create new metadata with those entries that satisfy the expression");
        addParamsLine("   count  <label>             : for each value of a given label create new metadata with the number of times the value appears");
        addParamsLine("   sum  <label1> <label2>   : group metadata by label1 and add quantities in label2");
        addParamsLine("   size                       : print Metadata size");
        addParamsLine("   labels                     : print Metadata labels");
        addParamsLine("   blocks                     : print blocks in file");
        addParamsLine("           alias -q;                                             ");

        addParamsLine("or --fill <labels> <fill_mode>                  : Fill a column values(should be of same type)");
        addParamsLine("   where <fill_mode>");
        addParamsLine("     constant  <value>                        : Fill with a constant value");
        addParamsLine("     lineal  <init_value> <step>              : Fill with a lineal serie starting at init_value with an step");
        addParamsLine("     rand_uniform  <a=0.> <b=1.>              : Follow a uniform distribution between a and b");
        addParamsLine("     rand_gaussian <mean=0.> <stddev=1.>      : Follow a gaussian distribution with mean and stddev");
        addParamsLine("     rand_student  <mean=0.> <stddev=1.> <df=3.> : Follow a student distribution with mean, stddev and df degrees of freedom.");
        addParamsLine("     expand                                   : Treat the column as the filename of an row metadata and expand values");
        addParamsLine("           alias -l;                                             ");

        addParamsLine("  [--print ]                    : Just print medata to stdout, or if -o is specified written to disk.");
        addParamsLine("                                : this option is useful for extrating data blocks inside a metadata.");
        addParamsLine("           alias -p;                                             ");

        addParamsLine(" [--mode <mode=overwrite>]   : Metadata writing mode.");
        addParamsLine("    where <mode>");
        addParamsLine("     overwrite   : Replace the content of the file with the Metadata");
        addParamsLine("     append      : Write the Metadata as a new block, removing the old one");

        addExampleLine(" Concatenate two metadatas. If label is not provided, by default is 'image'", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --set union mD2.doc  -o out.doc");
        addExampleLine(" Intersect two metadatas using label 'order_'", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --set intersection mD2.doc order_ -o out.doc");
        addExampleLine(" Combine columns from two metadatas. Be sure of both have same number of rows and also", false);
        addExampleLine(" there aren't common columns, in that case second metadata columns will be used", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --set merge mD2.doc -o out.doc");
        addExampleLine(" Sort the elements in metadata (using default label 'image').", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc -s sort -o out.doc");
        addExampleLine(" You can also add columns and 'filling' its values with different options", false);
        addExampleLine("By example, to add the column 'shiftX' with uniform random value between 0 and 10", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --fill shiftX rand_uniform 0 10 -o out.doc");
        addExampleLine("Or for initialize metadata columns 'shiftX' and 'shiftY' with a constant value of 5", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc -l \"shiftX shiftY\" constant 5 -o out.doc");
        addExampleLine("If you have columns that represent the filename of a metadata with other data (ex CTFParams)", false);
        addExampleLine("you cant 'expand' the column with the values in that metadta", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --fill CTFParams expand -o outExpanded.doc");
        addExampleLine("For check all options availables for 'filling' mode, use: ", false);
        addExampleLine ("   xmipp_metadata_utilities --help fill");
        addExampleLine(" write metadata as table in Sqlite3 database. (use xmipp_sqlite3 to visualize results)", false);
        addExampleLine ("   xmipp_metadata_utilities -i blocknameIn@mD1.doc -o blocknameOut@mD1.sqlite");
        addExampleLine(" write metadata as xml file.", false);
        addExampleLine ("   xmipp_metadata_utilities -i blocknameIn@mD1.doc -o blocknameOut@mD1.xml");
        addExampleLine(" Copy files in metadata to a location. The metadata will be also copied to new location", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --file copy /home/pepe/newLocation");
        addExampleLine(" Delete files in metadata.", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --file delete");
        addExampleLine(" Select elements in metadata that satisfy a given constrain.", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --query select \"angleRot > 10 AND anglePsi < 0.5\" -o out.doc");
        addExampleLine(" You can also modify your data using SQLite syntax expression", false);
        addExampleLine("  xmipp_metadata_utilities  -i a.doc --operate modify_values \"angleRot=2.*angleRot\" -o b.doc");
        addExampleLine("  xmipp_metadata_utilities  -i a.doc --operate modify_values \"angleRot=radians(angleRot)\" -o b.doc");
        addExampleLine("  xmipp_metadata_utilities  -i a.doc --operate modify_values \"angleRot=sin(radians(angleRot))\" -o b.doc");
        addExampleLine("  xmipp_metadata_utilities  -i a.doc --operate modify_values \"angleRot=sqrt(angleRot)\" -o b.doc");
        addExampleLine("  xmipp_metadata_utilities  -i a.doc --operate modify_values \"image=replace(image, 'xmp','spi')\" -o b.doc");
        addExampleLine("  xmipp_metadata_utilities  -i a.doc --operate modify_values \"image='new_prefix_dir/'||image\" -o b.doc");
        addExampleLine(" Count number of images per CTF", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc -q count CTFModel -o out.doc");
        addExampleLine(" images asigned a ctfgroup", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc -q sum defocusGroup count -o out.doc");
        addExampleLine(" Print the metadata Size", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --query size");
        addExampleLine(" Rename Column", false);
        addExampleLine ("   xmipp_metadata_utilities -i mD1.doc --operate rename_column \"weight wRobust\"");

    }

    void readParams()
    {
        fn_in = getParam("-i");
        // Prevent from reading the input metadata for some special cases
        bool readIn = !((checkParam("--file") && STR_EQUAL(getParam("--file"), "import_txt")) || //when importing from .txt files
                        (checkParam("--query") && STR_EQUAL(getParam("--query"), "blocks")));

        if (checkParam("--query") &&
            (STR_EQUAL(getParam("--query"), "size") ||
             STR_EQUAL(getParam("--query"), "labels")))
          mdIn.setMaxRows(1); //avoid parse the entire metadata

        if (readIn)
          mdIn.read(fn_in);
        doWrite = true;
        fn_out = checkParam("-o") ? getParam("-o") : fn_in;
        mode = MD_OVERWRITE;
        if (checkParam("--mode")
            && STR_EQUAL(getParam("--mode"), "append"))
            mode = MD_APPEND;
    }

    void doSet()
    {
        operation = getParam("--set", 0);
        md2.read(getParam("--set", 1));
        MDLabel label = MDL::str2Label(getParam("--set", 2));

        if (operation == "union")
            mdIn.unionDistinct(md2, label);
        else if (operation == "union_all")
            mdIn.unionAll(md2);
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
        else if (operation == "natural_join")
        {
            MetaData md;
            md.join(mdIn, md2, label,NATURAL);
            mdIn = md;
        }
        else if (operation == "merge")
            mdIn.merge(md2);
    }//end of function doSet

    void doOperate()
    {
        operation = getParam("--operate", 0);

        if ( operation == "keep_column")
        {
            MDL::str2LabelVector(getParam("--operate", 1), labels);
            mdIn.keepLabels(labels);
            std::cout << labels.size() << std::endl;
            //for (int i = 0; i < labels.size(); ++i)
            //    mdIn.addLabel(labels[i]);//removeLabel(labels[i]);
        }
        else if ( operation == "drop_column")
               {
                   MDL::str2LabelVector(getParam("--operate", 1), labels);
                   for (size_t i = 0; i < labels.size(); ++i)
                       mdIn.removeLabel(labels[i]);
               }
        else if ( operation == "rename_column")
               {
                   MDL::str2LabelVector(getParam("--operate", 1), labels);
				   mdIn.renameColumn(labels[0],labels[1]);
               }
        else if (operation == "modify_values")// modify_values
        {
            MDSql::activateMathExtensions();
            mdIn.operate(getParam("--operate", 1));
        }
        else if (operation == "expand")// modify_values
        {
          int factor = getIntParam("--operate", 1);
          MetaData md;
          for (int i = 0; i < factor; i++)
            md.unionAll(mdIn);

          mdIn = md;
        }else
        {
            MetaData md(mdIn);
            if (operation == "sort")
            {
            	String order=getParam("--operate",2);
                mdIn.sort(md, getParam("--operate", 1),order=="asc");
            }
            else if (operation == "randomize")
                mdIn.randomize(md);
            else if (operation == "bootstrap")
            {
            	std::vector<size_t> objId;
            	objId.resize(md.size());
            	size_t n=0;
            	FOR_ALL_OBJECTS_IN_METADATA(md)
            	objId[n++]=__iter.objId;
            	// md.getColumnValues(MDL_OBJID,objId); COSS: It should work, but it does not
            	int N_1=((int)objId.size())-1;
            	MDRow row;
            	MetaData mdAux;
            	FOR_ALL_OBJECTS_IN_METADATA(md)
            	{
            		md.getRow(row,objId[(size_t)rnd_unif(0,N_1)]);
            		mdAux.setRow(row,mdAux.addObject());
            	}
            	mdIn.sort(mdAux,MDL_IMAGE);
            }
            else if (operation == "random_subset")
            {
            	MetaData mdAux, mdAux2;
                mdAux.randomize(md);
                md.clear();
                mdAux2.selectPart(mdAux, 0, getIntParam("--operate", 1));
                mdAux.clear();
                mdIn.sort(mdAux2,MDL_IMAGE);
            }
        }
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
        else if (operation == "lineal")
            generator = new MDLinealGenerator(getDoubleParam("--fill", 2), getDoubleParam("--fill", 3));

        //Fill columns
        for (size_t i = 0; i < labels.size(); ++i)
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

        if (operation == "count")//note second label is a dummy parameter
        {
            label = MDL::str2Label(getParam("--query", 1));
            md2 = mdIn;
            mdIn.aggregate(md2, AGGR_COUNT,label,label,MDL_COUNT);
        }
        else if (operation == "sum")
        {
            label = MDL::str2Label(getParam("--query", 1));
            MDLabel label2;
            label2 = MDL::str2Label(getParam("--query", 2));
            md2 = mdIn;
            mdIn.aggregate(md2, AGGR_SUM,label,label2,MDL_SUM);
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
            std::cout << fn_in + " size is: " << mdIn.getParsedLines() << std::endl;
        }
        else if (operation == "labels")
        {
            doWrite = false;
            std::cout << fn_in + " has labels: " << std::endl;
            MDLabelVector labels = mdIn.getActiveLabels();
            for (size_t i = 0; i < labels.size(); ++i)
              std::cout << "  " << MDL::label2Str(labels[i]) << std::endl;
        }
        else if (operation == "blocks")
        {
            doWrite = false;
            StringVector blocks;
            std::cout << "Blocks in " << fn_in << ": " << std::endl;
            getBlocksInMetaDataFile(fn_in, blocks);
            for (size_t i = 0; i < blocks.size(); ++i)
                std::cout << blocks[i] << std::endl;
        }
    }//end of function doQuery

    void doFile()
    {
        operation = getParam("--file", 0);

        //        if (operation == "convert2db")
        //        {
        //            doWrite = false;
        //            MDSql::dumpToFile(fn_out);
        //        }
        //        else if (operation == "convert2xml")
        //        {
        //         fn_out.re
        //            doWrite = false;
        //            mdIn.writeXML(fn_out);
        //        }
        if (operation == "import_txt")
        {
            mdIn.readPlain(fn_in, getParam("--file", 1));
        }
        else
        {
            bool doDelete;
            FileName path, inFnImg, outFnImg, oldOutFnImg="";
            if (!(doDelete = operation == "delete"))//copy or move
            {
                path = getParam("--file", 1);
                if (!path.exists())
                    if (path.makePath() != 0)
                        REPORT_ERROR(ERR_IO_NOPERM, (String)"Cannot create directory "+ path);
            }
            label = MDL::str2Label(getParam("--file", doDelete ? 1 : 2));
            doWrite = !doDelete;

            int counter=FIRST_IMAGE;
            FOR_ALL_OBJECTS_IN_METADATA(mdIn)
            {

                mdIn.getValue(label, inFnImg, __iter.objId);
                bool isStack=inFnImg.isInStack();
                if (doDelete)
                {
                    if(isStack)
                        REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Cannot delete files from a stack");
                    else
                        remove(inFnImg.c_str());
                }
                else
                {
                    outFnImg = inFnImg.removeDirectories();
                    //outfilename for output sel
                    if (operation == "copy" && isStack)
                    {
                        if(oldOutFnImg != outFnImg && oldOutFnImg !="")
                        {
                            counter = FIRST_IMAGE;
                        }
                        oldOutFnImg = outFnImg;
                        outFnImg.compose(counter, outFnImg);
                    }
                    mdIn.setValue(label, outFnImg, __iter.objId);
                    //output file name for copy
                    if (operation == "copy" && isStack)
                    {
                        outFnImg=outFnImg.removePrefixNumber();
                        outFnImg = path + "/" + outFnImg;
                        outFnImg = integerToString(counter)+'@'+outFnImg;
                    }
                    else
                        outFnImg = path + "/" + outFnImg;
                    if (operation == "copy")
                    {
                        copyImage( inFnImg, outFnImg);
                        ++counter;
                    }
                    else if (operation == "move")
                    {
                        if(isStack)
                            REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Cannot move files from a stack");
                        else
                            rename(inFnImg.c_str(), outFnImg.c_str());
                    }
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

        if (checkParam("--print"))
            mdIn.write(std::cout);
        else if (doWrite)
            mdIn.write(fn_out, mode);
    }

}
;//end of class ProgMetaDataUtilities

RUN_XMIPP_PROGRAM(ProgMetadataUtilities)
