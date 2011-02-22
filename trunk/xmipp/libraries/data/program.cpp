/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (josem@cnb.csic.es)
 *
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

#include "program.h"
#include <stdlib.h>

void XmippProgram::init()
{
    progDef = new ProgramDef();
    this->defineParams();
    ///Add some common definitions to all Xmipp programs
    addParamsLine("== Common options ==");
    addParamsLine("[-v+ <verbose_level=1>] : Verbosity level, 0 means no output.");
    addParamsLine("alias --verbose;");
    addParamsLine("[-h+* <param=\"\">]      : If not param is supplied show this help message.");
    addParamsLine("                        : Otherwise, specific param help is showed,");
    addParamsLine("                        : param should be provided without the '-'");
    addParamsLine("alias --help;");
    addParamsLine("[--gui*]                 : Show a GUI to launch the program.");
    addParamsLine("[--more*]         : Show additional options.");

    ///This are a set of internal command for MetaProgram usage
    ///they should be hidden
    addParamsLine("==+++++ Internal section ==");
    addParamsLine("[--xmipp_write_definition* <dbname>] : Print metadata info about the program to sqlite database");
    addParamsLine("[--xmipp_write_wiki* ] : Print metadata info about the program in wiki format");

    progDef->parse();
}

bool XmippProgram::checkBuiltIns()
{
    ///If -more_options provided, show extended usage
    if (checkParam("--more"))
    {
        usage(1);
        return true;
    }
    ///If help requested, print usage message
    if (checkParam("--help"))
    {
        std::string helpParam = getParam("-h");
        if (helpParam != "")
        {
            std::string cmdHelp("-");
            cmdHelp += helpParam;
            if (existsParam(cmdHelp.c_str()))
                usage(cmdHelp);
            else
            {
                cmdHelp.insert(0, "-");
                if (existsParam(cmdHelp.c_str()))
                    usage(cmdHelp);
                else
                {
                    std::cerr << "Unrecognized param " << helpParam << " neither - or --" << std::endl;
                    usage();
                }
            }
        }
        else
            usage();
        return true;
    }
    if (checkParam("--xmipp_write_definition"))
    {
        writeToDB("programs.db");
        return true;
    }
    if (checkParam("--xmipp_write_wiki"))
    {
        createWiki();
        return true;
    }
    if (checkParam("--gui"))
    {
        createGUI();
        return true;
    }
    return false;
}

void XmippProgram::writeToDB(const FileName &dbName)
{
    XmippDB db;
    DbProgram progData;
    progData.name = name();
    progData.keywords = progDef->keywords;
    StringVector::const_iterator it;
    StringVector & desc = progDef->usageComments.comments;
    for (it = desc.begin(); it < desc.end(); ++it)
        progData.description += *it + "\n";
    db.beginTrans();
    db.insertProgram(&progData);
    db.commitTrans();
}

void XmippProgram::createGUI()
{
    TkPrinter * tk = new TkPrinter();
    tk->printProgram(*progDef);
    delete tk;
}

void XmippProgram::createWiki()
{
    WikiPrinter * wiki = new WikiPrinter();
    wiki->printProgram(*progDef, 3);
    delete wiki;
}

XmippProgram::XmippProgram()
{
    progDef = NULL;
    notRun = true;
    errorCode = 0;
}

XmippProgram::XmippProgram(int argc, char ** argv)
{
    notRun = true;
    errorCode = 0;
    init();
    read(argc, argv);
}

XmippProgram::~XmippProgram()
{
    delete progDef;
}

void XmippProgram::defineParams()
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "function 'defineParams'");
}

void XmippProgram::run()
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "function 'run'");
}

void XmippProgram::quit(int exit_code) const
{
    exit(exit_code);
}

void XmippProgram::readParams()
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "function 'readParams'");
}

void XmippProgram::read(int argc, char ** argv, bool reportErrors)
{
    if (progDef == NULL)
        init();

    setProgramName(argv[0]);

    notRun = true;
    errorCode = 0; //suppose no errors
    ///If not arguments are provided, show the GUI or console program help
    //this behavior will be defined with environment variable XMIPP_BEHAVIOR
    if (argc == 1)
    {
        char * var = getenv("XMIPP_GUI_ON");
        if (var != NULL)
            createGUI();
        else
            usage();
    }
    else
    {
        try
        {
            this->argc = argc;
            this->argv = argv;
            progDef->read(argc, argv, reportErrors);
            if (!checkBuiltIns())
            {
                verbose = getIntParam("--verbose");
                this->readParams();
                notRun = false;
            }
        }
        catch (XmippError xe)
        {
            ///If an input error, shows error message
            std::cerr << xe;
            std::cerr << "For more info use --help" << std::endl;
            errorCode = xe.__errno;
        }
    }
}

void XmippProgram::read(const String &argumentsLine)
{
    int argc;
    char ** argv=NULL;
    char * copy=NULL;

    generateCommandLine(argumentsLine, argc, argv, copy);
    read(argc, argv);

}

int XmippProgram::tryRun()
{
    try
    {
        if (!notRun)
            this->run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
        errorCode = xe.__errno;
    }
    return errorCode;
}

void XmippProgram::setProgramName(const char * name)
{
    progDef->name = name;
}

void XmippProgram::addUsageLine(const char * line, bool verbatim)
{
    progDef->usageComments.addComment(line,verbatim?1:0);
}
void XmippProgram::addExampleLine(const char * example, bool verbatim)
{
    progDef->examples.addComment(example, verbatim ? 1 : 0);
}
void XmippProgram::addSeeAlsoLine(const char * seeAlso)
{
    progDef->seeAlso = seeAlso;
}

void XmippProgram::clearUsage()
{
    progDef->usageComments.clear();
}

void XmippProgram::addParamsLine(const char * line)
{
    progDef->pLexer->addLine((std::string)line);
}

void XmippProgram::addImageFormatParam()
{
    addParamsLine(" -i <input_file>   : Input file: metadata, stack, volume or image.");
    addParamsLine("         :++ Supported read formats are:");
    addParamsLine("         :++ dm3 : Digital Micrograph 3.");
    addParamsLine("         :++ img : Imagic.");
    addParamsLine("         :++ inf,raw : RAW file with header INF file.");
    addParamsLine("         :++ mrc : CCP4.");
    addParamsLine("         :++ spe : Princeton Instruments CCD camera.");
    addParamsLine("         :++ spi, xmp : Spider.");
    addParamsLine("         :++ tif : TIFF.");
    addParamsLine("         :++ ser : tecnai imaging and analysis.");
    addParamsLine("         :++ raw#xDim,yDim,[zDim],offset,datatype,[r] : RAW image file without header file.");
    addParamsLine("         :++ where datatype can be: uint8,int8,uint16,int16,uint32,int32,long,float,double,cint16,cint32,cfloat,cdouble,bool");
    addParamsLine(" alias --input;");
}

void XmippProgram::addWhereImageFormat(const char * whereName)
{
    char whereLine[256];
    sprintf(whereLine, "       where <%s>", whereName);
    addParamsLine(whereLine);
    addParamsLine("         img : Imagic (Data types: uint8, int16, float* and cfloat).");
    addParamsLine("         inf : RAW file with header INF file (Data types: (u)int8, (u)int16 and float*).");
    addParamsLine("         raw : RAW file with header INF file (Data types: (u)int8, (u)int16 and float*).");
    addParamsLine("         mrc : CCP4 (Data types: int8, float* and cfloat).");
    addParamsLine("         spi : Spider (Data types: float* and cfloat).");
    addParamsLine("         xmp : Spider (Data types: float* and cfloat).");
    addParamsLine("         tif : TIFF. (Data types: uint8*, uint16, uint32 and float).");
    addParamsLine("         custom <ext> : Custom extension name, the real format will be Spider.");
}

void XmippProgram::addWhereRandomType(const char * randomName)
{
    char randomLine[256];
    sprintf(randomLine, "       where <%s>", randomName);
    addParamsLine(randomLine);
    addParamsLine("gaussian <stddev> <avg=0.>        :Gaussian distribution parameters");
    addParamsLine("student <df> <stddev> <avg=0.> :t-student distribution parameters");
    addParamsLine("uniform  <min> <max>           :Uniform distribution parameters");
}

void XmippProgram::addKeywords(const char * keywords)
{
    progDef->keywords += " ";
    progDef->keywords += keywords;
}

const char * XmippProgram::getParam(const char * param, int arg)
{
    return progDef->getParam(param, arg);
}

const char * XmippProgram::getParam(const char * param, const char * subparam, int arg)
{
    return progDef->getParam(param, subparam, arg);
}

int XmippProgram::getIntParam(const char * param, int arg)
{
    return textToInteger(progDef->getParam(param, arg));
}

int XmippProgram::getIntParam(const char * param, const char * subparam, int arg)
{
    return textToInteger(progDef->getParam(param, subparam, arg));
}

double XmippProgram::getDoubleParam(const char * param, int arg)
{
    return textToFloat(progDef->getParam(param, arg));
}

double XmippProgram::getDoubleParam(const char * param, const char * subparam, int arg)
{
    return textToFloat(progDef->getParam(param, subparam, arg));
}

void XmippProgram::getListParam(const char * param, StringVector &list)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((std::string)"Doesn't exists param: " + param));
    list.clear();
    for (int i = 0; i < paramDef->cmdArguments.size(); ++i)
        list.push_back(paramDef->cmdArguments[i]);
}

int XmippProgram::getCountParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((std::string)"Doesn't exists param: " + param));
    return paramDef->cmdArguments.size();
}

bool XmippProgram::checkParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((std::string)"Doesn't exists param: " + param));
    return paramDef->counter == 1;
}

bool XmippProgram::existsParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    return paramDef != NULL;
}



const char * XmippProgram::name() const
{
    return progDef->name.c_str();
}

void XmippProgram::usage(int verb) const
{
    ConsolePrinter cp;
    cp.printProgram(*progDef, verb);
}

void XmippProgram::usage(const std::string & param, int verb)
{
    ConsolePrinter cp;
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((std::string)"Doesn't exists param: " + param));
    cp.printParam(*paramDef, verb);
    quit(0);
}

void XmippProgram::show() const
{}

int XmippProgram::version() const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"");
}

int XmippProgram::runProgram(XmippProgram * program, const String &arguments, bool destroy)
{
    program->read(arguments);
    int retCode = program->tryRun();
    if (destroy)
        delete program;
    return retCode;
}

/// Empty constructor
XmippMetadataProgram::XmippMetadataProgram()
{
    oroot = oext = fn_out = fn_in = "";
    produces_an_output = false;
    each_image_produces_an_output = false;
    allow_time_bar = true;
    apply_geo = false;
    decompose_stacks = true;
    save_metadata_stack = false;
}

void XmippMetadataProgram::defineParams()
{
    addImageFormatParam();
    addParamsLine(" [--mode+ <mode=overwrite>]   : Metadata writing mode.");
    addParamsLine("    where <mode>");
    addParamsLine("     overwrite   : Replace the content of the file with the Metadata");
    addParamsLine("     append      : Write the Metadata as a new block, removing the old one");

    if (each_image_produces_an_output)
    {
        addParamsLine("  [-o <output_file=\"\">]  : Output file: metadata, stack, volume or image.");
        addParamsLine("   alias --output;");
        //      DO NOT LONGER SUPPORT --OEXT SINCE NOW IT DOES MIND CHANGE FORMAT
        //        addParamsLine("  [--oext <extension=\"\">] :  Output file format extension.");
        //        addExtensionWhere("extension");
        addParamsLine("  [--oroot <root=\"\">]     : Rootname of output individual images.");
        addParamsLine("                            : Output image format can be set adding extension after rootname as \":ext\".");
    }
    else if (produces_an_output)
    {
        addParamsLine("  [-o <output_file=\"\">]  : Output file: metadata, stack, volume or image.");
        addParamsLine("   alias --output;");
    }

    if (apply_geo)
    {
        addParamsLine("  [--dont_apply_geo]   : for 2D-images: do not apply transformation stored in the header");
    }
}

void XmippMetadataProgram::readParams()
{
    fn_in = getParam("-i");
    mode = metadataModeConvert(getParam("--mode"));

    if (produces_an_output)
        fn_out = checkParam("-o") ? getParam("-o") : "";

    if (each_image_produces_an_output)
    {
        fn_out = checkParam("-o") ? getParam("-o") : "";
        oroot = getParam("--oroot");
    }

    if (fn_out != fn_in && oroot == "")
    {
        FileName fn_stack_plain = fn_out.removeFileFormat();
        if (exists(fn_stack_plain))
            unlink(fn_stack_plain.c_str());
    }

    mdIn.read(fn_in, NULL, decompose_stacks);

    single_image = input_is_stack = false;
    if (!fn_in.isMetaData())
    {
        if (mdIn.size() == 1)
            single_image = true;
        else
            input_is_stack = true;
    }
    single_image = !fn_in.isMetaData() && (mdIn.size() == 1);

    if (mdIn.containsLabel(MDL_ENABLED))
        mdIn.removeObjects(MDValueEQ(MDL_ENABLED, -1));

    if (mdIn.isEmpty())
        REPORT_ERROR(ERR_MD_NOOBJ, "");

    if (apply_geo)
        apply_geo = !checkParam("--dont_apply_geo");

}

void XmippMetadataProgram::show()
{
    if (verbose==0)
        return;
    std::cout << "Input File: " << fn_in << std::endl;
    if (apply_geo)
        std::cout << "Applying transformation stored in header of 2D-image" << std::endl;
    if (each_image_produces_an_output)
    {
        if (fn_out != "")
            std::cout << "Output File: " << fn_out << std::endl;
        //      DO NOT LONGER SUPPORT --OEXT SINCE NOW IT DOES MIND CHANGE FORMAT
        //if (oext != "")
        //    std::cout << "Output Extension: " << oext << std::endl;
        if (oroot != "")
            std::cout << "Output Root: " << oroot << std::endl;
    }
}

void XmippMetadataProgram::preProcess()
{}

void XmippMetadataProgram::postProcess()
{}

void XmippMetadataProgram::startProcessing()
{
    //Show some info
    show();
    // Initialize progress bar
    time_bar_size = mdIn.size();
    if (allow_time_bar && verbose && !single_image)
        init_progress_bar(time_bar_size);
    time_bar_step = CEIL((double)time_bar_size / 60.0);
    time_bar_done = 0;
}

void XmippMetadataProgram::finishProcessing()
{
    if (allow_time_bar && verbose && !single_image)
        progress_bar(time_bar_size);

    if (!single_image && !mdOut.isEmpty())
    {
        if (oroot != "") // Out as independent images
            mdOut.write(fn_out);
        else if (save_metadata_stack && fn_out != "") // Out as stack
            mdOut.write(fn_out.withoutExtension().addExtension("sel"));
    }
}

void XmippMetadataProgram::showProgress()
{
    if (time_bar_done % time_bar_step == 0 && allow_time_bar && verbose && !single_image)
        progress_bar(time_bar_done);
}

size_t XmippMetadataProgram::getImageToProcess()
{
    if (time_bar_done == 0)
        iter = new MDIterator(mdIn);
    else
        iter->moveNext();

    ++time_bar_done;
    return iter->objId;
}

void XmippMetadataProgram::run()
{
    try
    {
        FileName fnImg, fnImgOut, baseName, pathBaseName, fullBaseName, oextBaseName;
        size_t objId;
        //Perform particular preprocessing
        preProcess();

        startProcessing();

        int kk = 0;

        if (oroot != "" )
        {
            if (oext == "")
                oext           = oroot.getFileFormat();
            oextBaseName   = oext;
            fullBaseName   = oroot.removeFileFormat();
            baseName       = fullBaseName.getBaseName();
            pathBaseName   = fullBaseName.getRoot();
        }

        //FOR_ALL_OBJECTS_IN_METADATA(mdIn)
        while ((objId = getImageToProcess()) != BAD_OBJID)
        {
            mdIn.getValue(MDL_IMAGE, fnImg, objId);

            if (fnImg == "")
                break;

            fnImgOut = fnImg;

            if (each_image_produces_an_output)
            {
                if (oroot != "" ) // Compose out name to save as independent images
                {
                    if (oext == "")
                        oextBaseName = fnImg.getFileFormat();

                    if (baseName != "")
                        fnImgOut.compose(fullBaseName, ++kk, oextBaseName);
                    else if (fnImg.isInStack())
                        fnImgOut.compose(pathBaseName + (fnImg.withoutExtension()).getDecomposedFileName(), ++kk, oextBaseName);
                    else
                        fnImgOut = pathBaseName + fnImg.withoutExtension()+ "." + oextBaseName;
                }
                else if (fn_out != "")
                {
                    if (single_image)
                        fnImgOut = fn_out;
                    else
                        fnImgOut.compose(++kk, fn_out); // Compose out name to save as stacks
                }
                else
                    fnImgOut = fnImg;

                newId = mdOut.addObject();
                mdOut.setValue(MDL_IMAGE, fnImgOut, newId);
                mdOut.setValue(MDL_ENABLED, 1, newId);
            }

            processImage(fnImg, fnImgOut, objId);

            showProgress();
        }

        // Generate name to save mdOut when output are independent images
        if (oroot != "" && fn_out == "")
        {
            if (baseName != "")
                fn_out = baseName.addExtension("sel");
            else
                fn_out = fn_in.withoutExtension() + "_" + oextBaseName + ".sel";
        }

        finishProcessing();

        postProcess();
    }
    catch (XmippError xe)
    {
        std::cout << xe;
        quit(xe.__errno);
    }
}


