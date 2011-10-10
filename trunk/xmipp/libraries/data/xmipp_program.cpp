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

#include "xmipp_program.h"
#include "metadata_extension.h"
#include <stdlib.h>

void XmippProgram::initComments()
{
    CommentList comments;
    comments.addComment("Verbosity level, 0 means no output.");
    defaultComments["-v"] = comments;
}

void XmippProgram::processDefaultComment(const char *param, const char *left)
{
    addParamsLine(((String)left+":"+defaultComments[param].comments[0]).c_str());
    int imax=defaultComments[param].comments.size();
    for (int i=1; i<imax; ++i)
        addParamsLine(((String)":"+defaultComments[param].comments[i]).c_str());
}

void XmippProgram::defineCommons()
{
    ///Add some common definitions to all Xmipp programs
    addParamsLine("== Common options ==");
    processDefaultComment("-v","[-v+ <verbose_level=1>]");
    addParamsLine("alias --verbose;");
    addParamsLine("[-h+* <param=\"\">]      : If not param is supplied show this help message.");
    addParamsLine("                         : Otherwise, specific param help is showed,");
    addParamsLine("                         : param should be provided without the '-'");
    addParamsLine("alias --help;");
    addParamsLine("[--gui*]                 : Show a GUI to launch the program.");
    addParamsLine("[--more*]                : Show additional options.");

    ///This are a set of internal command for MetaProgram usage
    ///they should be hidden
    addParamsLine("==+++++ Internal section ==");
    addParamsLine("[--xmipp_write_definition* <dbname>] : Print metadata info about the program to sqlite database");
    addParamsLine("[--xmipp_write_wiki* ] : Print metadata info about the program in wiki format");
    addParamsLine("[--xmipp_write_protocol* <scriptfile>] : Generate protocol header file");
    addParamsLine("[--xmipp_run_id <run_id>] : This is only meanful when execute throught protocols");
}

void XmippProgram::init()
{
    initComments();
    progDef = new ProgramDef();
    this->defineParams();
    this->defineCommons();
    progDef->parse();
}

bool XmippProgram::checkBuiltIns()
{
    ///If -more_options provided, show extended usage
    if (checkParam("--more"))
        usage(1);
    ///If help requested, print usage message
    else if (checkParam("--help"))
    {
        String helpParam = getParam("-h");
        if (helpParam != "")
        {
            String cmdHelp("-");
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
                    if (verbose)
                        std::cerr << "Unrecognized param " << helpParam << " neither - or --" << std::endl;
                    usage();
                }
            }
        }
        else
            usage();
    }
    else if (checkParam("--xmipp_write_definition"))
        writeToDB();
    else if (checkParam("--xmipp_write_wiki"))
        createWiki();
    else if (checkParam("--xmipp_write_protocol"))
        writeToProtocol();
    else if (checkParam("--gui"))
        createGUI();
    else
        return false;
    return true;
}

void XmippProgram::writeToDB()
{
    ProgramDb db;
    db.printProgram(*progDef);
}

void XmippProgram::writeToProtocol( )
{
    String scriptfile = getParam("--xmipp_write_protocol");
    ProtPrinter pp(scriptfile.c_str());
    pp.printProgram(*progDef);
}

void XmippProgram::createGUI()
{
    TkPrinter tk;
    tk.printProgram(*progDef);
}

void XmippProgram::createWiki()
{
    WikiPrinter wiki;
    wiki.printProgram(*progDef, 3);
}

XmippProgram::XmippProgram()
{
    //by defaul all programs have verbose = 1
    // this can be changed on mpi slaves node for no output at all
    verbose = 1;
    progDef = NULL;
    runWithoutArgs = doRun = false;
    errorCode = 0;
}

XmippProgram::XmippProgram(int argc, char ** argv)
{
    runWithoutArgs = doRun = false;
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

    doRun = false;
    errorCode = 0; //suppose no errors
    ///If not arguments are provided, show the GUI or console program help
    //this behavior will be defined with environment variable XMIPP_BEHAVIOR
    if (argc == 1)
    {
        if (runWithoutArgs)
            doRun = true;
        else
        {
            char * var = getenv("XMIPP_GUI_ON");
            if (var != NULL)
                createGUI();
            else
                usage();
        }
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
                if (verbose) //if 0, ignore the parameter, useful for mpi programs
                    verbose = getIntParam("--verbose");
                this->readParams();
                doRun = true;
            }
        }
        catch (XmippError &xe)
        {
            ///If an input error, shows error message
            if (verbose)
            {
                std::cerr << xe;
                std::cerr << "For more info use --help" << std::endl;
            }
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
        if (doRun)
            this->run();
    }
    catch (XmippError &xe)
    {
        std::cerr << xe;
        errorCode = xe.__errno;
    }
        std::cout << "Finished run" << std::endl;
    return errorCode;
}

void XmippProgram::setProgramName(const char * name)
{
    progDef->name = name;
}

void XmippProgram::addUsageLine(const char * line, bool verbatim)
{
    progDef->usageComments.addComment(line,verbatim);
}
void XmippProgram::addExampleLine(const char * example, bool verbatim)
{
    progDef->examples.addComment(example, verbatim);
}
void XmippProgram::addSeeAlsoLine(const char * seeAlso)
{
    if (progDef->seeAlso=="")
        progDef->seeAlso = seeAlso;
    else
    {
        progDef->seeAlso +=", ";
        progDef->seeAlso +=seeAlso;
    }
}

void XmippProgram::clearUsage()
{
    progDef->usageComments.clear();
}
void XmippProgram::addParamsLine(const String &line)
{
    progDef->pLexer->addLine(line);
}

void XmippProgram::addParamsLine(const char * line)
{
    progDef->pLexer->addLine((String)line);
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
        REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
    list.clear();
    for (int i = 0; i < paramDef->cmdArguments.size(); ++i)
        list.push_back(paramDef->cmdArguments[i]);
}

int XmippProgram::getCountParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
    return paramDef->cmdArguments.size();
}

bool XmippProgram::checkParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
    return paramDef->counter == 1;
}

bool XmippProgram::existsParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    return paramDef != NULL;
}


ParamDef * XmippProgram::getParamDef(const char * param) const
{
    return progDef->findParam(param);
}

const char * XmippProgram::name() const
{
    return progDef->name.c_str();
}

void XmippProgram::usage(int verb) const
{
    if (verbose)
    {
        ConsolePrinter cp;
        char * var = getenv("XMIPP_COLOR_OFF");
        if (var != NULL)
            cp.color = false;
        cp.printProgram(*progDef, verb);
    }
}

void XmippProgram::usage(const String & param, int verb)
{
    if (verbose)
    {
        ConsolePrinter cp;
        ParamDef * paramDef = progDef->findParam(param);
        if (paramDef == NULL)
            REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
        cp.printParam(*paramDef, verb);
    }
    quit(0);
}

void XmippProgram::show() const
{}

int XmippProgram::version() const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"");
}

/// Empty constructor
XmippMetadataProgram::XmippMetadataProgram()
{
    oroot = oext = fn_out = fn_in = "";
    apply_geo=false;
    allow_apply_geo = false;
    produces_an_output = false;
    each_image_produces_an_output = false;
    allow_time_bar = true;
    decompose_stacks = true;
    save_metadata_stack = false;
    delete_output_stack = true;
    remove_disabled = true;
    single_image = input_is_stack = output_is_stack = false;
    iter = NULL;
}

void XmippMetadataProgram::initComments()
{
    XmippProgram::initComments();

    CommentList comments;
    comments.addComment("Input file: metadata, stack, volume or image.");
    defaultComments["-i"]=comments;

    comments.clear();
    comments.addComment("Output file: metadata, stack, volume or image.");
    defaultComments["-o"]=comments;

    comments.clear();
    comments.addComment("Rootname of output individual images.");
    comments.addComment("Output image format can be set adding extension after rootname as \":ext\".");
    defaultComments["--oroot"]=comments;
}

void XmippMetadataProgram::defineParams()
{
    processDefaultComment("-i","-i <input_file>");
    addParamsLine("alias --input;");
    addParamsLine(" [--mode+ <mode=overwrite>]   : Metadata writing mode.");
    addParamsLine("    where <mode>");
    addParamsLine("     overwrite   : Replace the content of the file with the Metadata");
    addParamsLine("     append      : Write the Metadata as a new block, removing the old one");

    if (each_image_produces_an_output)
    {
        processDefaultComment("-o","[-o <output_file=\"\">]");
        addParamsLine("   alias --output;");
        processDefaultComment("--oroot","[--oroot <root=\"\">]");
    }
    else if (produces_an_output)
    {
        processDefaultComment("-o","[-o <output_file=\"\">]");
        addParamsLine("   alias --output;");
    }

    if (allow_apply_geo)
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

    mdIn.read(fn_in, NULL, decompose_stacks);

    if (mdIn.isEmpty())
        REPORT_ERROR(ERR_MD_NOOBJ, "Empty input Metadata.");

    mdInSize = mdIn.size();

    bool inputIsMetaData = fn_in.isMetaData();

    if (!inputIsMetaData)
    {
        if (mdInSize == 1)
            single_image = true;
        else
            input_is_stack = true;
    }
    /* Output is stack if, given a filename in fn_out, mdIn has multiple images.
     * In case no output name is given, then input is overwritten and we have to
     * check if it is stack. */
    output_is_stack = mdInSize > 1 && oroot.empty() && (!fn_out.empty() || input_is_stack);

    save_metadata_stack = save_metadata_stack && inputIsMetaData && output_is_stack;

    if (remove_disabled)
        mdIn.removeDisabled();

    // if input is volume do not apply geo
    if (allow_apply_geo && zdimOut == 1)
        apply_geo = !checkParam("--dont_apply_geo");

    // Only delete output stack in case we are not overwritting input
    delete_output_stack = (output_is_stack && delete_output_stack) ?
                          !(fn_out.empty() && oroot.empty()) : false;

    // If the output is a stack, create empty stack file in advance to avoid concurrent access to the header
    create_empty_stackfile = (each_image_produces_an_output && output_is_stack && !fn_out.empty());

    // if create, then we need to read the dimensions of the input stack
    ImgSize(mdIn, xdimOut, ydimOut, zdimOut, ndimOut);

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
        if (!fn_out.empty())
            std::cout << "Output File: " << fn_out << std::endl;
        if (!oroot.empty())
            std::cout << "Output Root: " << oroot << std::endl;
    }
}

void XmippMetadataProgram::preProcess()
{}

void XmippMetadataProgram::postProcess()
{}

void XmippMetadataProgram::startProcessing()
{
    if (delete_output_stack)
    {
        FileName fn_stack_plain = fn_out.removeFileFormat();
        fn_stack_plain.deleteFile();
    }

    if (create_empty_stackfile)
        createEmptyFile(fn_out, xdimOut, ydimOut, zdimOut, mdInSize, true, WRITE_OVERWRITE);

    //Show some info
    show();
    // Initialize progress bar
    time_bar_size = mdInSize;
    if (allow_time_bar && verbose && !single_image)
        init_progress_bar(time_bar_size);
    time_bar_step = CEIL((double)time_bar_size / 60.0);
    time_bar_done = 0;
}

void XmippMetadataProgram::finishProcessing()
{
    if (allow_time_bar && verbose && !single_image)
        progress_bar(time_bar_size);

    if (!single_image && !mdOut.isEmpty()&& !fn_out.empty())
    {
        if (produces_an_output || !oroot.empty()) // Out as independent images
            mdOut.write(fn_out);
        else if (save_metadata_stack) // Output is stack and also save its associated metadata
        {
            FileName outFileName;

            outFileName = fn_out.withoutExtension().addExtension("xmd");
            mdOut.write(outFileName);
        }
    }
}

void XmippMetadataProgram::showProgress()
{
    if (time_bar_done % time_bar_step == 0 && allow_time_bar && verbose && !single_image)
        progress_bar(time_bar_done);
}

bool XmippMetadataProgram::getImageToProcess(size_t &objId, size_t &objIndex)
{
    if (time_bar_done == 0)
        iter = new MDIterator(mdIn);
    else
        iter->moveNext();

    ++time_bar_done;
    objIndex = iter->objIndex;
    return ((objId = iter->objId) != BAD_OBJID);
}

void XmippMetadataProgram::run()
{
    FileName fnImg, fnImgOut, baseName, pathBaseName, fullBaseName, oextBaseName;
    size_t objId , outId;
    MDRow rowIn, rowOut;

    //Perform particular preprocessing
    preProcess();

    startProcessing();

    size_t objIndex = 0;

    if (!oroot.empty())
    {
        if (oext.empty())
            oext           = oroot.getFileFormat();
        oextBaseName   = oext;
        fullBaseName   = oroot.removeFileFormat();
        baseName       = fullBaseName.getBaseName();
        pathBaseName   = fullBaseName.getDir();
    }

    //FOR_ALL_OBJECTS_IN_METADATA(mdIn)
    while (getImageToProcess(objId, objIndex))
    {
        ++objIndex; //increment for composing starting at 1

        mdIn.getRow(rowIn, objId);
        rowIn.getValue(MDL_IMAGE, fnImg);

        if (fnImg.empty())
            break;

        fnImgOut = fnImg;

        if (each_image_produces_an_output)
        {
            if (!oroot.empty()) // Compose out name to save as independent images
            {
                if (oext.empty()) // If oext is still empty, then use ext of indep input images
                    oextBaseName = fnImg.getFileFormat();

                if (!baseName.empty() )
                    fnImgOut.compose(fullBaseName, objIndex, oextBaseName);
                else if (fnImg.isInStack())
                    fnImgOut.compose(pathBaseName + (fnImg.withoutExtension()).getDecomposedFileName(), objIndex, oextBaseName);
                else
                    fnImgOut = pathBaseName + fnImg.withoutExtension()+ "." + oextBaseName;
            }
            else if (!fn_out.empty() )
            {
                if (single_image)
                    fnImgOut = fn_out;
                else
                    fnImgOut.compose(objIndex, fn_out); // Compose out name to save as stacks
            }
            else
                fnImgOut = fnImg;

            if (keep_input_columns)
                rowOut = rowIn;
            else
                rowOut.clear();

            rowOut.setValue(MDL_IMAGE, fnImgOut);
            rowOut.setValue(MDL_ENABLED, 1);

        }

        processImage(fnImg, fnImgOut, rowIn, rowOut);

        if (each_image_produces_an_output)
            mdOut.addRow(rowOut);

        showProgress();
    }

    //free iterator memory
    delete iter;

    /* Generate name to save mdOut when output are independent images. It uses as prefix
     * the dirBaseName in order not overwriting files when repeating same command on
     * different directories. If baseName is set it is used, otherwise, input name is used.
     * Then, the suffix _oext is added.*/

    if (fn_out.empty() && !oroot.empty())
        if (!baseName.empty() )
            fn_out = findAndReplace(pathBaseName,"/","_") + baseName + "_" + oextBaseName + ".xmd";
        else
            fn_out = findAndReplace(pathBaseName,"/","_") + fn_in.getBaseName() + "_" + oextBaseName + ".xmd";

    finishProcessing();

    postProcess();
}

XmippProgramGeneric::XmippProgramGeneric()
{
    initComments();
    progDef = new ProgramDef();
    definitionComplete = false;
}

void XmippProgramGeneric::endDefinition()
{
    definitionComplete = true;
    this->defineCommons();
    progDef->parse();
}

void XmippProgramGeneric::read(int argc, char ** argv, bool reportErrors)
{
    if (!definitionComplete)
        endDefinition();
    XmippProgram::read(argc, argv, reportErrors);
}
//All the following are necesary to override the base class implementation
void XmippProgramGeneric::readParams()
{}
void XmippProgramGeneric::defineParams()
{}
void XmippProgramGeneric::show()
{}
void XmippProgramGeneric::run()
{}
