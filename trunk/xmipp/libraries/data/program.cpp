/***************************************************************************
 * Authors:     AUTHOR_NAME (josem@cnb.csic.es)
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
    addParamsLine("[--more*]         : Show additional options.");

    ///This are a set of internal command for MetaProgram usage
    ///they should be hidden
    addParamsLine("==+++++ Internal section ==");
    addParamsLine("[--xmipp_write_definition* <dbname>] : Print metadata info about the program");

    progDef->parse();
}

void XmippProgram::checkBuiltIns()
{
    ///If -more_options provided, show extended usage
    if (checkParam("--more"))
        usage(1);
    ///If help requested, print usage message
    if (checkParam("--help"))
    {
        std::string cmdHelp = (std::string)"-" + getParam("-h");
        if (cmdHelp == "-")
            usage();
        else
            usage(cmdHelp);
    }
    if (checkParam("--xmipp_write_definition"))
    {
        writeToDB("programs.db");
        exit(0);
    }
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

XmippProgram::XmippProgram()
{
    progDef = NULL;
}

XmippProgram::XmippProgram(int argc, char ** argv)
{
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

void XmippProgram::readParams()
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "function 'readParams'");
}

void XmippProgram::read(int argc, char ** argv, bool reportErrors)
{
    if (progDef == NULL)
        init();

    setProgramName(argv[0]);

    ///If not arguments are provided, show usage message
    if (argc == 1)
        usage();

    try
    {
        this->argc = argc;
        this->argv = argv;
        progDef->read(argc, argv, reportErrors);
        checkBuiltIns();
        verbose = getIntParam("--verbose");
        this->readParams();
    }
    catch (XmippError xe)
    {
        ///If an input error, shows error message and usage
        std::cerr << xe;
        usage();
    }
}

void XmippProgram::tryRun()
{
    try
    {
        this->run();
    }
    catch (XmippError xe)
    {
        std::cout << xe;
        exit(xe.__errno);
    }
}

void XmippProgram::setProgramName(const char * name)
{
    progDef->name = name;
}

void XmippProgram::addUsageLine(const char * line)
{
    progDef->usageComments.addComment(line);
}

void XmippProgram::clearUsage()
{
    progDef->usageComments.clear();
}

void XmippProgram::addParamsLine(const char * line)
{
    progDef->pLexer->addLine((std::string)line);
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
    return progDef->getParam(param, arg);
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

bool XmippProgram::checkParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((std::string)"Doesn't exists param: " + param));
    return paramDef->counter == 1;
}

const char * XmippProgram::name() const
{
    return progDef->name.c_str();
}

void XmippProgram::usage(int verb) const
{
    ConsolePrinter cp;
    cp.printProgram(*progDef, verb);
    exit(1);
}

void XmippProgram::usage(const std::string & param, int verb)
{
    ConsolePrinter cp;
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((std::string)"Doesn't exists param: " + param));
    cp.printParam(*paramDef, verb);
    exit(0);
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
    produces_an_output = false;
    each_image_produces_an_output = false;
    allow_time_bar = true;
    apply_geo = false;
}

void XmippMetadataProgram::defineParams()
{
    addParamsLine(" -i <metadata>   : Input file: metadata, stack, volume or image.");
    addParamsLine(" alias --input;");

    if (each_image_produces_an_output)
    {
        addParamsLine("  [-o <output_file=\"\">]  : Output file: metadata, stack, volume or image.");
        addParamsLine("   alias --output;");
        addParamsLine("  [-oext <extension=\"\">] :  Output file format extension.");
        addParamsLine("  [-oroot <root=\"\">]     : Rootname of output individual images.");
    }
    else if (produces_an_output)
    {
        addParamsLine("  [-o <output_file=\"\">]  : Output file: metadata, stack, volume or image.");
        addParamsLine("   alias --output;");
    }

    if (apply_geo)
    {
        addParamsLine("  [-dont_apply_geo]   : for 2D-images: do not apply transformation stored in the header");
    }
}

void XmippMetadataProgram::readParams()
{
    fn_in = getParam("-i");

    if (produces_an_output)
        fn_out = checkParam("-o") ? getParam("-o") : fn_in;

    if (each_image_produces_an_output)
    {
        fn_out = checkParam("-o") ? getParam("-o") : fn_in;
        oext = getParam("-oext");
        oroot = getParam("-oroot");
    }

    if (fn_in.isMetaData())
    {
        mdIn.read(fn_in, NULL);
        if (mdIn.containsLabel(MDL_ENABLED))
            mdIn.removeObjects(MDValueEQ(MDL_ENABLED, -1));

        if (mdIn.isEmpty())
            REPORT_ERROR(ERR_MD_NOOBJ, "");
    }
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
        if (oext != "")
            std::cout << "Output Extension: " << oext << std::endl;
        if (oroot != "")
            std::cout << "Output Root: " << oroot << std::endl;
    }
}

void XmippMetadataProgram::preProcess()
{}

void XmippMetadataProgram::postProcess()
{}

void XmippMetadataProgram::run()
{
    try
    {
        //Perform particular preprocessing
        preProcess();

        //Show some info
        show();
        if (!fn_in.isMetaData())
        {
            fnImg=fn_in;
            fnImgOut = fn_out;
            processImage();
        }
        else
        {
            // Initialize progress bar
            int i = 0, mdSize = mdIn.size();

            if (allow_time_bar && verbose)
                init_progress_bar(mdSize);
            int istep = CEIL((double)mdSize / 60.0);

            FOR_ALL_OBJECTS_IN_METADATA(mdIn)
            {
                mdIn.getValue(MDL_IMAGE, fnImg);

                if (fnImg == "")
                    break;

                if (each_image_produces_an_output)
                {
                    fnImgOut = fnImg;
                    if (oroot != "")
                        fnImgOut = oroot + fnImgOut.withoutExtension();
                    if (oext != "")
                        fnImgOut = fnImgOut.withoutExtension() + "." + oext;
                    if (fnImgOut!=fnImg)
                    {
                        mdOut.addObject();
                        mdOut.setValue(MDL_IMAGE,fnImgOut);
                    }
                }

                processImage();

                if (i++ % istep == 0 && allow_time_bar && verbose)
                    progress_bar(i);
            }

            if (allow_time_bar && verbose)
                progress_bar(mdSize);

            if (!mdOut.isEmpty())
                mdOut.write(fn_out);
        }
        postProcess();
    }
    catch (XmippError xe)
    {
        std::cout << xe;
        exit(xe.__errno);
    }
}
