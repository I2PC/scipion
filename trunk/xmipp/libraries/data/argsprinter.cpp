/***************************************************************************
 * Authors:     J.M de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#include "argsprinter.h"
#include "xmipp_filename.h"
#include "xmipp_program.h"
#include "xmipp_color.h"

//-------------------   PRINTER IMPLEMENTATIONS   --------------------------------
void Printer::printToken(ArgToken * token)
{

    std::cerr << "token: '" << token->lexeme
    << "' type: " << ArgToken::typeString(token->type)
    << " line: " << token->line + 1
    << " pos: " << token->start + 1 << std::endl;
}

//--------- CONSOLE PRINTER -----------------------
#define COLOR(x, c) (color ? colorString(x, c) : String(x))

ConsolePrinter::ConsolePrinter(std::ostream & out, bool color)
{
    this->pOut = &out;
    this->color = color;
}

void ConsolePrinter::printProgram(const ProgramDef &program, int v)
{
    //print program name and usage
    *pOut << COLOR("PROGRAM", RED) << std::endl << "   " << program.name << std::endl;
    if (program.usageComments.size() > 0)
    {
        *pOut << COLOR("USAGE", RED) << std::endl;
        for (size_t i = 0; i < program.usageComments.size(); ++i)
            if (program.usageComments.visibility[i] <= v)
                *pOut << "   " << program.usageComments.comments[i] << std::endl;
    }
    //print see also
    if (!program.seeAlso.empty())
    {
        *pOut << COLOR("SEE ALSO", RED) << std::endl;
        *pOut << "   " << program.seeAlso << std::endl;
    }

    //print sections and params
    if (program.sections.size() > 0)
    {
        *pOut << COLOR("OPTIONS", RED) << std::endl;
        for (size_t i = 0; i < program.sections.size(); ++i)
            printSection(*program.sections[i], v);
    }
    //print examples
    if (program.examples.size() > 0)
    {
        *pOut << COLOR("EXAMPLES", RED) << std::endl;
        for (size_t i = 0; i < program.examples.size(); ++i)
            if (program.examples.visibility[i] <= v)
            {
                if (program.examples.wikiVerbatim[i])
                    *pOut << "      " << COLOR(program.examples.comments[i].c_str(), BLUE) << std::endl;
                else
                    *pOut << "   " << program.examples.comments[i] << std::endl;
            }

    }
}

void ConsolePrinter::printSection(const SectionDef &section, int v)
{
    if (section.visible <= v)
    {
        *pOut << std::endl;
        if (section.name.length() > 0)
            *pOut << COLOR(section.name.c_str(), RED) << std::endl;
        for (size_t i = 0; i < section.params.size(); ++i)
            printParam(*section.params[i], v);
    }
}

void ConsolePrinter::printRequiresList(StringVector requires)
{
    if (!requires.empty())
    {
        *pOut << " ( requires ";
        for (size_t i = 0; i < requires.size(); ++i)
            *pOut << requires[i] << " ";
        *pOut << ")";
    }
}

void ConsolePrinter::printParam(const ParamDef &param, int v)
{
    if (param.visible <= v)
    {
        if (param.orBefore)
            *pOut << "   OR" << std::endl;

        *pOut << "   ";
        int pColor = BLUE;
        if (!param.notOptional)
        {
            *pOut << "[";
            pColor = GREEN;
        }
        *pOut << COLOR(param.name.c_str(), pColor);
        //print alias
        for (size_t i = 0; i < param.aliases.size(); ++i)
            *pOut << ", " << param.aliases[i];
        //print arguments
        for (size_t i = 0; i < param.arguments.size(); ++i)
        {
            *pOut << " ";
            printArgument(*param.arguments[i], v);
        }
        if (!param.notOptional)
            *pOut << "]";

        printRequiresList(param.requires);
        *pOut << std::endl;
        printCommentList(param.comments, v);

        for (size_t i = 0; i < param.arguments.size(); ++i)
        {
            ArgumentDef &arg = *param.arguments[i];
            if (!arg.subParams.empty())
            {
                *pOut << "      where <" << arg.name << "> can be:" << std::endl;
                for (size_t j = 0; j < arg.subParams.size(); ++j)
                {
                    *pOut << "        " << arg.subParams[j]->name;
                    for (size_t k = 0; k < arg.subParams[j]->arguments.size(); ++k)
                    {
                        *pOut << " ";
                        printArgument(*(arg.subParams[j]->arguments[k]));
                    }
                    printRequiresList(arg.subParams[j]->requires);
                    *pOut << std::endl;
                    printCommentList(arg.subParams[j]->comments);

                }
            }
        }

    }
}

void ConsolePrinter::printArgument(const ArgumentDef & argument, int v)
{
    *pOut << "<" << argument.name;
    if (argument.hasDefault)
        *pOut << "=" << argument.argDefault;
    *pOut << ">";
}

void ConsolePrinter::printCommentList(const CommentList &comments, int v)
{
    for (size_t i = 0; i < comments.size(); ++i)
        if (comments.visibility[i] <= v)
            *pOut << "          " << comments.comments[i] << std::endl;
}

//-------------------   TK PRINTER IMPLEMENTATIONS   --------------------------------

TkPrinter::TkPrinter()
{
    FileName dir = xmippBaseDir();
    dir.append("/applications/scripts/program_gui/program_gui.py");
    output = popen(dir.c_str(), "w");
}

TkPrinter::~TkPrinter()
{
    pclose(output);
}

void TkPrinter::printProgram(const ProgramDef &program, int v)
{
    //    *pOut << "PROGRAM" << std::endl << "   " << program.name << std::endl;
    fprintf(output, "XMIPP %d.%d - %s\n", XMIPP_MAJOR, XMIPP_MINOR, program.name.c_str());
    size_t numberOfComments = 0;
    for (size_t i = 0; i < program.usageComments.size(); ++i)
        if (program.usageComments.visibility[i] <= v)
            ++numberOfComments;
    //Send number of usage lines
    fprintf(output, "%d\n", numberOfComments);
    if (numberOfComments > 0)
    {
        for (size_t i = 0; i < program.usageComments.size(); ++i)
            if (program.usageComments.visibility[i] <= v)
                fprintf(output, "%s\n", program.usageComments.comments[i].c_str());
    }

    for (size_t i = 0; i < program.sections.size(); ++i)
        printSection(*program.sections[i], v);
}

void TkPrinter::printSection(const SectionDef &section, int v)
{
    if (section.visible <= v)
    {
        //Just ignore in the GUI this section
        if (section.name == " Common options ")
            return;

        //if (section.name.length() > 0)
        fprintf(output, "section = self.addSection('%s');\n", section.name.c_str());
        bool first_group = true;
        for (size_t i = 0; i < section.params.size(); ++i)
        {
            if (section.params[i]->visible <= v)
            {
                if (!section.params[i]->orBefore)
                {
                    const char * single = (i < section.params.size()-1 && section.params[i+1]->orBefore) ? "False" : "True";

                    if (!first_group)
                        fprintf(output, "section.addGroup(group);\n");
                    else
                        first_group = false;
                    fprintf(output, "group = ParamsGroup(section, %s);\n", single);

                }
                printParam(*section.params[i], v);
            }
        }
        //close last open group
        if (!first_group)
            fprintf(output, "section.addGroup(group);\n");
    }
}

void TkPrinter::printParam(const ParamDef &param, int v)
{
    if (param.visible <= v)
    {
        //Independent params are some kind of special ones
        if (param.independent)
            return;

        int n_args = param.arguments.size();

        fprintf(output, "param = ParamWidget(group, \"%s\");\n", param.name.c_str());
        if (param.notOptional)
            fprintf(output, "param.notOptional = True; \n");
        for (size_t i = 0; i < param.arguments.size(); ++i)
        {
            printArgument(*param.arguments[i], v);
        }
        //Add comments to the help
        for (size_t i = 0; i < param.comments.size(); ++i)
            //if (param.comments.visibility[i] <= v)
            fprintf(output, "param.addCommentLine('''%s''');\n", param.comments.comments[i].c_str());
        //End with options of the param
        fprintf(output, "param.endWithOptions();\n");

    }
}

void TkPrinter::printArgument(const ArgumentDef & argument, int v)
{
    static String paramStr = "param";
    fprintf(output, "%s.addOption(\"%s\", \"%s\", %d);\n",
            paramStr.c_str(), argument.name.c_str(), argument.argDefault.c_str(), argument.subParams.size());
    if (argument.subParams.size() > 0)
    {

        for (size_t j = 0; j < argument.subParams.size(); ++j)
        {
            fprintf(output, "subparam = param.addSubParam(\"%s\");\n",
                    argument.subParams[j]->name.c_str());
            for (size_t i = 0; i < argument.subParams[j]->comments.size(); ++i)
                fprintf(output, "subparam.addCommentLine('''%s''');\n", argument.subParams[j]->comments.comments[i].c_str());
            paramStr = "subparam";
            for (size_t k = 0; k < argument.subParams[j]->arguments.size(); ++k)
            {
                printArgument(*(argument.subParams[j]->arguments[k]));
            }
            paramStr = "param";
        }
    }
}

//--------- WIKI  PRINTER -----------------------
WikiPrinter::WikiPrinter(std::ostream & out)
{
    this->pOut = &out;
}

void WikiPrinter::printProgram(const ProgramDef &program, int v)
{
    //print program name and usage
    *pOut << "---+ !!" << program.name << " (v" << XMIPP_MAJOR <<"." << XMIPP_MINOR << ")" << std::endl;
    *pOut << "%TOC%" << std::endl;
    //print usage
    if (program.usageComments.size() > 0)
    {
        *pOut << "---++ Usage" << std::endl;
        for (size_t i = 0; i < program.usageComments.size(); ++i)
            if (program.usageComments.wikiVerbatim[i])
                *pOut << "   <pre>" << program.usageComments.comments[i] << "</pre>\n";
            else
                *pOut << "   " << program.usageComments.comments[i] << std::endl;
    }
    if (!program.seeAlso.empty())
    {
        *pOut << std::endl << "*See also* %BR%" << std::endl;
        StringVector links;
        splitString(program.seeAlso, ",", links);
        for (int i = 0; i < links.size(); ++i)
            *pOut << "[[" << links[i] << "_v" << XMIPP_MAJOR << "][" << links[i] <<"]]  ";
        *pOut << "%BR%" << std::endl;
    }
    //print sections and params
    if (program.sections.size() > 0)
    {
        *pOut << std::endl << "*Parameters*" << std::endl;
        for (size_t i = 0; i < program.sections.size(); ++i)
            printSection(*program.sections[i], v);
    }
    //print examples
    if (program.examples.size() > 0)
    {
        *pOut << "---++ Examples and notes" << std::endl;
        bool verbatim = false;
        for (size_t i = 0; i < program.examples.size(); ++i)
        {
            if (program.examples.wikiVerbatim[i])
            {
                if (!verbatim)
                {
                    *pOut << "<pre>" << std::endl;
                    verbatim = true;
                }
            }
            else
            {
                if (verbatim)
                {
                    *pOut << "</pre>" << std::endl;
                    verbatim = false;
                }
            }
            *pOut << "   " << program.examples.comments[i] << std::endl;
        }
        if (verbatim)
            *pOut << "</pre>" << std::endl;
    }
    //print user comments
    *pOut << "---++ User's comments" << std::endl;
    *pOut << "%COMMENT{type=\"tableappend\"}%" << std::endl;
}

void WikiPrinter::printSection(const SectionDef &section, int v)
{
    if (section.name != " Common options "
        && section.visible <= v)
    {
        *pOut << std::endl;
        String name = section.name;
        trim(name);
        if (name.length() > 0)
            *pOut << "_" << name << "_" << std::endl;
        for (size_t i = 0; i < section.params.size(); ++i)
            printParam(*section.params[i], v);
    }
}

void WikiPrinter::printRequiresList(StringVector requires)
{
    if (!requires.empty())
    {
        *pOut << " ( requires ";
        for (size_t i = 0; i < requires.size(); ++i)
            *pOut << requires[i] << " ";
        *pOut << ")";
    }
}

void WikiPrinter::printParam(const ParamDef &param, int v)
{
    //* =%BLUE%-i [selfile] %ENDCOLOR%= This file contains all the images that are to build the 3D reconstruction
    //* =%GREEN% -o [output file root name] %ENDCOLOR%= If you don't supply this parameter, the same as the input selection one is taken without extension. If you give, for instance, =-o art0001= the following files are created:
    if (param.visible <= v)
    {
        *pOut << "   $";

        if (param.orBefore)
            *pOut << " or";

        String color = param.notOptional ? "BLUE" : "GREEN";

        *pOut << " =%" << color << "%" << param.name;
        //print alias
        for (size_t i = 0; i < param.aliases.size(); ++i)
            *pOut << ", " << param.aliases[i];
        //print arguments
        for (size_t i = 0; i < param.arguments.size(); ++i)
        {
            *pOut << " ";
            printArgument(*param.arguments[i], v);
        }
        *pOut << " %ENDCOLOR%=";
        printRequiresList(param.requires);
        *pOut <<": " ;
        printCommentList(param.comments, v);

        if (param.comments.size() == 0)
            *pOut << "%BR%" << std::endl;

        for (size_t i = 0; i < param.arguments.size(); ++i)
        {
            ArgumentDef &arg = *param.arguments[i];
            if (!arg.subParams.empty())
            {
                *pOut << "      where &lt;" << arg.name << "&gt; can be:" << std::endl;
                for (size_t j = 0; j < arg.subParams.size(); ++j)
                {
                    *pOut << "      * %MAROON% " << arg.subParams[j]->name;
                    for (size_t k = 0; k < arg.subParams[j]->arguments.size(); ++k)
                    {
                        *pOut << " ";
                        printArgument(*(arg.subParams[j]->arguments[k]));
                    }
                    *pOut << " %ENDCOLOR%" << std::endl;
                    printRequiresList(arg.subParams[j]->requires);
                    //                    *pOut << std::endl;
                    printCommentList(arg.subParams[j]->comments);

                }
            }
        }

    }
}

void WikiPrinter::printArgument(const ArgumentDef & argument, int v)
{
    *pOut << "&lt;" << argument.name;
    if (argument.hasDefault)
        *pOut << "=" << argument.argDefault;
    *pOut << "&gt;";
}

void WikiPrinter::printCommentList(const CommentList &comments, int v)
{
    *pOut << "          ";
    for (size_t i = 0; i < comments.size(); ++i)
        if (comments.visibility[i] <= v)
            *pOut << comments.comments[i] << " ";
    *pOut << "%BR%" << std::endl;
}

//-------------------   PROTOCOL PRINTER IMPLEMENTATIONS   --------------------------------

ProtPrinter::ProtPrinter()
{
    FileName dir = xmippBaseDir();
    dir.append("/test_header.py");
    output = fopen(dir.c_str(), "w+");
}

ProtPrinter::~ProtPrinter()
{
    fclose(output);
}

void printBegin(FILE * output, const ProgramDef &program)
{
    fprintf(output, "# -------------------------------------------------------------------------------\n");
    fprintf(output, "# Protocol header automatically generated for program: %s\n", program.name.c_str());
    size_t numberOfComments = program.usageComments.size();
    for (size_t i = 0; i < numberOfComments; ++i)
    {
        fprintf(output, "#   %s\n", program.usageComments.comments[i].c_str());
        printf("#   %s\n", program.usageComments.comments[i].c_str());
    }
    fprintf(output,
            "# -------------------------------------------------------------------------------\n"
            "# {begin_of_header}\n"
            "#------------------------------------------------------------------------------------------\n"
            "# {section}{has_question} Comment\n"
            "#------------------------------------------------------------------------------------------\n"
            "# Display comment\n"
            "DisplayComment = False\n\n"
            "# {text} Write a comment:\n"
            "\"\"\"\nDescribe your run here...\n\"\"\"\n"
            "#-----------------------------------------------------------------------------\n"
            "# {section} Run parameters\n"
            "#-----------------------------------------------------------------------------\n"
            "# Run name:\n"
            "\"\"\" This will identify your protocol run. It need to be unique for each protocol. You could have run1, run2 for protocol X, but not two\n"
            "run1 for it. This name together with the protocol output folder will determine the working dir for this run.\n\"\"\"\n"
            "RunName = \"run_001\"\n\n");
}

void printEnd(FILE * output)
{
    fprintf(output,
            "# {hidden} Show expert options"
            "\"\"\"If True, expert options will be displayed\"\"\"\n"
            "ShowExpertOptions = False\n"
            "#------------------------------------------------------------------------------------------\n"
            "# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE\n"
            "#------------------------------------------------------------------------------------------------\n\n"
            "from protocol_dummy import *\n\n"
            "if __name__ == '__main__':\n"
            "    protocolMain(ProtDummy)\n");
}

void ProtPrinter::printProgram(const ProgramDef &program, int v)
{
    printBegin(output, program);
    for (size_t i = 0; i < program.sections.size(); ++i)
        printSection(*program.sections[i], v);
    printEnd(output);
}

void ProtPrinter::printSection(const SectionDef &section, int v)
{
    //Just ignore in the GUI this section
    if (section.name == " Common options ")
        return;

    bool expert = section.visible > v;
    fprintf(output,
            "#------------------------------------------------------------------------------------------\n"
            "# {section}{has_question} %s %s\n", (expert ? "{expert}" : ""), section.name.c_str());
    fprintf(output, "# Show\nshow_%s = True", section.name.c_str());
    for (size_t i = 0; i < section.params.size(); ++i)
        printParam(*section.params[i], v);
}

void ProtPrinter::printParam(const ParamDef &param, int v)
{
    //Independent params are some kind of special ones
    if (param.independent)
        return;
    bool expert = param.visible > v;
    label = param.name;
    parent_name = removeChar(param.name, '-');

    size_t n_args = param.arguments.size();
    for (size_t i = 0; i < n_args; ++i)
    {
        printArgument(*param.arguments[i], v);
        condition = "";
    }
}

void ProtPrinter::addCondition(const String &newcondition)
{
  if (!condition.empty())
    condition += ",";
  condition += newcondition;
}

void ProtPrinter::printArgument(const ArgumentDef & argument, int v)
{
    String tags = "";
    if (argument.subParams.size() > 0)
    {
        tags = "{list_combo}(" + argument.subParams[0]->name;
        for (size_t j = 1; j < argument.subParams.size(); ++j)
            tags += "," + argument.subParams[j]->name;
        tags += ")";
    }
    if (!condition.empty())
    {
      tags += formatString("{condition}(%s)", condition.c_str());
    }

    label += "   " + argument.name;
    fprintf(output, "# %s %s\n", tags.c_str(), label.c_str());
    String tmp = parent_name;
    parent_name = formatString("%s_%s", parent_name.c_str(), argument.name.c_str());
    std::cerr << "DEBUG_JM: label: " << label << std::endl;

    if (label.find('-') == 0) //In this case is the first argument of the param and should print help
    {
        ParamDef * parent = (ParamDef*)argument.parent;
        fprintf(output, "\"\"\"\n", parent->name.c_str());
        for (size_t i = 0; i < parent->comments.size(); ++i)
            fprintf(output, "%s\n", parent->comments.comments[i].c_str());
        fprintf(output, "\"\"\"\n");
    }
    fprintf(output, "%s = \"%s\"\n\n",
            parent_name.c_str(), argument.argDefault.c_str());

    if (argument.subParams.size() > 0)
    {
        for (size_t j = 1; j < argument.subParams.size(); ++j)
        {
            String condBackup = condition;
            addCondition(formatString("%s=%s", parent_name.c_str(), argument.subParams[j]->name.c_str()));
            for (size_t k = 0; k < argument.subParams[j]->arguments.size(); ++k)
            {
                label = argument.subParams[j]->name;
                printArgument(*(argument.subParams[j]->arguments[k]));
            }
        }
    }
    parent_name = tmp;
    label = "";
}

