/***************************************************************************
 * Authors:     J.M.de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#ifndef ARGSPARSER_H_
#define ARGSPARSER_H_

//#include <cstring>
//#include <iostream>
//#include <sstream>
//#include <vector>
#include <map>
//#include <algorithm>
//#include <stdio.h>
//#include <stdlib.h>
#include "xmipp_error.h"
#include "xmipp_strings.h"


typedef std::vector<String> StringVector;

//TODO (MARIANA) Please give more documentation and in a good structure e.g. @name (see args.h as example)
/* MARIANA: I defined the name of this group. Please define it as you want. */

/** @defgroup Arguments1 Arguments parser
 *  @ingroup DataLibrary
 * @{
 */

/** Type of tokens for lexical analysis */
typedef enum {
    TOK_ID,     //identifier
    TOK_OPT,    // -ID or --ID
    TOK_INT,    //integer number
    TOK_FLOAT,  //float number
    TOK_STR,    //string value
    TOK_MINUS,  // - or --
    TOK_PLUS,   // +
    TOK_EQ,     // =
    TOK_COMM,   // Comment, from : until end of line
    TOK_LAN,    // <
    TOK_RAN,    // >
    TOK_LBRA,   // [
    TOK_RBRA,   // ]
    TOK_SEMI,   // semicolor ;
    TOK_COMMA,  // ,
    TOK_ETC,    // ...
    TOK_END,    // end of input
    ///Reserved words
    TOK_WHERE,  // 'WHERE' keyword
    TOK_ALIAS,  // 'ALIAS' keyword
    TOK_OR,     // 'OR' keyword
    TOK_REQUIRES,// 'REQUIRES' keyword
    TOK_SECTION // section defined by == Section ==
} ArgTokenType;

class ProgramDef;

/** Just a simple struct to hold information about tokens */
class ArgToken
{
public:
    ArgTokenType type; ///< Type of the token
    String lexeme; ///< the string literal value of the token
    int line; ///< line where token was found
    int start, end; ///< start and end position of the lexeme.
    /// this info will be used by parser and printers
    /// 0 - means visible
    /// 1 - less visible
    /// while great is the number is less visible
    int visibility;
    /// Some special mark to tokens
    bool starred;

    static const char * typeString(ArgTokenType type);
};

/** This class will split the input stream into tokens.
 * The tokens will be further used by the Parser
 * to build the syntax tree.
 */
class ArgLexer
{
private:
    StringVector input;
    size_t line;
    size_t pos; ///< reading position of the input
    size_t offset; ///< the offset from pos of each token
    ArgToken * pToken; ///< pointer to the input token.
    ///Dictionary for identify reserved words
    std::map<String, ArgTokenType> reservedWords;

    //Some utils functions
    void readSpaces();
    void readDigits();
    void readId();
    void setupToken(ArgTokenType type);
    void checkVisibility();
    void checkIndependent();
    void nextLine();

public:
    /** Constructor */
    ArgLexer();
    /** Destructor */
    ~ArgLexer();
    /** Add input lines to the lexer */
    void addLine(const String &line);
    /** Function to parse a new token.
     * If the token is TOK_END will return false
     * and true otherwise. The current token will be changed.
     */
    bool nextToken();
    ArgToken * currentToken() const;
    ArgTokenType lookahead() const;

};

/** Just a class for holding comments */
class CommentList
{
public:
    StringVector comments;
    std::vector<int> visibility;
    std::vector<bool> wikiVerbatim;

    void addComment(const String &comment, int visible = 0, bool wikiVerbatim=false);
    void addComment(const char * comment, bool verbatim=false);
    void clear();
    size_t size() const;
};

/** Following classes represent the Abstract Syntax Tree
 * for the language of definition of a program.
 */

/** Class representing the nodes of the tree.
 * All nodes will have the parse method, which need
 * a ArgLexer to ask for tokens.
 * Also a 'consume' method to use the terminal symbols
 */
class ASTNode
{
public:
    ASTNode(ArgLexer *lexer = NULL, ASTNode * parent = NULL);
    virtual ~ASTNode()
    {
    }
    ;

    ASTNode * parent;
    ArgLexer * pLexer;
    ArgToken token;
    String name;
    int visible;

    virtual bool parse() = 0; //abstract function
    virtual void check(std::stringstream & errors) = 0; //abstract function
    virtual bool consume(ArgTokenType type);
    ArgTokenType lookahead() const;
    bool lookahead(ArgTokenType type) const;
    ArgToken * currentToken() const;
    void nextToken();
    bool parseCommentList(CommentList &comments);
    void error(String msg);
    void unexpectedToken(String msg = "");
};

class ParamDef;

class ArgumentDef: public ASTNode
{
public:
    String argDefault;
    bool isList;
    bool isType;
    std::vector<ParamDef*> subParams;
    bool hasDefault;

    ArgumentDef(ArgLexer *lexer, ASTNode * parent);
    ~ArgumentDef();
    virtual bool parse();
    virtual void check(std::stringstream & errors)
    {
    }
    /// This function will take an index and check if there are enougth arguments
    // to pass to this parameter and increase the index
    bool acceptArguments(std::stringstream &errors, size_t &argIndex, std::vector<const char *> &cmdArguments);
};

/** Class representing the definition of an param
 * An param definition is in the form:
 * */
class ParamDef: public ASTNode
{
public:
    bool notOptional; //contradictory param not paramal :)
    bool orBefore;
    bool independent;
    std::vector<ArgumentDef*> arguments;
    std::vector<const char *> cmdArguments;
    std::vector<ParamDef*> *exclusiveGroup;
    int counter; ///< for count the number of times it appears in command line

    CommentList comments;
    StringVector aliases;
    StringVector requires;

    //Empty constructor
    ParamDef(ArgLexer *lexer, ASTNode * parent);
    ~ParamDef();
    bool parseParamList(ArgTokenType startToken, ProgramDef * prog, StringVector &paramList, bool addName);
    bool parseArgumentList();
    virtual bool parse();
    bool checkRequires(std::stringstream & errors, ProgramDef * prog);
    virtual void check(std::stringstream & errors);
    bool containsArgument(const String & argName);
    ArgumentDef * findArgument(const String & argName);
    bool containsAlias(const String & alias);
};

class SectionDef: public ASTNode
{
public:
    CommentList comments;
    std::vector<ParamDef*> params; ///< All params defined for the program.

    SectionDef(ArgLexer * lexer, ASTNode * parent);
    ~SectionDef();
    virtual bool parse();
    virtual void check(std::stringstream & errors)
    {
    }
    ///Add a param to the section
    void addParamDef(ParamDef * param);
};

class ProgramDef: public ASTNode
{
private:
  std::vector<ParamDef*> *exclusiveGroup;
public:
    std::vector<SectionDef*> sections;
    CommentList usageComments; ///< comments of usage
    CommentList examples; ///< examples of use
    std::map<String, ParamDef*> paramsMap; ///< Dictionary with all params and alias names
    StringVector pendingRequires; ///< This is for checking that requires names exists
    String keywords;
    String seeAlso;
    ///This flag is used to check if an independent option was found like: --more, --help
    ///that avoid others options restrictions.
    bool singleOption;

    ProgramDef();
    ~ProgramDef();
    virtual bool parse();
    virtual void check(std::stringstream & errors);
    ParamDef * findParam(const String &param);
    /** Find a param and if not provided in cmd line, fill with its defaults values */
    ParamDef * findAndFillParam(const String &param);
    const char * getParam(const char * paramName, size_t paramNumber = 0);
    const char * getParam(const char * paramName, const char * subParam, size_t paramNumber = 0);
    void addParamName(const String & name, ParamDef *param);
    void addParamRequires(const String &name);
    void addParamExclusiveGroup(ParamDef * param);
    ///clear read arguments
    void clear();
    /// Read and validate commmand line
    void read(int argc, const char ** argv, bool reportErrors = true);
    /// Add a section to the program definition
    /// and return a pointer to it, usefull for manually
    SectionDef * addSection(String sectionName, int visibility = 0);

};

/** @} */
#endif /* ARGSPARSER_H_ */
