/***************************************************************************
*
* Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef ARGS_H
#define ARGS_H

#include "xmipp_strings.h"
#include "xmipp_error.h"

template <typename T>
class Matrix1D;

/** @defgroup Arguments Arguments parsing
 *  @ingroup DataLibrary
 *
 * This set of functions are designed for make easier checking and reading
 * parameters from a string. The single value readings don't modify the value of
 * the input string, while the multiple value readings do. In general, the
 * reading is oriented to be done with a chain of strtok (a general C function),
 * for this reason the string itself is modified by the succesive calls to
 * strtok.
 *
 * The second group (list reading) uses tokens internally, while in the single
 * read functions the token must be given by hand. Anyway, it is not so
 * difficult to read a set of fields from a line with all the normal checks we
 * would like (existence of a parameter, checking for numerical correctness...)
 *
 * The following is an example of its use. This code tries to read a number, a
 * character and a list of numbers from a line. In the example you can also see
 * how to write code for the different error handling methods.
 *
 * @code
 * std::string line;
 * int key, param_no;
 * std::vector< float > data;
 *
 * try
 * {
 * key = textToFloat(firstToken(line), 1602, "Error reading key");
 * param_no = textToInt(nextToken(), 1602, "Error reading number parameters");
 * readFloatList(NULL, param_no, data, 1602, "Error reading doc file line");
 * }
 * catch (XmippError XE)
 * {
 *     std::cout << XE;
 *     DL.line_type = 0;
 *     REPORT_ERROR(1602, "Line is discarded");
 * }
 * @endcode
 */
//@{
/** @name Read lists
 *
 * These functions try to read N values of the desired type into the given
 * structure (either a STL vector of any numerical type by adding the read
 * values at the end or a Matrix1D of any numerical type and then the values are
 * written at PHYSICAL positions 0 ... N-1, the Matrix1D must be big enough to
 * hold all the data since it is not resized internally.
 *
 * If it is not possible to read all parameters an exception is thrown. You can
 * provide the exception in the function call but there are default values.
 * These functions are based in those for reading single values, so it might be
 * possible that these other functions throw other exceptions.
 *
 * The different elements of the list are selected using the tokenizing
 * functions (different elements in the string are separated by spaces), so
 * after the application of this function to the input string, this is modified
 * and NULL characters are introduced as delimiters of the elements.
 *
 * The following code is an example of doing so:
 *
 * @code
 * getline(in_stream, line);
 * readFloatList(line, 10, v1); // Read 10 values from line in v1
 * readFloatList(NULL, 10, v2); // Read NEXT!! 10 values in v2
 * @endcode
 */
//@{
/** List to STL vector.
 */
template <typename T>
void readFloatList(const char* str,
                   int N, std::vector< T >& v)
{
    T valueF;
    char* token;

    token = firstToken(str);
    for (int i = 0; i < N; i++)
    {
        if (token == NULL)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "Cannot convert string into a list of numbers");

        valueF = (T) textToFloat(token);
        v.push_back(valueF);

        if (i != N - 1)
            token = nextToken();
    }
}

/** List to STL vector.
 */
template <typename T>
void readFloatList(const std::string& str,
				   size_t& i,
                   int N,
                   std::vector< T >& v)
{
    T valueF;
    std::string token;

    token = nextToken(str, i);
    for (int j = 0; j < N; j++)
    {
        if (token == "")
            REPORT_ERROR(ERR_VALUE_INCORRECT, "Cannot convert string into list of floats");
        valueF = (T) textToFloat(token.c_str());
        v.push_back(valueF);

        if (j != N - 1)
            token = nextToken(str, i);
    }
}

/** Read list into a Matrix1D.
 */
template <typename T>
void readFloatList(const char* str,
                   int N,
                   Matrix1D< T >& v,
                   int _errno = 2105,
                   std::string errmsg = "Error reading floating list",
                   int exit = 0)
{
    T valueF;
    char* token;

    token = firstToken(str);
    for (int i = 0; i < N; i++)
    {
        if (token == NULL)
        {
            // CO: Should not report other error than the required one
            // std::cout << "Read float list: Number of true parameters doesn't coincide\n";
            REPORT_ERROR(_errno, errmsg);
        }

        try
        {
            valueF = (T) textToFloat(token);
        }
        catch (XmippError)
        {
            REPORT_ERROR(_errno, errmsg);
        }

        v(i) = valueF;
        if (i != N - 1)
            token = nextToken();
    }
}
//@}

/** @name Functions for parsing the command line
 *
 * These functions help you to manage the command line parameters
 */
//@{
/** Get parameters from the command line.
 *
 * This function assumes that the command line is structured in such a way that
 * for each parameter a block like "-param <param_value>" is defined. The label
 * "param" can be substituted by any other one you like. If the parameter is
 * optional then this function allows you to define a default value. If no
 * default value is supplied and the parameter is not specified in the command
 * line, then an exception is thrown. You may change the default exception.
 *
 * You may also indicate that in case of error no exception is raised and force
 * the program to abort (use the exit variable).
 *
 * @code
 * m_param = textToFloat(getParameter(argc, argv, "-m"));
 *
 * // Get compulsory parameter "-m"
 * m_param = textToFloat(getParameter(argc, argv, "-m","2.65"));
 *
 * // Optional parameter, if no parameter is given it takes 2.65 by default
 * m_param = textToFloat(getParameter(argc, argv, "-m", NULL, 6001, "-m parameter not \
 *     found. I'm going out", TRUE);
 *
 * // Compulsory parameter, if not found give an special error message and exit
 * // the program
 *
 * @endcode
 */
const char* getParameter(int argc,
                   const char** argv,
                   const char* param,
                   const char* option = NULL);

/** Get two float parameters after a flag from the command line.
 *
 * An exception is thrown if there are not enough parameters after the flag, if
 * the message is empty then "Not enough parameters after <param>" is shown. The
 * default values must be provided. TRUE is returned if the two values have been
 * found
 */
bool getTwoDoubleParams(int argc,
                        const char** argv,
                        const char* param,
                        double& v1,
                        double& v2,
                        double v1_def,
                        double v2_def);

/** Get 3 float parameters after a flag from the command line.
 *
 * An exception is thrown if there are not enough parameters after the flag, if
 * the message is empty then "Not enough parameters after <param>" is shown. The
 * default values must be provided. TRUE is returned if the two values have been
 * found
 */
bool getThreeDoubleParams(int argc,
                          const char** argv,
                          const char* param,
                          double& v1,
                          double& v2,
                          double& v3,
                          double v1_def,
                          double v2_def,
                          double v3_def);

/** Get boolean parameters from the command line.
 *
 * This function assumes that the command line is structured in such a way that
 * for each parameter a block like "-param" is defined. The label "param" can be
 * substituted by any other one you like. It might be used to look for a boolean
 * parameter, for instance:
 *
 *     -verbose means that verbose functionality is set (TRUE)
 *
 * @code
 * verbose = checkParameter(argc, argv, "-verbose"));
 *
 * // checks if "-verbose" was supplied in the command line. If -verbose was
 * // supplied the function returns TRUE (1), otherwise returns FALSE (0)
 * @endcode
 */
bool checkParameter(int argc, const char** argv, const char* param);

/** Returns the position where the given parameter is in the command line.
 *
 * This function assumes that the command line is structured in such a way that
 * for each parameter a block like "-param" is defined. The label "param" can be
 * substituted by any other one you like. It returns -1 if the parameter is not
 * found. It is used to look for parameters having a list of values behind, for
 * instance:
 *
 *     -ang rot tilt psi
 *
 * @code
 * i = paremeterPosition(argc, argv, "-ang"));
 *
 * // This condition checks if 3 arguments where introduced after -ang parameter
 * // (assuming that the -ang argument is the last one in the string)
 * if (i+3 >= argc)
 *     EXIT_ERROR(1, "Not enough parameters behind -ang");
 *
 * ang1 = argv[i+1];
 * ang2 = argv[i+2];
 * ang3 = argv[i+3];
 * @endcode
 */
int paremeterPosition(int argc, const char** argv, const char* param);

/** Return the number of components of a vector argument.
 *
 * A vector argument is defined as [x,y,z,...]. It returns 0 if the string does
 * not contain a vector
 */
int numComponents(const std::string& str);

/** Get float vector.
 *
 * A vector is defined as a "[x,y,z, ...]" of any dimension (by default 2D
 * vectors are read). The vector must not contain blank spaces.
 *
 * @code
 * a = getVectorParameter(argc, argv, "-a", 3); // a 3D vector called "-a".
 * @endcode
 *
 * The vector is internally resized properly. If the dimension is -1 then all
 * the vector components are read disregarding their dimensionality, ie,-1 is
 * used to read vectors of an unknown dimension. If the parameter is not found
 * when no dimensionality is given an empty vector is returned but no exception
 * is thrown. If there is no dimensionality and a single parameter is behind the
 * flag then no brackets are needed
 */
Matrix1D< double > getVectorParameter(int argc,
                                      const char** argv,
                                      const char* param,
                                      int dim = 2);

/** Get float vector.
 *
 * Same as the previous function but from a file.
 */
Matrix1D<double> getVectorParameter(FILE *fh, const char *param, int dim=2);

/** Generate argc and argv for a string.
*
* Given a string this function makes a copy of the string and divides it into
* tokens such that they can be used as argc and argv, as if it were a command
* line.
*
* The string copy remains in "copy" and it can be freed by disposing this
* variable.
*
* argvp[0] (normally the program name) is set to any value, in this case to
* "autom", standing for "automatically generated".
*
* argcp==0 if now valid command line is provided, ie, the line is empty or only
* with blanks.
*
* Next time the function is called it checks that argv and copy are empty
* (pointing to NULL), if they aren't then firstly the associated memory is
* freed.
*
* @code
* int argcp;
* char** argvp;
* char* copy;
*
* copy = NULL;
* argvp = NULL;
*
* string command_line = "-i input_file -o output_file";
*
* generateCommandLine(command_line, argcp, &argvp, &copy);
*
* if (argcp != 0)
* read_parameters(argcp, argvp);
* @endcode
*/
void generateCommandLine(const std::string& command_line,
                         int& argcp,
                         char**& argvp,
                         char*& copy);

/** Generate articial command line from a file.
*
* The copy variable must be destroyed outside by "delete copy". This function
* takes "input_file=<input_file>" and turns it into "-input_file <input_file>"
* The appropiate argc, argv are also returned.
*
* Returns TRUE if the parameter is found in the file, and FALSE if it is not
*/
bool generateCommandLine(FILE* fh,
                         const char* param,
                         int& argcp,
                         char**& argvp,
                         char*& copy);

/** Get parameter from file.
 *
 * Parameters are supposed to be identified with an =, so any line which doesn't
 * contain an = character cannot contain a parameter. The left part of the "="
 * is the identifier, and the right part the value. The value is returned
 * without any extra space, if it is compound of several words, then the spaces
 * in between are simplified to a single blank space.
 *
 * The file position inside the file is not moved and comments are allowed
 * starting with "#" and ";".
 *
 * Parameter skip controls the number of these parameters to skip before
 * returning the value, ie, if there are several "parameter=" tags in a file and
 * you want the first one then you should skip 0, if you want the second the
 * skip=1, ...
 *
 * The meaning and use of the exit, errors and optional value is the same as in
 * the command line getParameter
 */
std::string getParameter(FILE* fh,
                         const char* param,
                         int skip = 0,
                         const char* option = NULL);

/** Check if a parameter is present in a file.
 *
 * The same as the previous function, but this function only reports if a
 * parameter is present or not in a file. Notice that boolean parameters must be
 * defined as "parameter=". If after the parameter, "no" comes then this
 * function returns FALSE
 */
bool checkParameter(FILE* fh, const char* param);
//@}
//@}
#endif
