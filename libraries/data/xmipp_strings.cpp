/**************************************************************************
 *
 * Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <regex.h>
#include <stdarg.h>
#include "xmipp_strings.h"
#include "xmipp_error.h"
#include "xmipp_macros.h"
#include "gcc_version.h"

String removeChar( const String& str, char character )
{
    String temp;

    for( unsigned int i = 0 ; i < str.length( ) ; i++ )
    {
        if ( str[ i ] != character )
            temp += str[ i ];
    }

    return temp;
}

String unescape( const String& str )
{
    String temp;

    for( unsigned int i = 0 ; i < str.length( ) ; i++ )
    {
        char current_char = str[ i ];

        if( current_char != '\n' && current_char != '\t' &&
            current_char != '\v' && current_char != '\b' &&
            current_char != '\r' && current_char != '\f' &&
            current_char != '\a' )
        {
            temp += str[ i ];
        }
    }

    return temp;
}

String simplify( const String& str )
{
    String temp;

    // First, unescape string
    String straux = unescape( str );

    // Remove spaces from the beginning
    int pos = straux.find_first_not_of( ' ' );
    straux.erase( 0, pos );

    // Trim the rest of spaces
    for( unsigned int i = 0 ; i < straux.length( ) ; )
    {
        temp += straux[ i ];

        if ( straux[ i ] == ' ' )
        {
            while( straux[ i ] == ' ' )
            {
                i++;
            }
        }
        else
        {
            i++;
        }
    }

    // Remove space left at the end of the string
    // if needed
    if( temp[ temp.size( ) - 1 ] == ' ' )
    {
        temp.resize( temp.size() - 1 );
    }

    return temp;
}

/* Trim all spaces from the begining and the end */
void trim(String& str)
{
    String::size_type pos = str.find_last_not_of(' ');

    if (pos != String::npos)
    {
        str.erase(pos + 1);
        pos = str.find_first_not_of(' ');
        if (pos != String::npos)
            str.erase(0, pos);
    }
    else
        str.clear();
}

float textToFloat(const char* str)
{
    if (str == NULL)
        REPORT_ERROR(ERR_MEM_NULLPOINTER, "Cannot be converted into float");

    float retval;
    int ok = sscanf(str, "%f", &retval);

    if (ok)
        return retval;
    REPORT_ERROR(ERR_VALUE_INCORRECT, "Conversion to float error");

    return 0;
}

int textToInteger(const char* str)
{
    if (str == NULL)
        REPORT_ERROR(ERR_MEM_NULLPOINTER, "Cannot be converted into int");

    int retval;
    int ok = sscanf(str, "%d", &retval);

    if (ok)
        return retval;
    REPORT_ERROR(ERR_VALUE_INCORRECT, "Conversion to int error");

    return 0;
}

size_t textToSizeT(const char * str)
{
    if (str == NULL)
        REPORT_ERROR(ERR_MEM_NULLPOINTER, "Cannot be converted into int");

    long unsigned int retval;
    int ok = sscanf(str, "%lu", &retval);

    if (ok)
        return retval;
    REPORT_ERROR(ERR_VALUE_INCORRECT, "Conversion to size_t error");

    return 0;
}

long long textToLongLong(const char* str)
{
    if (str == NULL)
        REPORT_ERROR(ERR_MEM_NULLPOINTER, "Cannot be converted into long long int");

    long long int retval;
    int ok = sscanf(str, "%lld", &retval);

    if (ok)
        return retval;
    REPORT_ERROR(ERR_VALUE_INCORRECT, "Conversion to long long int error");

    return 0;
}

int bestPrecision(float F, int _width)
{
    // If it is 0
    if (F == 0)
        return 1;

    // Otherwise
    int exp = FLOOR(log10(ABS(F)));
    int advised_prec;

    if (exp >= 0)
        if (exp > _width - 3)
            advised_prec = -1;
        else
            advised_prec = _width - 2;
    else
    {
        advised_prec = _width + (exp - 1) - 3;
        if (advised_prec <= 0)
            advised_prec = -1;
    }

    if (advised_prec < 0)
        advised_prec = -1; // Choose exponential format

    return advised_prec;
}

String floatToString(float F, int _width, int _prec)
{
#if GCC_VERSION < 30300
    char aux[15];
    std::ostrstream outs(aux, sizeof(aux));
#else

    std::ostringstream outs;
#endif

    outs.fill(' ');

    if (_width != 0)
        outs.width(_width);

    if (_prec == 0)
        _prec = bestPrecision(F, _width);

    if (_prec == -1 && _width > 7)
    {
        outs.precision(_width - 7);
        outs.setf(std::ios::scientific);
    }
    else
        outs.precision(_prec);

#if GCC_VERSION < 30301

    outs << F << std::ends;
#else

    outs << F;
#endif

#if GCC_VERSION < 30300

    return String(aux);
#else

    String retval = outs.str();
    int i = retval.find('\0');

    if (i != -1)
        retval = retval.substr(0, i);

    return retval;
#endif
}

String integerToString(int I, int _width, char fill_with)
{
    char aux[15];

    // Check width
    int width = _width;
    int Iaux = ABS(I);

    if (SGN(I) < 0)
        width--;

    if (width == 0)
        do
        {
            Iaux /= 10;
            width++;
        }
        while (Iaux != 0);

    // Fill the number with the fill character
    for (int i = 0; i < width; i++)
        aux[i] = fill_with;

    // Start filling the array
    aux[width--] = '\0';
    Iaux = ABS(I);
    do
    {
        int digit = Iaux % 10;
        Iaux /= 10;
        aux[width--] = '0' + digit;
    }
    while (Iaux != 0);

    if (SGN(I) < 0)
        return static_cast< String >("-")  + aux;
    else
        return static_cast< String >(aux);
}

int charToInt(const char* str)
{
    if (str == NULL)
        REPORT_ERROR(ERR_MEM_NULLPOINTER, "Cannot convert NULL into int");

    char readval;
    int ok = sscanf(str, "%c", &readval);

    if (ok)
        return readval - 48;
    REPORT_ERROR(ERR_VALUE_INCORRECT, "Conversion to int error");

    return 0;
}

String stringToString(const String& str, size_t _width)
{
    if (_width == 0)
        return str;

    if (_width < str.length())
        return str.substr(0, _width);

    String aux = str;
    return aux.append(_width - str.length(), ' ');
}

void checkAngle(const String& str)
{
    if (str == "rot")
        return;

    if (str == "tilt")
        return;

    if (str == "psi")
        return;

    REPORT_ERROR(ERR_VALUE_INCORRECT,
                 static_cast< String >(
                     "checkAngle: Not recognized angle type: " + str));
}

String removeSpaces(const String& _str)
{
    String retval;
    int first = _str.find_first_not_of("\n \t");
    int last = _str.find_last_not_of("\n \t");
    bool after_blank = false;

    for (int i = first; i <= last; i++)
    {
        if (_str[i] == ' ' || _str[i] == '\n' || _str[i] == '\t')
        {
            if (!after_blank)
                retval += _str[i];

            after_blank = true;
        }
        else
        {
            retval += _str[i];
            after_blank = false;
        }
    }

    return retval;
}

// Remove quotes ===========================================================
void removeQuotes(char **_str)
{
    String retval = *_str;
    if (retval.length() == 0)
        return;
    char c = retval[0];
    if (c == '\"' || c == '\'')
        retval = retval.substr(1, retval.length() - 1);
    c = retval[retval.length()-1];
    if (c == '\"' || c == '\'')
        retval = retval.substr(0, retval.length() - 1);
    free(*_str);
    *_str = strdup(retval.c_str());
}

// Split a string ==========================================================
int splitString(const String& input,
                const String& delimiter,
                StringVector & results,
                bool includeEmpties)
{
    results.clear();
    size_t delimiterSize = delimiter.size();
    if (input.size()== 0 || delimiterSize == 0)
        return 0;

    size_t newPos, iPos = 0;
    String emptyString;
    while ((newPos = input.find(delimiter, iPos))!=String::npos)
    {
        if (newPos==iPos)
        {
            if (includeEmpties)
                results.push_back(emptyString);
        }
        else
            results.push_back(input.substr(iPos, newPos-iPos));
        iPos = newPos+delimiterSize;
    }
    if (iPos>=input.size())
    {
        if (includeEmpties)
            results.push_back(emptyString);
    }
    else
        results.push_back(input.substr(iPos, String::npos));
    return results.size();
}

// To lower ================================================================
void toLower(char *_str)
{
    int i = 0;
    while (_str[i] != '\0')
    {
        if (_str[i] >= 'A' && _str[i] <= 'Z')
            _str[i] += 'a' -'A';
        i++;
    }
}

void toLower(String &_str)
{
    int i = 0;
    while (_str[i] != '\0')
    {
        if (_str[i] >= 'A' && _str[i] <= 'Z')
            _str[i] += 'a' -'A';
        i++;
    }
}

// Next token ==============================================================
String nextToken(const String &str, size_t &i)
{
    String retval;
    if (i >= str.length())
        return retval;
    int j = str.find_first_not_of(" \t\n", i);
    if (j == -1)
        return retval;
    int k = str.find_first_of(" \t\n", j + 1);
    if (k == -1)
        k = str.length();
    retval = str.substr(j, k - j + 1);
    i = k + 1;
    return retval;
}

// Get word ================================================================
char *firstWord(char *str)
{
    char *token;

    // Get token
    if (str != NULL)
        token = firstToken(str);
    else
        token = nextToken();

    // Check that there is something
    if (token == NULL)
        REPORT_ERROR(ERR_VALUE_EMPTY, "Empty token");

    return token;
}

// Tokenize a C++ string ===================================================
void tokenize(const String& str, StringVector& tokens,
              const String& delimiters)
{
    tokens.clear();
    // Skip delimiters at beginning.
    String::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    String::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (String::npos != pos || String::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

// Find and replace =======================================================
String findAndReplace(const String& tInput, const String &tFind,
                      const String &tReplace)
{
    size_t uFindLen = tFind.length();

    if( uFindLen == 0 )
        return tInput;

    size_t uPos = 0;
    size_t uReplaceLen = tReplace.length();
    String tOut=tInput;
    for( ;(uPos = tOut.find(tFind, uPos)) != String::npos; )
    {
        tOut=tOut.replace( uPos, uFindLen, tReplace );
        uPos += uReplaceLen;
    }
    return tOut;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/* Tokenizer  for char arrays. It does NOT modify
 * the input array
 *  src is a pointer to the array beginning.
 *  It may be modified to trim the token
 * _end is a pointer to the end of the token
 * sep is an array with the valid separator
 */

char   *memtok(char **src,  char **_end, const char *sep)
{
    char   *start = *src;
    char   *end;

    /*
     * Skip over leading delimiters.
     */
    start += strspn(start, sep);
    if (*start == 0)
    {
        *src = start;
        return (0);
    }

    /*
     * Separate off one token.
     */
    end = start + strcspn(start, sep);
    *src = start;
    *_end=end;
    return (start);
}


#undef _memmem

/* Return the first occurrence of NEEDLE in HAYSTACK. Taken from GNU C Library */
void * _memmem ( const void *haystack, size_t haystack_len, const void *needle, size_t needle_len)
{
    const char *begin;
    const char *const last_possible = (const char *) haystack + haystack_len - needle_len;

    if (needle_len == 0)
        /* The first occurrence of the empty string is deemed to occur at
           the beginning of the string.  */
        return (void *) haystack;

    /* Sanity check, otherwise the loop might search through the whole
       memory.  */
    if (haystack_len < needle_len)
        return NULL;

    for (begin = (const char *) haystack; begin <= last_possible; ++begin)
        if (begin[0] == ((const char *) needle)[0] &&
            !memcmp ((const void *) &begin[1],
                     (const void *) ((const char *) needle + 1),
                     needle_len - 1))
            return (void *) begin;

    return NULL;
}

/* Obtain an string from a format in the way of printf works
 *
 */
String formatString(const char * format, ...)
{
    char formatBuffer[1024];
    va_list args;
    va_start(args, format);
    vsprintf (formatBuffer, format, args);
    String result(formatBuffer);
    va_end (args);

    return result;
}

/* Obtain an string from a format in the way of printf works
 *
 */
void formatStringFast(String &str, const char * format, ...)
{
    char formatBuffer[1024];

    va_list args;
    va_start(args, format);
    vsprintf (formatBuffer, format, args);
    str=formatBuffer;
    va_end (args);
}

/* Matches regular expression */
bool matchRegExp(const String &inputString, const String &pattern)
{
    // Construct regular expression
    regex_t re;
    if (regcomp(&re, pattern.c_str(), REG_EXTENDED|REG_NOSUB) != 0)
        REPORT_ERROR(ERR_ARG_INCORRECT,(String)"Pattern cannot be parsed:"+pattern);

    // Check if the string matches the pattern
    int status = regexec(&re, inputString.c_str(), (size_t) 0, NULL, 0);
    regfree(&re);
    if (status != 0)
        return false;
    return true;
}
#include <sstream>

String WordWrap(const String &inputString, size_t lineLength)
{

    if(inputString.size() <= lineLength)
        return ((String)"# " + inputString + "\n");
    std::istringstream iss(inputString);
    std::ostringstream ss;
    std::string line;

    line.clear();
    do
    {
        std::string word;
        iss >> word;

        if (line.length() + word.length() > lineLength)
        {
            ss << "# "  << line << std::endl;
            line.clear();
        }
        line += word + " ";

    }
    while (iss);

    if (!line.empty())
    {
        ss << "# "  << line << std::endl;
    }
    return ss.str();
}

String escapeForRegularExpressions(const String &str)
{
	String aux;
	aux=findAndReplace(str,"+","\\+");
	return aux;
}
