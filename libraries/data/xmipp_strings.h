/***************************************************************************
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

#ifndef XMIPP_STRINGS_H
#define XMIPP_STRINGS_H

#include <vector>
#include <string>
#include <string.h>

//TODO: For now just a typedef, I think that would be worth to write an String class
typedef std::string String;
typedef std::vector<String> StringVector;

/// @defgroup StringUtilities String utilities
/// @ingroup DataLibrary
//@{

/** Macro to test if to string are the same */
#define STR_EQUAL(str1, str2) (strcmp((str1), (str2)) == 0)

//@name String processing
//@{

/** Removes all occurrences of 'character' from the string no matter
where they are */
String removeChar( const String& str, char character );

/** Removes escaped symbols ESC+n, t, v, b, r, f, and a */
String unescape( const String& str );

/** Best precision for a float number.
 *
 * This function returns the best precision to be used in a "printf" format if
 * this number is to fit in a given width. It returns -1 if the exponential
 * format is advised.
 *
 * @code
 * template<typename T>
 * std::ostream& operator<<(std::ostream& out, const T& val)
 * {
 *     int i,j;
 *
 *     if (val.xdim == 0)
 *         out << "NULL matrix" << std::endl;
 *     else
 *         out << std::endl;
 *
 *     T aux = ABSnD(val);
 *     int prec = bestPrecision(aux.max(), 10);
 *
 *     for (i=STARTINGY(val); i<=FINISHINGY(val); i++)
 *     {
 *         for (j=STARTINGX(val); j<=FINISHINGX(val); j++)
 *         {
 *             out << floatToString((float) val(i,j), 10, prec) << ' ';
 *         }
 *         out << std::endl;
 *     }
 *
 *     return out;
 * }
 *
 * @endcode
 */
int bestPrecision(float F, int _width);

/** String (char*) to float conversion.
 *
 * @code
 * float key = textToFloat(firstToken(line), 1602, "Error reading key");
 * @endcode
 */
float textToFloat(const char* str);

/** String (String) to integer conversion. */
inline float textToFloat(const String& str)
{
	return textToFloat(str.c_str());
}

/** String (char*) to integer conversion.
 *
 * @code
 * int param_no = textToInteger(nextToken(), 1602, "Error reading number parameters")
 * @endcode
 */
int textToInteger(const char* str);

/** String (String) to size_t conversion. */
size_t textToSizeT(const char * str);

/** String (String) to integer conversion. */
inline int textToInteger(const String& str)
{
	return textToInteger(str.c_str());
}

/** String (char*) to long long integer conversion.
 *
 * @code
 * long long param_no = textToLongLong(nextToken(), 1602, "Error reading number
 *     parameters")
 * @endcode
 */
long long textToLongLong(const char* str);

/** Float to string conversion.
 *
 * If precision==0 the precision is automatically computed in such a way that
 * the number fits the width (the exponential format might be chosen). If
 * precision==-1 then the exponential format is forced. If width==0 then the
 * minimum width is used.
 *
 * @code
 * REPORT_ERROR(ERR_VALUE_INCORRECT, "Value not recognised " + floatToString(val));
 * @endcode
 */
String floatToString(float F, int _width = 8, int _prec = 0);

/** Integer to string conversion.
 *
 * If width==0 then writes the number with the number of digits needed. The
 * fill_with field indicates which is the filling character for the left
 * positions.
 *
 * @code
 * REPORT_ERROR(ERR_VALUE_INCORRECT, "Error reading key " + integerToString(key));
 * @endcode
 */
String integerToString(int I, int _width = 0, char fill_with = '0');

/** Character to integer conversion.
 *
 * Takes a character and produces a number according to its ASCII code minus 48.
 * For instance, ASCII=48 produces number 0, ASCII=49 produces 1, ..., ASCII=57
 * produces 9, ASCII=58 produces 10!!, ... This is used when you have codified
 * numbers greater than 9 in a single character.
 *
 * @code
 * int param_no = textToInt(token, 1602, "Error reading number parameters");
 * @endcode
 */
 int charToInt(const char* str);

/** String to string with given length conversion.
 *
 * The output string will have the information of the input one with the given
 * width. If the width is smaller than the string length then the string is
 * truncated and if it is greater the string is right padded with spaces. If
 * width==0 then the same string is returned.
 */
String stringToString(const String& str, size_t _width = 0);

/** Check angle.
 *
 * If the argument is not "rot", "tilt" nor "psi" an exception is thrown
 */
void checkAngle(const String& str);

/** To lower.
 *
 * All characters between A-Z are brought to a-z. Result is rewritten on input
 * string
 */
void toLower(char* _str);

/** To lower, for STL strings.
 */
void toLower(String& _str);

/** Removes white spaces from the beginning and the end of the string
as well as escaped characters
and simplifies the rest of groups of white spaces of the string to
a single white space */
String simplify( const String& str );

/** Remove trailing spaces */
void trim(String& str);

/** Remove consecutive spaces.
 *
 * All consecutive spaces are replaced by a single one and starting and
 * finishing spaces are removed
 */
String removeSpaces(const String& _str);

/** Remove quotes.
 *
 * This function removes the first character if it is a double or single quote,
 * as well as the last character. The char pointer might be moved.
 *
 * @code
 * char str[10] = "\"Hello\"";
 * (&str);
 * @endcode
 */
void removeQuotes(char** _str);

/** Replace.
 * This function replaces in the string all the occurrences of tFind and
 * replaces it with tReplace.
 */
String findAndReplace(const String& tInput, const String &tFind,
	const String &tReplace);
//@}

/** @name Tokenization
 *
 * These functions allow to split a string into small pieces separated by blank
 * spaces, giving you a pointer to the different word each time. The different
 * elements from the string are selected using strtok, so after the application
 * of this function to the input string, this is modified and NULL characters
 * are introduced as delimiters of the elements. This is useful in most
 * situations since after reading a list you might go on reading more things,
 * but you must be aware of it. Here goes an example of doing so:
 *
 * @code
 * std::cout << "Whole  line: " << line << std::endl;
 * std::cout << "First  word: " << firstToken(line) << std::endl;
 * std::cout << "Second word: " << nextToken() << std::endl;
 * std::cout << "Third  word: " << nextToken() << std::endl;
 * ...
 * @endcode
 *
 * When there are no more words, both functions return a NULL pointer. Here we
 * make a distinction between tokens (words that might be empty) and words
 * (words that cannot be empty, if they are then an exception or an exit error
 * is thrown).
 *
 * For STL there is another way. You supply a string object and a vector of
 * strings is returned with all the elements
 */
//@{
/** Split a STL string given some delimiter.
 *
 * Returns a the number of tokens found. The tokens are in the variable results.
 */
int splitString(const String& input,
                const String& delimiter,
                StringVector & results,
                bool includeEmpties = false);

/** Returns first token (char*).
 *
 * @code
 * char line[80];
 *
 * std::cout << "First  word: " << firstToken(line) << std::endl;
 * @endcode
 */
inline char* firstToken(const char* str)
{
    return strtok((char*) str, " \t\n");
}

/** Returns first token (STL).
 *
 * @code
 * String line;
 *
 * std::cout << "First  word: " << firstToken(line) << std::endl;
 * @endcode
 */
inline char* firstToken(const String& str)
{
    return strtok((char*) str.c_str(), " \t\n");
}

/** Returns next token.
 *
 * This functions returns the next word of the line we have given last as
 * parameter to firstToken.
 *
 * @code
 * char line[80];
 * ...
 * firstToken(line);
 * std::cout << "Second  word: " << nextToken(line) << std::endl;
 *
 * stf::string line;
 * ...
 * firstToken(line);
 * std::cout << "Second  word: " << nextToken(line) << std::endl;
 * @endcode
 */
inline char* nextToken()
{
    return strtok((char*) NULL, " \t\n");
}

/** Returns next token.
 *
 * It reads from position i. Returns (in i) the following position to search on.
 * When there are no more tokens. It returns "".
 */
String nextToken(const String& str, size_t& i);

/** Get non empty string (char*).
 *
 * This function returns the first word found in the given line disregarding the
 * leading blanks. If no word is found then an exception or an exit error is
 * produced. After calling this function the first blank character after the
 * word is substituted by a NULL character (as it uses the function firstToken.
 * Further word readings should use the function read_nextWord
 */
char* firstWord(char* str);

/** Get non empty string (STL).
 *
 * Same as the previous function but for STL strings
 */
inline char* firstWord(String& str)
{
    // FIXME C-style cast
    return firstWord((char*) str.c_str());
}

/** Get next non empty string.
 *
 * This is the same as the nextToken, but an exception is thrown or an exit
 * error produced if the word is empty
 */
inline char* nextWord()
{
    return firstWord((char*) NULL);
}

/** Tokenize a string and return a list of tokens
 *
 */
void tokenize(const String& str,
              StringVector & tokens,
              const String& delimiters = " \t");

/** Tokenizer  for char arrays. Similar to strtok but does NOT modify
 * the input array
 *  src is a pointer to the array beginning.
 *  It may be modified to trim the token
 * _end is a pointer to the end of the token
 * sep is an array with the valid separator
 @code
    char   *start;
    char   *str;
    char   *end;
    start=(char *) malloc(128);
    strcpy(start,"Filtrose    el\tmozo\n en una zahurda lobrega, las paredes\
    enhollinadas");
    fprintf(stderr,"Mstart=%d",(void *) start);
    char aux[64];

    str = mystrtok(&start, &end," \t\r\n");
    strncpy(aux,str,end-str); aux[end-str]=0;
    printf(">%s< %d\n", aux,(void *) end);
    str = mystrtok(&start, &end," \t\r\n");
    strncpy(aux,str,end-str); aux[end-str]=0;
    printf(">%s< %d\n", aux,(void *) end);
    str = mystrtok(&start, &end," \t\r\n");
    strncpy(aux,str,end-str); aux[end-str]=0;
    printf(">%s< %d\n", aux,(void *) end);
    str = mystrtok(&start, &end," \t\r\n");
    strncpy(aux,str,end-str); aux[end-str]=0;
    printf(">%s< %d\n", aux,(void *) end);
 @endcode
 */
char   *memtok(char **src,  char **_end, const char *sep);

/** Memory string search, taken from GNU C Library */
void * _memmem ( const void *haystack, size_t haystack_len, const void *needle, size_t needle_len);
//@}

/** Obtain an string from a format in the way of printf works
 * Example: formatString("vectorHeader@%s",fn_out.c_str())
 */
String formatString(const char * format, ...);

/** Obtain an string from a format in the way of printf works
 *
 */
void formatStringFast(String &str, const char * format, ...);

/** True if the inputString matches the regular expression in pattern */
bool matchRegExp(const String &inputString, const String &pattern);
/** split long comments in several lines starting with # */
String WordWrap(const String &inputString, size_t lineLength);

/** Escape for regular expressions.
 * Escape characters that could be misunderstood by regular expressions
 */
String escapeForRegularExpressions(const String &str);

//@}
#endif
