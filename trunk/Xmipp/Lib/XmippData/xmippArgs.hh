/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/
#ifndef _XMIPPARGS_HH
   #define _XMIPPARGS_HH
#include <stdio.h>
#include <vector>
#include "xmippFuncs.hh"
#include "xmippMatrices1D.hh"

/**@name Arguments
    This set of functions are designed for make easier the checkings and
    parameter reading from a string. The single value readings don't
    modify the value of the input string, while the multiple value readings
    do. In general, the reading is oriented to be done with a chain of
    strtok (a general C function), for this reason the string itself is
    modified by the succesive applications of strtok.
    The second group (list reading) makes the use of tokens internally,
    while in the single read functions the token must be given by hand.
    Anyway, it is not so difficult to read a set of fields from the line
    with all the normal checkings we would like (existence of a parameter,
    checking for numerical correctness, ...)
    
    Here goes an example of use, the following piece of code try to read
    a number, a character and a list of numbers from a line. The example
    is taken from the \Ref{Document File} line reading. In the example
    you can also see how to write code for the different error handling
    methods.

    \begin{verbatim}
    string           line;
    int              key, param_no;
    vector<float>    data;

    #ifndef _NO_EXCEPTION
    try {
    #endif
       key=AtoF(first_token(line),1602,"Error reading key");
       param_no=CtoI(next_token(),1602,"Error reading number parameters");
       read_float_list(NULL,param_no,data,1602,"Error reading doc file line");
    #ifndef _NO_EXCEPTION
    } catch (Xmipp_error XE) {
       cout << XE;
       DL.line_type=0;
       REPORT_ERROR(1602,"Line is discarded");
    }
    #endif
   
    \end{verbatim}
*/
//@{
// Type conversions ========================================================
/**@name Type conversions.
    All these functions try to produce a value of the desired type from
    the given string (char *) or to produce a string from a value. If it is not
    possible (there is nothing in the input string,
    there's something which is not a float value), then the function throws
    an exception. There are default values for the exception thrown but
    you can provide a different exception in the parameter line. */
//@{
/** Best precision for a float number.
    This function returns the best precision to be used in a "printf"
    format if this number is to fit in a given width. It returns -1 if the
    exponential format is advised.
    \\Ex:
    \begin{verbatim}
    template <class T>
       ostream& operator << (ostream& ostrm, const mT& v) {
       int i,j;

       if (v.xdim==0)
          ostrm << "NULL matrix\n";
       else
          ostrm << endl;
          mT aux=ABSnD(v);
          int prec=best_prec(aux.max(),10);

          for (i=STARTINGY(v); i<=FINISHINGY(v); i++) {
             for (j=STARTINGX(v); j<=FINISHINGX(v); j++) {
                ostrm << FtoA((float)v(i,j),10,prec) << ' ';
             }
             ostrm << endl;
          }

       return ostrm;
    }
    \end{verbatim}
*/
int best_prec(float F, int _width);

/** String (char *) --> Double.
    \\Ex: key=AtoD(first_token(line),1602,"Error reading key"); */
double AtoD(const char *str, int _errno=2101, string errmsg="Error in AtoD",
      int exit=0);

/** String (char *) --> Float.
    \\Ex: key=AtoF(first_token(line),1602,"Error reading key"); */
float AtoF(const char *str, int _errno=2101, string errmsg="Error in AtoF",
      int exit=0);

/** String (Stl) --> Float.
    \\Ex: key=AtoF(str,1602,"Error reading key"); */
inline float AtoF(const string str, int _errno=2101, string errmsg="Error in AtoF",
      int exit=0) {return AtoF(str.c_str(),_errno,errmsg,exit);}

/** String (char *) --> Integer.
    \\Ex: param_no=AtoI(next_token(),1602,"Error reading number parameters")*/
int   AtoI(const char *str, int _errno=2102, string errmsg="Error in AtoI",
      int exit=0);

/** String (Stl) --> Integer.
    \\Ex: param_no=AtoI(str,1602,"Error reading number parameters")*/
inline int AtoI(const string str, int _errno=2102, string errmsg="Error in AtoI",
      int exit=0) {return AtoI(str.c_str(),_errno,errmsg,exit);}

/** String (char *) --> Long Long Integer.
    \\Ex: param_no=AtoLL(next_token(),1602,"Error reading number
    parameters")*/
long long AtoLL(const char *str, int _errno=2102, string errmsg="Error in AtoL",
      int exit=0);

/** Float --> String.
    If precision==0 then the precision is automatically computed
    such that the number fits the width (the exponential format
    might be chosen). If precision==-1 then the exponential format
    is forced. If width==0 then the minimum width is used.
    \\Ex: REPORT_ERROR(1602,"Value not recognised "+FtoA(val)); */
string FtoA(float F, int _width=8, int _prec=0);

/** Integer --> String.
    If width is equal 0 then it writes the number with
    the number of digits needed. The fill_with field indicates
    which is the filling character for the left positions.
    \\Ex: REPORT_ERROR(1602,"Error reading key "+ItoA(key)); */
string ItoA(int I, int _width=0, char fill_with='0');

/** Character --> Integer.
    This function takes a character and produces a number according to its
    ASCII code minus 48. For instance, ASCII=48 produces number 0, ASCII=49
    produces 1, ..., ASCII=57 produces 9, ASCII=58 produces 10!!, ...
    This is used when you have codified numbers greater than 9 in a single
    character.
    \\ Ex: param_no=CtoI(token,1602,"Error reading number parameters")*/
int   CtoI(const char *str, int _errno=2103, string errmsg="Error in CtoI",
      int exit=0);

/** String --> String with given length.
    The output string will have the information of the input one with the
    given width. If the width is smaller than the string length then the
    string is truncated and if it is greater the string is right padded
    with spaces. If width==0 then the same string is returned. */
string AtoA(const string &str, int _width=0);

/** Check angle.
    If the argument is not "rot", "tilt" nor "psi" an exception is thrown */
void check_angle_descr(const string &str);

/** To lower.
    All characters between A-Z are brought to a-z. Result is rewritten
    on input string*/
void tolower(char *_str);

/** To lower for C++ strings.*/
void tolower(string &_str);

/** Remove consecutive spaces.
    All consecutive spaces are replaced by a single one and
    starting and finishing spaces are removed*/
string remove_spaces(const string &_str);

/** Remove quotes.
    This function removes the first character if it is a double or single
    quote, as well as the last character. The char pointer might be moved.
    \\Ex: char str[10]="\"Hello\""; remove_quotes(&str); */
void remove_quotes(char **_str);

//@}

// Tokenization ============================================================
/**@name Tokenization

   These functions allow to split a string into small pieces separated by
   blank spaces, giving you a pointer to the different words each time.
   The different elements from the string are selected using strtok, so after
   the application of this function to the input string, this is modified
   and NULL characters are introduced as delimiters of the elements.
   This is useful in most situations since after reading a list you might
   go on reading more things, but you must be aware of it. Here goes an
   example of doing so:
   
   \begin{verbatim}
   cout << "Whole  line: " << line << endl;
   cout << "First  word: " << first_token(line) << endl;
   cout << "Second word: " << next_token() << endl;
   cout << "Third  word: " << next_token() << endl;
   ...
   \end{verbatim}
   
   When there aren't more words, both functions return a NULL pointer.
   Here we make a distinction between tokens (words that might be empty)
   and words (words that cannot be empty, if they are then an exception
   or an exit error is thrown).
   
   For C++, there is anothere way. You supply a string object and a vector
   of strings is returned with all the elements.*/
//@{
   /** Split a STL string given some delimiter.
       Returns a the number of tokens found. The tokens are in the variable
       results.*/
int splitString(const string& input, 
    const string& delimiter, vector<string>& results, 
    bool includeEmpties=false);

   /** Returns first token (char *).
       \\Ex: char line[80]; ...; cout << "First  word: "
             << first_token(line) << endl;*/
inline char * first_token(const char * str) {return strtok((char *)str," \n");}

   /** Returns first token (string).
       \\Ex: string line; ...; cout << "First  word: "
             << first_token(line) << endl;*/
inline char * first_token(const string &str) {return strtok((char *)str.c_str()," \n");}
   
   /** Returns next token.
       This functions returns the next word of the line we have given
       last as parameter to first_token.
       \\Ex: char line[80]; ...; first_token(line); cout << "Second  word: "
             << next_token(line) << endl;
       \\Ex: string line; ...; first_token(line); cout << "Second  word: "
             << next_token(line) << endl;*/
inline char * next_token() {return strtok((char *)NULL," \n");}

   /** Returns next token.
       It reads from position i. Returns (in i) the following position to search on.
       When there are no more tokens. It returns "".*/
       string next_token(const string &str, int &i);

   /** Get non empty string (char *).
       This function returns the first word found in the given line disregarding
       the leading blanks. If no word is found then an exception or an exit error
       is produced. After calling this function the first blank character after
       the word is substituted by a NULL character (as it uses the function
       \Ref{first_token}). Further word readings should use the function
       \Ref{read_next_word}. */
   char * first_word(char *str, int _errno=2106, 
    string errmsg="first word: String not found", int exit = 0);

   /** Get non empty string (string).
       Same as the previous function but for STL strings. */
inline char * first_word(string &str, int _errno=2106, 
       string errmsg="first word: String not found", int exit = 0)
       {return first_word((char *)str.c_str(),_errno, errmsg, exit);}

   /** Get next non empty string.
       This is the same as the next_token, but an exception is thrown
       or an exit error produced if the word is empty. */
inline char * next_word(int _errno=2106, 
       string errmsg="next word: String not found", int exit = 0)
       {return first_word((char *)NULL,_errno, errmsg, exit);}

void tokenize(const string& str,
              vector<string>& tokens,const string& delimiters = " \t");

//@}

// Read lists ==============================================================
/**@name Read lists.
   These functions try to read N values of the desired type into the
   given structure (either a STL vector of any numerical type by adding 
   the read values at the end or a matrix1D of any numerical type and then 
   the values are written at PHYSICAL positions 0 ... N-1, the matrix1D must 
   be big enough to hold all the data since it is not resized internally.
   
   If it is not possible to read all parameters an exception is thrown.
   You can provide the exception in the function call but there are
   default values. These functions are based in those for reading
   single values, so it might be possible that these other functions
   might throw other exceptions.
   
   The different elements of the list are selected using the tokenizig
   functions (different elements in the string are separated by spaces),
   so after the application of this function to the input string, this
   is modified and NULL characters are introduced as delimiters of the
   elements.  Here goes an example of doing so:

   \begin{verbatim}
       getline(in_stream,line);
       read_float_list(line,10,v1);   // Read 10 values from line in v1
       read_float_list(NULL,10,v2);   // Read NEXT!! 10 values in v2
   \end{verbatim} */
//@{

/** list --> STL vector. */
template <class T> 
void read_float_list(const char *str, int N, vector<T> &v,
   int _errno=2105, string errmsg="Error reading list",
   int exit=0) {
   T  valueF;
   char     *token;

   token=first_token(str);
   for (int i=0; i<N; i++) {
      if (token==NULL) {
         //CO: Should not report other error than the required one
	 // cout << "Read float list: Number of true parameters doesn't coincide\n";
         REPORT_ERROR(_errno,errmsg);
      }

      try {
         valueF=(T) AtoF(token);
      }
      catch (Xmipp_error) {REPORT_ERROR(_errno,errmsg);}

      v.push_back(valueF);
      if (i!=N-1) token=next_token();
   }
}

/** list --> STL vector. */
template <class T> 
void read_float_list(const string &str, int &i, int N, vector<T> &v,
   int _errno=2105, string errmsg="Error reading list",
   int exit=0) {
   T  valueF;
   string token;

   token=next_token(str,i);
   for (int j=0; j<N; j++) {
      if (token=="") {
         //CO: Should not report other error than the required one
	 // cout << "Read float list: Number of true parameters doesn't coincide\n";
         REPORT_ERROR(_errno,errmsg);
      }

      try {
         valueF=(T) AtoF(token.c_str());
      }
      catch (Xmipp_error) {REPORT_ERROR(_errno,errmsg);}

      v.push_back(valueF);
      if (j!=N-1) token=next_token(str,i);
   }
}

/** list --> matrix1D. */
template <class T> 
void read_float_list(const char *str, int N, matrix1D<T> &v,
   int _errno=2105, string errmsg="Error reading floating list",
   int exit=0) {
   T  valueF;
   char *token;

   token=first_token(str);
   for (int i=0; i<N; i++) {
      if (token==NULL) {
         //CO: Should not report other error than the required one
	 // cout << "Read float list: Number of true parameters doesn't coincide\n";
         REPORT_ERROR(_errno,errmsg);
      }

      try {
         valueF= (T) AtoF(token);
      }
      catch (Xmipp_error) {REPORT_ERROR(_errno,errmsg);}

      DIRECT_VEC_ELEM(v,i)=valueF;
      if (i!=N-1) token=next_token();
   }
}

//@}

/**@name Command line functions
   These functions help you to manage the command line parameters. */
//@{
/** Get parameters from the command line.
    This function assumes that the command line is structured in such a
    way that for each parameter a block like "-param <param_value>" is
    defined. The label "param" can be substituted by any other one you
    like. If the parameter is optional then this function allows you
    to define a default value. If no deafult value is supplied and the
    parameter is not specified in the command line, then an exception
    is thrown. You may change the default exception. You may also
    indicate that in case of error no exception is raised and force the
    program to abort (use the exit variable).
    \\Ex:
    \begin{verbatim}
    m_param = AtoF(get_param(argc, argv, "-m"));
        // Get compulsory parameter "-m"
    m_param = AtoF(get_param(argc, argv, "-m","2.65"));
        // Optional parameter, if no parameter is given it takes 2.65
        // by default
    m_param = AtoF(get_param(argc, argv, "-m", NULL, 6001,
        "-m parameter not found. I'm going out",TRUE);
        // Compulsory parameter, if not found give an special
        // error message and exit the program
    \end{verbatim}*/
char *get_param(int argc, char **argv, const char *param,
                const char *option = NULL, 
	        int _errno=-1, 
		string errmsg="",
		int exit = 0);

/** Get 2 float parameters after a flag from the command line.
    An exception is thrown if there are not enough parameters after the flag,
    if the message is empty then "Not enough parameters after <param>" is shown.
    The default values must be provided. TRUE is returned if the two values
    have been found.*/
bool get_2_double_params(int argc, char **argv, const char *param,
   double &v1, double &v2, double v1_def, double v2_def,
   int _errno=2104, 
   string errmsg="",
   int exit = 0);

/** Get 3 float parameters after a flag from the command line.
    An exception is thrown if there are not enough parameters after the flag,
    if the message is empty then "Not enough parameters after <param>" is shown.
    The default values must be provided. TRUE is returned if the two values
    have been found.*/
bool get_3_double_params(int argc, char **argv, const char *param,
   double &v1, double &v2, double &v3,
   double v1_def, double v2_def, double v3_def,
   int _errno=2104, 
   string errmsg="",
   int exit = 0);

/** Get boolean parameters from the command line.
    This function assumes that the command line is structured in such a
    way that for each parameter a block like "-param" is
    defined. The label "param" can be substituted by any other one you
    like. It might be used to look for a boolean parameter, for instance:
    -verbose means that verbose functionality is set (TRUE).
    
    \\Ex:
    \begin{verbatim}
    verbose = check_param(argc, argv, "-verbose"));
        // checks if "-verbose" was supplied in the command line. If -verbose
	was supplied the function returns TRUE (1), otherwise returns FALSE (0) 
    \end{verbatim}*/
bool check_param(int argc, char **argv, const char *param);

/** Returns the position where the given parameter is in the command line.
    This function assumes that the command line is structured in such a
    way that for each parameter a block like "-param" is
    defined. The label "param" can be substituted by any other one you
    like. It returns -1 if the parameter is not found.
    It is used to look for parameters having a list of values behind,
    for instance:
    
    -ang rot tilt psi
    
    \\Ex:
    \begin{verbatim}
    i = position_param(argc, argv, "-ang"));
    // This condition checks if 3 arguments where introduced after -ang
    // parameter (assuming that the -ang argument is the last one in the string)
    if (i+3>=argc) EXIT_ERROR(1,"Not enough parameters behind -ang");
    ang1=argv[i+1];
    ang2=argv[i+2];
    ang3=argv[i+3];
    \end{verbatim}*/
int position_param(int argc, char **argv, const char *param);


/** Return the number of components of a vector argument.
    A vector argument is defined as [x,y,z,...].
    It returns 0 if the string does not contain a vector */
int component_no(const string &str);



/** Get float vector.
    A vector is defined as a "[x,y,z, ...]" of any dimension (by default
    2D vectors are read). The vector must not contain blank spaces.
    \\Ex: a=get_vector_param(argc,argv,"-a",3); a 3D vector called "-a".
    The vector is internally resized properly. If the dimension is -1 then
    all the vector components are read disregarding their dimensionality, ie,
    -1 is used to read vectors of an unknown dimension. If the parameter is
    not found when no dimensionality is given an empty vector is returned
    but no exception is thrown. If there is no dimensionality and a single
    parameter is behind the flag then no brackets are needed. */
matrix1D<double> get_vector_param(int argc, char **argv, const char *param,
   int dim=2, int _errno=-1, 
   string errmsg="",
   int exit = 0);

/** Get specific command line.
    This function allows to pass command line parameters to several
    programs which are integrated in a single program. It assumes
    that the command line is a compound line like this
    \begin{verbatim}
    whole_process prog1 ( -i in_file1 -o out_file1 ) prog2 ( -i out_file1 )
    \end{verbatim}
    In this example, a single program "whole_process" is compund of 2
    programs ("prog1" and "prog2"). The parameters for prog1 are within
    the brackets after prog1 (notice that spaces around brackets are very
    important), and then the output of program 1 is the input of
    program 2.
    
    This function returns new argc and argv such that prog1 parameters
    could be read from this new command line.
    \\Ex:
    \begin{verbatim}
    int argcp;
    char **argvp;
    specific_command_line("prog1",argc,argv,argcp,&argvp);
    if (argcp==0)
       EXIT_ERROR(1, "No parameters for prog1");
    prog1_parameters.read(argcp,argvp);
    \end{verbatim}
    
    If the program name is not found or the corresponding brackets then
    argcp==0 */
void specific_command_line(const string &prog_name, int argc, char **argv,
    int &argcp, char ***argvp);

/** Generate argc and argv for a string.
    Given a string this function makes a copy of the string and divides
    it into tokens such that they can be used as argc and argv, as
    if it were a command line. The string copy remains in "copy"
    and it can be freed by disposing this variable.
    argvp[0] (normally the program name) is set to any value,
    in this case to "autom", standing fro "automatically generated".
    
    argcp==0 if now valid command line is provided, ie, the line is
    empty or only with blanks.

    Next time the function is called it checks that argv and copy
    are empty (pointing to NULL), if they aren't then firstly the
    associated memory is freed.
    \\Ex:
    \begin{verbatim}
    int argcp;
    char **argvp, *copy;
    copy=NULL; argvp=NULL;
    string command_line="-i input_file -o output_file";
    generate_command_line(command_line,argcp,&argvp,&copy)
    if (argcp!=0)
       read_parameters(argcp,argvp);
    \end{verbatim} */
void generate_command_line(const string &command_line, int &argcp,
    char ** &argvp, char* &copy);

/** Generate articial command line from a file.
    The copy variable must be destroyed outside by "delete copy".
    This function takes "input_file=<input_file>" and turns it into
    "-input_file <input_file>" The appropiate argc, argv are also returned.
    
    Returns TRUE if the parameter is found in the file, and FALSE if it is
    not.*/
    bool generate_command_line(FILE *fh, const char *param, int &argcp,
       char ** &argvp, char* &copy);

/** Get parameter from file.
    Parameters are supposed to be identified with an =, so any line
    which doesn't contain an = character cannot contain a parameter.
    The left part of the "=" is the identifier, and the right part the
    value. The value is returned without any extra space, if it is 
    compound of several words, then the spaces in between are simplified
    to a single blank space.
    The file position inside the file is not moved and comments are
    allowed starting with "#" and ";".
    
    Parameter skip controls the number of these parameters to skip before
    returning the value, ie, if there are several "parameter=" tags in
    a file and you want the first one then you should skip 0, if you
    want the second the skip=1, ...
    
    The meaning and use of the exit, errors and optional value is the same
    as in the command line get_param.*/
    string get_param(FILE *fh, const char *param, int skip=0,
       const char *option=NULL, int _errno=-1, 
       string errmsg="",
       int exit = 0);

/** Check if a parameter is present in a file.
    The same as the previous function, but this function only reports
    if a parameter is present or not in a file. Notice that boolean
    parameters must be defined as "parameter=". If after the parameter, "no"
    comes then this function returns FALSE*/
    bool check_param(FILE *fh, const char *param);

/** Get float vector from a file.
    The same as before but reading is done from a file. */
matrix1D<double> get_vector_param(FILE *fh, const char *param,
   int dim=2, int _errno=-1, 
   string errmsg="",
   int exit = 0);

//@}
//@}

#endif
