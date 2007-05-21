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

#include "args.h"

#include <cstdio>
#include <cmath>
#include <cerrno>

#define GCC_VERSION (__GNUC__ * 10000 \
   + __GNUC_MINOR__ * 100 \
   + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
   #include <sstream>
#else
   #include <strstream>
#endif

/* NOTE: not a very safe implemenation but standard c functions do not retrieve
 * more than 6 significative digits */
double AtoD(const char* str, int _errno, std::string errmsg, int exit)
{
    double retval;
    int ok;

    if (str == NULL)
        if (exit)
            EXIT_ERROR(_errno, errmsg);
        else
            REPORT_ERROR(_errno, errmsg);

    ok = sscanf(str, "%lf", &retval);

    if (ok)
        return retval;

    if (exit)
        EXIT_ERROR(_errno, errmsg);
    else
        REPORT_ERROR(_errno, errmsg);

    return 0;
}

float AtoF(const char* str, int _errno, std::string errmsg, int exit)
{
    float retval;
    int ok;

    if (str == NULL)
        if (exit)
            EXIT_ERROR(_errno, errmsg);
        else
            REPORT_ERROR(_errno, errmsg);

    ok = sscanf(str, "%f", &retval);

    if (ok)
        return retval;

    if (exit)
        EXIT_ERROR(_errno,errmsg);
    else
        REPORT_ERROR(_errno,errmsg);

    return 0;
}

int AtoI(const char* str, int _errno, std::string errmsg, int exit)
{
    int retval;
    int ok;

    if (str == NULL)
        if (exit)
            EXIT_ERROR(_errno, errmsg);
        else
            REPORT_ERROR(_errno, errmsg);

    ok = sscanf(str, "%d", &retval);

    if (ok)
        return retval;

    if (exit)
        EXIT_ERROR(_errno, errmsg);
    else
        REPORT_ERROR(_errno, errmsg);

    return 0;
}

long long AtoLL(const char* str, int _errno, std::string errmsg, int exit)
{
    long long int retval;
    int ok;

    if (str == NULL)
        if (exit)
            EXIT_ERROR(_errno, errmsg);
        else
            REPORT_ERROR(_errno, errmsg);

    ok = sscanf(str, "%lld", &retval);

    if (ok)
        return retval;

    if (exit)
        EXIT_ERROR(_errno,errmsg);
    else
        REPORT_ERROR(_errno,errmsg);
    return 0;
}

int best_prec(float F, int _width)
{
    // If it is 0
    if (F == 0)
        return 1;

    // Otherwise
    int exp = FLOOR(log10(ABS(F)));
    int advised_prec;

    if (exp >= 0)
        if (exp > _width-3)
            advised_prec = -1;
        else
            advised_prec = _width-2;
    else
    {
        advised_prec = _width + (exp-1) - 3;
        if (advised_prec <= 0)
            advised_prec = -1;
    }

    if (advised_prec < 0)
        advised_prec = -1; // Choose exponential format

    return advised_prec;
}

std::string FtoA(float F, int _width, int _prec)
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
        _prec = best_prec(F, _width);

    if (_prec == -1 && _width > 7)
    {
        outs.precision(_width - 7);
        outs.setf(std::ios::scientific);
    }
    else
        outs.precision(_prec);

#if GCC_VERSION < 30301
    outs << F << ends;
#else
    outs << F;
#endif

#if GCC_VERSION < 30300
    return std::string(aux);
#else
    std::string retval = outs.str();
    int i = retval.find('\0');

    if (i != -1)
        retval = retval.substr(0, i);

    return retval;
#endif
}

std::string ItoA(int I, int _width, char fill_with)
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
    for (int i=0; i < width; i++)
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
        return static_cast< std::string >("-")  + aux;
    else
        return static_cast< std::string >(aux);
}

int CtoI(const char* str, int _errno, std::string errmsg, int exit)
{
    char readval;
    int ok;

    if (str == NULL)
        if (exit)
            EXIT_ERROR(_errno, errmsg);
        else
            REPORT_ERROR(_errno, errmsg);

    ok = sscanf(str, "%c", &readval);

    if (ok)
        return readval - 48;

    if (exit)
        EXIT_ERROR(_errno, errmsg);
    else
        REPORT_ERROR(_errno, errmsg);

    return 0;
}

string AtoA(const std::string& str, int _width)
{
    if (_width == 0)
        return str;

    if (_width < str.length())
        return str.substr(0, _width);

    std::string aux = str;
    return aux.append(_width - str.length(), ' ');
}

void check_angle_descr(const std::string& str)
{
    if (str == "rot")
        return;

    if (str == "tilt")
        return;

    if (str == "psi")
        return;

    REPORT_ERROR(1,
        static_cast< std::string >(
        "check_angle_descr: Not recognized angle type: " + str));
}

std::string remove_spaces(const std::string& _str)
{
    std::string retval;
    int first = _str.find_first_not_of("\n \t");
    int last = _str.find_last_not_of("\n \t");
    bool after_blank = false;
    int imax = _str.length();

    for (int i=first; i<=last; i++)
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
void remove_quotes(char **_str)
{
    string retval=*_str;
    if (retval.length()==0)
        return;
    char c=retval[0];
    if (c=='\"' || c=='\'')
        retval=retval.substr(1,retval.length()-1);
    c=retval[retval.length()-1];
    if (c=='\"' || c=='\'')
        retval=retval.substr(0,retval.length()-1);
    free(*_str);
    *_str=strdup(retval.c_str());
}

// Split a string ==========================================================
int splitString(const std::string& input,
                const std::string& delimiter,
                std::vector< std::string >& results,
                bool includeEmpties)
{
    results.clear();
    int iPos = 0;
    int newPos = -1;
    int sizeS2 = static_cast< int >(delimiter.size());
    int isize = static_cast< int >(input.size());

    if (isize == 0 || sizeS2 == 0)
        return 0;

    std::vector< int > positions;
    newPos = input.find(delimiter, 0);

    if (newPos < 0)
        return 0;

    int numFound = 0;
    while (newPos >= iPos)
    {
        numFound++;
        positions.push_back(newPos);
        iPos = newPos;
        newPos = input.find(delimiter, iPos + sizeS2);
    }

    if (numFound == 0)
        return 0;

    for (int i=0; i <= static_cast< int >(positions.size()); i++)
    {
        string s("");
        if (i==0)
            s=input.substr(i, positions[i]);
        int offset = positions[i-1] + sizeS2;
        if (offset<isize)
        {
            if (i==positions.size())
                s=input.substr(offset);
            else if (i>0)
                s=input.substr(positions[i-1]+sizeS2,
                               positions[i]-positions[i-1]-sizeS2);
        }
        if (includeEmpties || s.size()>0)
            results.push_back(s);
    }
    return numFound;
}

// To lower ================================================================
void tolower(char *_str)
{
    int i=0;
    while (_str[i]!='\0')
    {
        if (_str[i]>='A' && _str[i]<='Z')
            _str[i] += 'a'-'A';
        i++;
    }
}

void tolower(string &_str)
{
    int i=0;
    while (_str[i]!='\0')
    {
        if (_str[i]>='A' && _str[i]<='Z')
            _str[i] += 'a'-'A';
        i++;
    }
}

// Next token ==============================================================
string next_token(const string &str, int &i)
{
    string retval;
    if (i>=str.length())
        return retval;
    int j=str.find_first_not_of(" \n",i);
    if (j==-1)
        return retval;
    int k=str.find_first_of(" \n",j+1);
    if (k==-1)
        k=str.length();
    retval=str.substr(j,k-j+1);
    i=k+1;
    return retval;
}

// Get word ================================================================
char *first_word(char *str, int _errno, string errmsg, int exit)
{
    char *token;

    // Get token
    if (str!=NULL)
        token = first_token(str);
    else
        token = next_token();

    // Check that there is something
    if (token==NULL)
        if (exit)
            EXIT_ERROR(_errno, errmsg);
        else
            REPORT_ERROR(_errno, errmsg);

    return token;
}

// Tokenize a C++ string ===================================================
void tokenize(const string& str, vector<string>& tokens,
              const string& delimiters)
{
    tokens.clear();
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

// Get parameters from the command line ====================================
char *get_param(int argc, char **argv, const char *param,
                const char *option, int _errno, string errmsg, int exit)
{
    int i = 0;

    while ((i < argc) && (strcmp(param, argv[i])))
        i++;
    if (i < argc-1)
        return(argv[i+1]);
    else
        if (option == NULL)
        {
            if (_errno==-1)
            {
                _errno=2104;
                errmsg=(string)"Argument "+param+" not found or invalid argument";
            }
            if (exit)
                EXIT_ERROR(_errno,errmsg);
            else
                REPORT_ERROR(_errno,errmsg);
        }

    return((char *) option);
}

// Get 2 parameters ========================================================
bool get_2_double_params(int argc, char **argv, const char *param,
                         double &v1, double &v2, double v1_def, double v2_def,
                         int _errno, string errmsg, int exit)
{
    bool retval;
    int i=position_param(argc,argv,param);
    if (i!=-1)
    {
        if (i+2>=argc)
        {
            string msg;
            if (errmsg=="")
                msg=(string)"Not enough arguments after "+*param;
            else
                msg=errmsg;
            if (exit)
                EXIT_ERROR(_errno,errmsg);
            else
                REPORT_ERROR(_errno,errmsg);
        }
        v1=AtoF(argv[i+1]);
        v2=AtoF(argv[i+2]);
        retval=true;
    }
    else
    {
        v1=v1_def;
        v2=v2_def;
        retval=false;
    }
    return retval;
}

// Get 3 parameters ========================================================
bool get_3_double_params(int argc, char **argv, const char *param,
                         double &v1, double &v2, double &v3,
                         double v1_def, double v2_def, double v3_def,
                         int _errno, string errmsg, int exit)
{
    bool retval;
    int i=position_param(argc,argv,param);
    if (i!=-1)
    {
        if (i+3>=argc)
        {
            string msg;
            if (errmsg=="")
                msg=(string)"Not enough arguments after "+*param;
            else
                msg=errmsg;
            if (exit)
                EXIT_ERROR(_errno,errmsg);
            else
                REPORT_ERROR(_errno,errmsg);
        }
        v1=AtoF(argv[i+1]);
        v2=AtoF(argv[i+2]);
        v3=AtoF(argv[i+3]);
        retval=true;
    }
    else
    {
        v1=v1_def;
        v2=v2_def;
        v3=v3_def;
        retval=false;
    }
    return retval;
}

// Checks if a boolean parameter was included the command line =============
bool check_param(int argc, char **argv, const char *param)
{
    int i = 0;

    while ((i < argc) && (strcmp(param, argv[i])!=0))
        i++;

    if (i < argc)
        return(true);
    else
        return(false);
}

// Position of a parameter in the command line =============================
int position_param(int argc, char **argv, const char *param)
{
    int i = 0;

    while ((i < argc) && (strcmp(param, argv[i])))
        i++;

    if (i < argc-1)
        return  i;
    else
        return -1;
}

// Number of components ====================================================
int component_no(const string &str)
{
    int imax=str.length();
    int retval=0;
    if (str[0]!='[' && str[imax-1]!=']')
        return retval;
    for (int i=0; i<imax; i++)
        if (str[i]==',')
            retval++;
    return retval+1;
}

// Get float vector ========================================================
matrix1D<double> get_vector_param(int argc, char **argv, const char *param,
                                  int dim, int _errno,
                                  string errmsg,
                                  int exit)
{
    matrix1D<double> aux;
    bool count_dimensionality=(dim==-1);

    // Find and form vector
    int pos=position_param(argc,argv,param);
    if (pos==-1 || pos+1==argc)
        if (count_dimensionality)
            return aux;
        else
        {
            if (_errno==-1)
            {
                _errno=2104;
                errmsg=(string)"Argument "+param+" not found or invalid argument";
            }
            if (exit)
                EXIT_ERROR(_errno,errmsg);
            else
                REPORT_ERROR(_errno,errmsg);
        }
    pos++;
    if (*(argv[pos])!='[')
    {
        try
        {
            double d=AtoF(argv[pos]);
            aux.resize(1);
            aux(0)=d;
            return aux;
        }
        catch (Xmipp_error XE)
        {
            if (_errno==-1)
            {
                _errno=2104;
                errmsg=(string)"Argument "+param+" not found or invalid argument";
            }
            if (exit)
                EXIT_ERROR(_errno,errmsg);
            else
                REPORT_ERROR(_errno,errmsg);
        }
    }

    string vector;
    bool finished=false;
    while (!finished)
    {
        vector += argv[pos];
        if (vector[vector.length()-1]==']')
            finished=true;
        if (++pos==argc && !finished)
        {
            if (_errno==-1)
            {
                _errno=2104;
                errmsg=(string)"Argument "+param+" not found or invalid argument";
            }
            if (exit)
                EXIT_ERROR(_errno,errmsg);
            else
                REPORT_ERROR(_errno,errmsg);
        }
    }

    // Remove brackets
    vector = vector.substr(1,vector.length()-2);

    // Count dimensionality
    int start_copy=0, end_copy;
    if (count_dimensionality)
    {
        dim=0;
        start_copy=0;
        do
        {
            end_copy=vector.find(',',start_copy);
            if (end_copy==-1)
                break;
            start_copy=end_copy+1;
            dim++;
        }
        while (1);
        dim++;
    }

    // Read diferent vector elements
    int i=0;
    start_copy=0;
    aux.resize(dim);
    while (i<dim-1)
    {
        // Find colon
        end_copy=vector.find(',',start_copy);
        // Store number
        aux(i)=AtoF(vector.substr(start_copy,end_copy));

        // Prepare for next iteration
        i++;
        start_copy=end_copy+1;
    }

    // Copy last element
    aux(i)=AtoF(vector.substr(start_copy,vector.length()));

    return aux;
}

// Specific command line ===================================================
void specific_command_line(const string &prog_name, int argc, char **argv,
                           int &argcp, char ***argvp)
{
    int i=1;
    // Look for the program name
    while (i<argc)
        if (prog_name==argv[i])
            break;
        else
            i++;
    // If not found, exit
    if (i==argc)
    {
        argcp=0;
        return;
    }

    // Check that following token is '('
    i++;
    if (i==argc || strcmp(argv[i],"(")!=0)
    {
        argcp=0;
        return;
    }

    // The specific command line starts here
    *argvp=&argv[i];
    argcp=1;

    // Look for ')'
    i++;
    if (i==argc)
    {
        argcp=0;
        return;
    }
    while (i<argc)
    {
        if (strcmp(argv[i],")")!=0)
        {
            i++;
            argcp++;
        }
        else
            break;
    }
    // If ')' not found exit
    if (i==argc)
    {
        argcp=0;
        return;
    }
}

// Generate command line ===================================================
#define INSIDE_WORD  1
#define OUTSIDE_WORD 2
void generate_command_line(const string &command_line, int &argcp,
                           char ** &argvp, char* &copy)
{
    int L=command_line.length();

    // Some initialization
    if (L==0)
    {
        argcp=0;
        return;
    }
    if (command_line[0]=='\n')
    {
        argcp=0;
        return;
    }

    // Check that argvp and copy are empty
    if (argvp!=NULL)
        delete argvp;
    if (copy!=NULL)
        delete copy;

    // Copy command line
    copy = new char[L+1];
    int i=0;
    while (i<L && command_line[i]!='\n')
    {
        copy[i]=command_line[i];
        i++;
    }
    L=i;
    copy[L]='\0';

    // Now count how many different words are there
    int words;
    int state;
    if (copy[0]==' ')
    {
        state=OUTSIDE_WORD;
        words=0;
    }
    else
    {
        state=INSIDE_WORD;
        words=1;
    }
    i=1;
    while (i<L)
    {
        if (state==OUTSIDE_WORD && copy[i]!=' ')
        {
            state=INSIDE_WORD;
            words++;
        }
        if (state==INSIDE_WORD && copy[i]==' ')
            state=OUTSIDE_WORD;
        i++;
    }

    // Resize argv and cut words
    argvp = new char *[words+1];
    argvp[0]=new char[6];
    strcpy(argvp[0],"autom");
    if (copy[0]==' ')
    {
        state=OUTSIDE_WORD;
        argcp=1;
    }
    else
    {
        state=INSIDE_WORD;
        argvp[1]=&(copy[0]);
        argcp=2;
    }
    i=1;
    while (i<L)
    {
        if (state==OUTSIDE_WORD && copy[i]!=' ')
        {
            state=INSIDE_WORD;
            argvp[argcp]=&(copy[i]);
            argcp++;
        }
        if (state==INSIDE_WORD && copy[i]==' ')
        {
            state=OUTSIDE_WORD;
            copy[i]='\0';
        }
        i++;
    }
}

// Generate command line from file =========================================
bool generate_command_line(FILE *fh, const char *param, int &argcp,
                           char ** &argvp, char* &copy)
{
    long actual_pos=ftell(fh);
    fseek(fh,0,SEEK_SET);

    char    line[201];
    char   *retval;
    bool found=false;

    // Read lines
    while (fgets (line, 200,fh) != NULL && !found)
    {
        if (line[0]==0)
            continue;
        if (line[0]=='#')
            continue;
        if (line[0]==';')
            continue;
        if (line[0]=='\n')
            continue;

        int i=0;
        while (line[i]!=0 && line[i]!='=')
            i++;
        if (line[i]=='=')
        {
            line[i]=0;
            if (strcmp(line,param)==0)
            {
                retval=line+i+1;
                found=true;
                break;
            }
        }
    }
    fseek(fh,actual_pos,SEEK_SET);
    if (!found)
        return false;

    string artificial_line;
    artificial_line=(string)"-"+param+" "+retval;
    generate_command_line(artificial_line, argcp, argvp, copy);
    return true;
}

// Get "parameter" from file ===============================================
string get_param(FILE *fh, const char *param, int skip, const char *option,
                 int _errno, string errmsg, int exit)
{
    long actual_pos=ftell(fh);
    fseek(fh,0,SEEK_SET);

    char    line[201];
    string  retval;
    bool found=false;
    int skipped=0;

    // Read lines
    while (fgets (line, 200,fh) != NULL && !found)
    {
        if (line[0]==0)
            continue;
        if (line[0]=='#')
            continue;
        if (line[0]==';')
            continue;
        if (line[0]=='\n')
            continue;

        int i=0;
        char *line_wo_spaces=line;
        while (*line_wo_spaces==' ' || *line_wo_spaces=='\t')
            line_wo_spaces++;
        while (line_wo_spaces[i]!=0 && line_wo_spaces[i]!='=')
            i++;
        if (line_wo_spaces[i]=='=')
        {
            line_wo_spaces[i]=0;
            if (strcmp(line_wo_spaces,param)==0)
            {
                if (skipped==skip)
                {
                    retval=line_wo_spaces+i+1;
                    found=true;
                    break;
                }
                else
                    skipped++;
            }
        }
    }
    fseek(fh,actual_pos,SEEK_SET);
    if (!found)
        if (option == NULL)
        {
            if (_errno==-1)
            {
                _errno=2104;
                errmsg=(string)"Argument "+param+" not found or invalid argument";
            }
            if (exit)
                EXIT_ERROR(_errno,errmsg);
            else
                REPORT_ERROR(_errno,errmsg);
        }
        else
            return option;
    else
        return remove_spaces(retval);
    return "";
}

// Check "parameter" from file =============================================
bool check_param(FILE *fh, const char *param)
{
    long actual_pos=ftell(fh);
    fseek(fh,0,SEEK_SET);

    char    line[201];
    bool found=false;
    string retval;

    // Read lines
    while (fgets (line, 200,fh) != NULL)
    {
        if (line[0]==0)
            continue;
        if (line[0]=='#')
            continue;
        if (line[0]=='\n')
            continue;

        int i=0;
        while (line[i]!=0 && line[i]!='=')
            i++;
        if (line[i]=='=')
        {
            line[i]=0;
            if (strcmp(line,param)==0)
            {
                retval=line+i+1;
                found=true;
                break;
            }
        }
    }
    fseek(fh,actual_pos,SEEK_SET);
    return found && retval!="No" && retval!="NO" && retval!="no";
}

// Get vector param from file ==============================================
matrix1D<double> get_vector_param(FILE *fh, const char *param,
                                  int dim, int _errno,  string errmsg, int exit)
{
    int    argcp;
    char **argvp=NULL;
    char  *copy=NULL;
    matrix1D<double> retval;
    if (!generate_command_line(fh, param, argcp, argvp, copy))
        if (dim==-1)
            return retval;
        else if (exit)
            EXIT_ERROR(_errno,errmsg);
        else
            REPORT_ERROR(_errno,errmsg);
    else
    {
        retval=get_vector_param(argcp, argvp, ((string)"-"+param).c_str(), dim,
                                _errno, errmsg, exit);
        delete copy;
        return retval;
    }
    return retval;
}
