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
#include "../xmippArgs.hh"
#include <stdio.h>
#include <math.h>

#define GCC_VERSION (__GNUC__ * 10000 \
   + __GNUC_MINOR__ * 100 \
   + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
   #include <sstream>
#else
   #include <strstream.h>
#endif

// Char --> Float ========================================================
float AtoF(const char *str, int _errno, string errmsg, int exit) _THROW {
   float retval;
   int   ok;
   if (str==NULL)
      if (exit) EXIT_ERROR(_errno,errmsg);
      else      REPORT_ERROR(_errno,errmsg);
   ok=sscanf(str,"%f",&retval);
   if (ok) return retval;
   if (exit) EXIT_ERROR(_errno,errmsg);
   else      REPORT_ERROR(_errno,errmsg);
}

// Char --> Integer ========================================================
int AtoI(const char *str, int _errno, string errmsg, int exit) _THROW {
   int retval;
   int   ok;
   if (str==NULL)
      if (exit) EXIT_ERROR(_errno,errmsg);
      else      REPORT_ERROR(_errno,errmsg);
   ok=sscanf(str,"%d",&retval);
   if (ok) return retval;
   if (exit) EXIT_ERROR(_errno,errmsg);
   else      REPORT_ERROR(_errno,errmsg);
}

// Char --> Long long Integer ==============================================
long long AtoLL(const char *str, int _errno, string errmsg, int exit) _THROW {
   long long int retval;
   int   ok;
   if (str==NULL)
      if (exit) EXIT_ERROR(_errno,errmsg);
      else      REPORT_ERROR(_errno,errmsg);
   ok=sscanf(str,"%lld",&retval);
   if (ok) return retval;
   if (exit) EXIT_ERROR(_errno,errmsg);
   else      REPORT_ERROR(_errno,errmsg);
}

// Best precission =========================================================
int best_prec(float F, int _width) {
   int exp=FLOOR(log10(ABS(F)));
   int advised_prec;
   if (exp>=0)
      if (exp>_width-3) advised_prec=-1;
      else              advised_prec=_width-2;
   else {
      advised_prec=_width+(exp-1)-3;
      if (advised_prec<=0) advised_prec=-1;
   }
      
   if (advised_prec<0) advised_prec=-1; // Choose exponential format
   return advised_prec;
}

// Float --> String ========================================================
string FtoA(float F, int _width, int _prec) {
   #if GCC_VERSION < 30300
      char aux[15];
      ostrstream outs(aux,sizeof(aux));
   #else
      ostringstream outs;
   #endif
   outs.fill(' ');
   if (_width!=0) outs.width(_width);
   if (_prec==0) _prec=best_prec(F,_width);
   if (_prec==-1 && _width>7)
      {outs.precision(_width-7); outs.setf(ios::scientific);}
   else
      outs.precision(_prec);
   #if GCC_VERSION < 30301
      outs << F << ends;
   #else
      outs << F;
   #endif 
   #if GCC_VERSION < 30300
      return (string)aux;
   #else
      string retval=outs.str();
      int i=retval.find('\0');
      if (i!=-1) retval=retval.substr(0,i);
      return retval;
   #endif
}

// Integer --> String ======================================================
string ItoA(int I, int _width, char fill_with) {
   char aux[15];
   
   // Check width
   int width=_width;
   int Iaux=I;
   if (width==0)
      do {
         Iaux/=10;
	 width++;
      } while (Iaux!=0);
   
   // Fill the number with the fill character
   for (int i=0; i<width; i++) aux[i]=fill_with;

   // Start filling the array
   aux[width--]='\0';
   Iaux=I;
   do {
      int digit=Iaux%10;
      Iaux/=10;
      aux[width--]='0'+digit;
   } while (Iaux!=0);
   
   return (string)aux;
}

// Character --> Integer ===================================================
int CtoI(const char *str, int _errno, string errmsg, int exit) _THROW {
   char  readval;
   int   ok;
   if (str==NULL)
      if (exit) EXIT_ERROR(_errno,errmsg);
      else      REPORT_ERROR(_errno,errmsg);
   ok=sscanf(str,"%c",&readval);
   if (ok) return readval-48;
   if (exit) EXIT_ERROR(_errno,errmsg);
   else      REPORT_ERROR(_errno,errmsg);
}

// String --> String with length ===========================================
string AtoA(const string &str, int _width) {
   if (_width==0) return str;
   if (_width<str.length()) return str.substr(0,_width);
   string aux=str; return aux.append(_width-str.length(),' ');
}

// Check angle description =================================================
void check_angle_descr(const string &str) _THROW {
   if (str=="rot")  return;
   if (str=="tilt") return;
   if (str=="psi")  return;
   REPORT_ERROR(1,(string)"check_angle_descr: Not recognized angle type: "+str);
}

// Remove spaces ===========================================================
string remove_spaces(const string &_str) {
   string retval;
   int first=_str.find_first_not_of("\n \t");
   int last=_str.find_last_not_of("\n \t");
   bool after_blank=FALSE;
   int imax=_str.length();
   for (int i=first; i<=last; i++) {
      if (_str[i]==' ' || _str[i]=='\n' || _str[i]=='\t') {
         if (!after_blank) retval += _str[i];
         after_blank=TRUE;
      }
      else {
         retval += _str[i];
         after_blank=FALSE;
      }
   }
   return retval;
}

// Remove quotes ===========================================================
void remove_quotes(char **_str) {
   string retval=*_str;
   if (retval.length()==0) return;
   char c=retval[0];
   if (c=='\"' || c=='\'') retval=retval.substr(1,retval.length()-1);
   c=retval[retval.length()-1];
   if (c=='\"' || c=='\'') retval=retval.substr(0,retval.length()-1);
   free(*_str);
   *_str=strdup(retval.c_str());
}

// To lower ================================================================
void tolower(char *_str) {
   int i=0;
   while (_str[i]!='\0') {
      if (_str[i]>='A' && _str[i]<='Z') _str[i] += 'a'-'A';
      i++;
   }
}

// Next token ==============================================================
string next_token(const string &str, int &i) {
   string retval;
   if (i>=str.length()) return retval;
   int j=str.find_first_not_of(" \n",i);
   if (j==-1) return retval;
   int k=str.find_first_of(" \n",j+1);
   if (k==-1) k=str.length();
   retval=str.substr(j,k-j+1);
   i=k+1;
   return retval;
}

// Get word ================================================================
char *first_word(char *str, int _errno, string errmsg, int exit) _THROW {
   char *token;

// Get token
   if (str!=NULL) token = first_token(str);
   else           token = next_token();

// Check that there is something
   if (token==NULL)
      if (exit) EXIT_ERROR(_errno, errmsg);
      else      REPORT_ERROR(_errno, errmsg);

   return token;
}

// Tokenize a C++ string ===================================================
void tokenize(const string& str, vector<string>& tokens,
   const string& delimiters) {
   tokens.clear();
   // Skip delimiters at beginning.
   string::size_type lastPos = str.find_first_not_of(delimiters, 0);
   // Find first "non-delimiter".
   string::size_type pos     = str.find_first_of(delimiters, lastPos);

   while (string::npos != pos || string::npos != lastPos) {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
   }
}

// Float list --> Vector ===================================================
template <class T> 
void read_float_list(const char *str, int N, vector<T> &v,
   int _errno, string errmsg, int exit) _THROW {
   T  valueF;
   char     *token;

   token=first_token(str);
   for (int i=0; i<N; i++) {
      if (token==NULL) {
         //CO: Should not report other error than the required one
	 // cout << "Read float list: Number of true parameters doesn't coincide\n";
         REPORT_ERROR(_errno,errmsg);
      }

      #ifndef _NO_EXCEPTION
         try {
      #endif
            valueF=(T) AtoF(token);
      #ifndef _NO_EXCEPTION
         }
         catch (Xmipp_error) {REPORT_ERROR(_errno,errmsg);}
      #endif

      v.push_back(valueF);
      if (i!=N-1) token=next_token();
   }
}

template <class T> 
void read_float_list(const string &str, int &i, int N, vector<T> &v,
   int _errno=2105, string errmsg="Error reading list",
   int exit=0) _THROW {
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

// list --> Matrix1D =================================================
template <class T> 
void read_float_list(char *str, int N, matrix1D<T> &v,
   int _errno, string errmsg, int exit) _THROW {
   T  valueF;
   char *token;

   token=first_token(str);
   for (int i=0; i<N; i++) {
      if (token==NULL) {
         //CO: Should not report other error than the required one
	 // cout << "Read float list: Number of true parameters doesn't coincide\n";
         REPORT_ERROR(_errno,errmsg);
      }

      #ifndef _NO_EXCEPTION
         try {
      #endif
            valueF= (T) AtoF(token);
      #ifndef _NO_EXCEPTION
         }
         catch (Xmipp_error) {REPORT_ERROR(_errno,errmsg);}
      #endif

      DIRECT_VEC_ELEM(v,i)=valueF;
      if (i!=N-1) token=next_token();
   }
}



// Get parameters from the command line ====================================
char *get_param(int argc, char **argv, const char *param,
   const char *option, int _errno, string errmsg, int exit) _THROW {
  int i = 0;

  while ((i < argc) && (strcmp(param, argv[i])))
    i++;
  if (i < argc-1) return(argv[i+1]);
  else
    if (option == NULL) {
      if (_errno==-1) {
         _errno=2104;
         errmsg=(string)"Argument "+param+" not found or invalid argument";
      }
      if (exit) EXIT_ERROR(_errno,errmsg);
      else      REPORT_ERROR(_errno,errmsg);
    }

  return((char *) option);
}

// Get 2 parameters ========================================================
bool get_2_double_params(int argc, char **argv, const char *param,
   double &v1, double &v2, double v1_def, double v2_def,
   int _errno, string errmsg, int exit) _THROW {
   bool retval;
   int i=position_param(argc,argv,param);
   if (i!=-1) {
      if (i+2>=argc) {
         string msg;
         if (errmsg=="") msg=(string)"Not enough arguments after "+*param;
         else            msg=errmsg;
         if (exit) EXIT_ERROR(_errno,errmsg);
         else      REPORT_ERROR(_errno,errmsg);
      }
      v1=AtoF(argv[i+1]);
      v2=AtoF(argv[i+2]);
      retval=TRUE;
   } else {
      v1=v1_def;
      v2=v2_def;
      retval=FALSE;
   }
   return retval;
}

// Get 3 parameters ========================================================
bool get_3_double_params(int argc, char **argv, const char *param,
   double &v1, double &v2, double &v3, 
   double v1_def, double v2_def, double v3_def,
   int _errno, string errmsg, int exit) _THROW {
   bool retval;
   int i=position_param(argc,argv,param);
   if (i!=-1) {
      if (i+3>=argc) {
         string msg;
         if (errmsg=="") msg=(string)"Not enough arguments after "+*param;
         else            msg=errmsg;
         if (exit) EXIT_ERROR(_errno,errmsg);
         else      REPORT_ERROR(_errno,errmsg);
      }
      v1=AtoF(argv[i+1]);
      v2=AtoF(argv[i+2]);
      v3=AtoF(argv[i+3]);
      retval=TRUE;
   } else {
      v1=v1_def;
      v2=v2_def;
      v3=v3_def;
      retval=FALSE;
   }
   return retval;
}

// Checks if a boolean parameter was included the command line =============
int check_param(int argc, char **argv, const char *param)
{
  int i = 0;
 
  while ((i < argc) && (strcmp(param, argv[i])!=0))
    i++;
 
  if (i < argc)  return(TRUE);
  else           return(FALSE);
}

// Position of a parameter in the command line =============================
int position_param(int argc, char **argv, const char *param) {
  int i = 0;
 
  while ((i < argc) && (strcmp(param, argv[i])))
    i++;
 
  if (i < argc-1) return  i;
  else            return -1;
}

// Number of components ====================================================
int component_no(const string &str) {
   int imax=str.length();
   int retval=0;
   if (str[0]!='[' && str[imax-1]!=']') return retval;
   for (int i=0; i<imax; i++)
       if (str[i]==',') retval++;
   return retval+1;
}

// Get float vector ========================================================
matrix1D<double> get_vector_param(int argc, char **argv, const char *param,
   int dim, int _errno, 
   string errmsg,
   int exit) {
   matrix1D<double> aux;
   bool count_dimensionality=(dim==-1);

   // Find and form vector
   int pos=position_param(argc,argv,param);
   if (pos==-1 || pos+1==argc)
      if (count_dimensionality) return aux;
      else {
         if (_errno==-1) {
            _errno=2104;
            errmsg=(string)"Argument "+param+" not found or invalid argument";
         }
         if (exit) EXIT_ERROR(_errno,errmsg);
         else      REPORT_ERROR(_errno,errmsg);
      }
   pos++;
   if (*(argv[pos])!='[') {
      try {
         double d=AtoF(argv[pos]);
         aux.resize(1); aux(0)=d;
         return aux;
      } catch (Xmipp_error XE) {
         if (_errno==-1) {
            _errno=2104;
            errmsg=(string)"Argument "+param+" not found or invalid argument";
         }
         if (exit) EXIT_ERROR(_errno,errmsg);
         else      REPORT_ERROR(_errno,errmsg);
      }
   }
      
   string vector;
   bool finished=FALSE;
   while (!finished) {
      vector += argv[pos];
      if (vector[vector.length()-1]==']') finished=TRUE;
      if (++pos==argc && !finished) {
         if (_errno==-1) {
            _errno=2104;
            errmsg=(string)"Argument "+param+" not found or invalid argument";
         }
         if (exit) EXIT_ERROR(_errno,errmsg);
         else      REPORT_ERROR(_errno,errmsg);
      }
   }

   // Remove brackets
   vector = vector.substr(1,vector.length()-2);

   // Count dimensionality
   int start_copy=0, end_copy;
   if (count_dimensionality) {
      dim=0;
      start_copy=0;
      do {
         end_copy=vector.find(',',start_copy);
         if (end_copy==-1) break;
         start_copy=end_copy+1;
         dim++;
      } while (1);
      dim++;
   }

   // Read diferent vector elements
   int i=0;
   start_copy=0;
   aux.resize(dim);
   while (i<dim-1) {
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
   int &argcp, char ***argvp) {
   int i=1;
   // Look for the program name
   while (i<argc)
      if (prog_name==argv[i]) break;
      else i++;
   // If not found, exit
   if (i==argc) {argcp=0; return;}
   
   // Check that following token is '('
   i++; if (i==argc || strcmp(argv[i],"(")!=0) {argcp=0; return;}

   // The specific command line starts here
   *argvp=&argv[i];
   argcp=1;
   
   // Look for ')'
   i++; if (i==argc){argcp=0; return;}
   while (i<argc) {
      if (strcmp(argv[i],")")!=0) {i++; argcp++;}
      else break;
   }
   // If ')' not found exit
   if (i==argc) {argcp=0; return;}
}

// Generate command line ===================================================
#define INSIDE_WORD  1
#define OUTSIDE_WORD 2
void generate_command_line(const string &command_line, int &argcp,
   char ** &argvp, char* &copy) {
   int L=command_line.length();

   // Some initialization
   if (L==0) {argcp=0; return;}
   if (command_line[0]=='\n') {argcp=0; return;}

   // Check that argvp and copy are empty
   if (argvp!=NULL) delete argvp;
   if (copy!=NULL)  delete copy;

   // Copy command line
   copy = new char[L+1];
   int i=0;
   while (i<L && command_line[i]!='\n') {
      copy[i]=command_line[i];
      i++;
   }
   L=i; copy[L]='\0';
   
   // Now count how many different words are there
   int words;
   int state;
   if (copy[0]==' ') {state=OUTSIDE_WORD; words=0;}
   else              {state=INSIDE_WORD;  words=1;}
   i=1;
   while (i<L) {
      if (state==OUTSIDE_WORD && copy[i]!=' ') {
         state=INSIDE_WORD; words++;
      }
      if (state==INSIDE_WORD && copy[i]==' ')
         state=OUTSIDE_WORD;
      i++;
   }
   
   // Resize argv and cut words
   argvp = new char *[words+1];
   argvp[0]=new char[6]; strcpy(argvp[0],"autom");
   if (copy[0]==' ') {state=OUTSIDE_WORD;                      argcp=1;}
   else              {state=INSIDE_WORD;  argvp[1]=&(copy[0]); argcp=2;}
   i=1;
   while (i<L) {
      if (state==OUTSIDE_WORD && copy[i]!=' ') {
         state=INSIDE_WORD; argvp[argcp]=&(copy[i]); argcp++;
      }
      if (state==INSIDE_WORD && copy[i]==' ') {
         state=OUTSIDE_WORD; copy[i]='\0';
      }
      i++;
   }
}

// Generate command line from file =========================================
bool generate_command_line(FILE *fh, const char *param, int &argcp,
   char ** &argvp, char* &copy) {
   long actual_pos=ftell(fh);
   fseek(fh,0,SEEK_SET);
   
   char    line[201];
   char   *retval;
   bool found=FALSE;

   // Read lines
   while (fgets (line, 200,fh) != NULL && !found) {
      if (line[0]==0)    continue;
      if (line[0]=='#')  continue;
      if (line[0]==';')  continue;
      if (line[0]=='\n') continue;
      
      int i=0; while (line[i]!=0 && line[i]!='=') i++;
      if (line[i]=='=') {
         line[i]=0;
         if (strcmp(line,param)==0) {
            retval=line+i+1; found=TRUE; break;
         }
      }
   }
   fseek(fh,actual_pos,SEEK_SET);
   if (!found) return FALSE;

   string artificial_line;
   artificial_line=(string)"-"+param+" "+retval;
   generate_command_line(artificial_line, argcp, argvp, copy);
   return TRUE;
}

// Get "parameter" from file ===============================================
string get_param(FILE *fh, const char *param, int skip, const char *option,
   int _errno, string errmsg, int exit) _THROW {
   long actual_pos=ftell(fh);
   fseek(fh,0,SEEK_SET);
   
   char    line[201];
   string  retval;
   bool found=FALSE;
   int skipped=0;

   // Read lines
   while (fgets (line, 200,fh) != NULL && !found) {
      if (line[0]==0)    continue;
      if (line[0]=='#')  continue;
      if (line[0]==';')  continue;
      if (line[0]=='\n') continue;
      
      int i=0; char *line_wo_spaces=line;
      while (*line_wo_spaces==' ' || *line_wo_spaces=='\t') line_wo_spaces++;
      while (line_wo_spaces[i]!=0 && line_wo_spaces[i]!='=') i++;
      if (line_wo_spaces[i]=='=') {
         line_wo_spaces[i]=0;
         if (strcmp(line_wo_spaces,param)==0) {
            if (skipped==skip) {retval=line_wo_spaces+i+1; found=TRUE; break;}
            else skipped++;
         }
      }
   }
   fseek(fh,actual_pos,SEEK_SET);
   if (!found)
      if (option == NULL) {
         if (_errno==-1) {
            _errno=2104;
            errmsg=(string)"Argument "+param+" not found or invalid argument";
         }
         if (exit) EXIT_ERROR(_errno,errmsg);
         else      REPORT_ERROR(_errno,errmsg);
      } else return option;
   else return remove_spaces(retval);
}

// Check "parameter" from file =============================================
bool check_param(FILE *fh, const char *param) {
   long actual_pos=ftell(fh);
   fseek(fh,0,SEEK_SET);
   
   char    line[201];
   bool found=FALSE;
   string retval;

   // Read lines
   while (fgets (line, 200,fh) != NULL) {
      if (line[0]==0)    continue;
      if (line[0]=='#')  continue;
      if (line[0]=='\n') continue;
      
      int i=0; while (line[i]!=0 && line[i]!='=') i++;
      if (line[i]=='=') {
         line[i]=0;
         if (strcmp(line,param)==0) {retval=line+i+1; found=TRUE; break;}
      }
   }
   fseek(fh,actual_pos,SEEK_SET);
   return found && retval!="No" && retval!="NO" && retval!="no";
}

// Get vector param from file ==============================================
matrix1D<double> get_vector_param(FILE *fh, const char *param,
   int dim, int _errno,  string errmsg, int exit) _THROW {
   int    argcp;
   char **argvp=NULL;
   char  *copy=NULL;
   matrix1D<double> retval; 
   if (!generate_command_line(fh, param, argcp, argvp, copy))
      if (dim==-1) return retval;
      else if (exit) EXIT_ERROR(_errno,errmsg);
      else      REPORT_ERROR(_errno,errmsg);
   else {
      retval=get_vector_param(argcp, argvp, ((string)"-"+param).c_str(), dim,
         _errno, errmsg, exit);
      delete copy;
      return retval;
   }
}

// Instantiation ===========================================================
template <class T>
void instantiateArgsTemplate(T t) {
  vector<T> v1;
  string line;
  int i;
  read_float_list(line.c_str(), 1, v1);
  read_float_list(line, i, 1, v1);
}

void instantiateArgs() {
   float f1; 
   instantiateArgsTemplate(f1);
   double f2; 
   instantiateArgsTemplate(f2);
   int f3; 
   instantiateArgsTemplate(f3);   
}
