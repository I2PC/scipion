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
#include "../xmippFuncs.hh"
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <complex.h>
#include <fstream.h>
#include <typeinfo>

/* Numerical functions ----------------------------------------------------- */
// Solve second degree equation. ax^2+bx+c=0 -------------------------------
int solve_2nd_degree_eq(float a, float b, float c, float &x1, float &x2,
   float prec) {
   // Degenerate case?
   if (ABS(a)<prec)
      if (ABS(b)<prec) return -1;
      else             {x1=-c/b; return 1;}

   // Normal case
   float d=b*b-4*a*c;
   if   (d<0)     return 0;
   else {
      x1=(-b+sqrt(d))/(2*a);
      x2=(-b-sqrt(d))/(2*a);
      return 2;
   }
}

/* Gaussian value ---------------------------------------------------------- */
double gaussian1D(double x, double sigma, double mu) {
   x -= mu;
   return 1/sqrt(2*PI*sigma*sigma)*exp(-0.5*((x/sigma)*(x/sigma)));
}

double gaussian2D(double x, double y, double sigmaX, double sigmaY,
   double ang, double muX, double muY) {
   // Express x,y in the gaussian internal coordinates
   x -= muX;
   y -= muY;
   double xp= cos(ang)*x+sin(ang)*y;
   double yp=-sin(ang)*x+cos(ang)*y;
   
   // Now evaluate
   return 1/sqrt(2*PI*sigmaX*sigmaY)*exp(-0.5*((xp/sigmaX)*(xp/sigmaX)+
      (yp/sigmaY)*(yp/sigmaY)));
}

/* Print a boolean value --------------------------------------------------- */
void print(ostream &o, const bool b) {
   if (b) o << "TRUE"; else o << "FALSE";
}
/* Print a value in binary ------------------------------------------------- */
template <class T>
void printb(ostream &o, T value) {
    char buf [CHAR_BIT * sizeof(T) + 1];
    size_t i;

    for (i = 0; i < CHAR_BIT * sizeof(T); ++i)
    {
        buf[i] = '0' + (value & 1);
        value >>= 1;
    }
    buf[i] = 0;

    o << buf;
}

/* Random functions -------------------------------------------------------- */
int idum;
// Uniform distribution ....................................................
void  init_random_generator()     {idum=-1; ran1(&idum);}
void  randomize_random_generator(){static  unsigned int seed;
                                   int rand_return;
                                   
                                   srand(seed);
                                   rand_return=rand();

                                   time_t t; 
                                   time(&t); 
                                   rand_return=abs(rand_return);
                                   idum=( -(int)(t%10000)
                                          -(int)(rand_return%10000));
                                   ran1(&idum);
                                   seed=(unsigned int)rand_return;}
                                 
float rnd_unif()                  {return ran1(&idum);}
float rnd_unif(float a, float b)  {if (a==b) return a;
                                   else return a+(b-a)*ran1(&idum);}

// Gaussian distribution ...................................................
float rnd_gaus()                 {return gasdev(&idum);}
float rnd_gaus(float a, float b) {if (b==0) return a;
                                  else return b*gasdev(&idum)+a;}
float gaus_within_x0(float x0, float mean, float stddev) {
   float z0=(x0-mean)/stddev;
   return erf(ABS(z0)/sqrt(2.0));
}

float gaus_outside_x0(float x0, float mean, float stddev) {
   float z0=(x0-mean)/stddev;
   return erfc(ABS(z0)/sqrt(2.0));
}

float gaus_up_to_x0(float x0, float mean, float stddev) {
   if      (x0 >mean) return 1.0-gaus_outside_x0(x0,mean,stddev)/2;
   else if (x0==mean) return 0.5;
   else               return gaus_outside_x0(x0,mean,stddev)/2;
}

float gaus_from_x0(float x0, float mean, float stddev) {
   if      (x0 >mean) return gaus_outside_x0(x0,mean,stddev)/2;
   else if (x0==mean) return 0.5;
   else               return 1.0-gaus_outside_x0(x0,mean,stddev)/2;
}

float gaus_outside_probb(float p, float mean, float stddev) {
   // Make a Bolzano search for the right value
   float p1, p2, pm, x1, x2, xm;
   x1=mean;
   x2=mean+5*stddev;
   do {
      xm=(x1+x2)/2;
      p1=gaus_outside_x0(x1,mean,stddev);
      p2=gaus_outside_x0(x2,mean,stddev);
      pm=gaus_outside_x0(xm,mean,stddev);
      if (pm>p) x1=xm; else x2=xm;
   } while (ABS(pm-p)/p>0.005);
   return xm;
}

// See Numerical Recipes, Chap. 6.3
float student_within_t0(float t0, float degrees_of_freedom) {
   return 1-betai(degrees_of_freedom/2,0.5,
      degrees_of_freedom/(degrees_of_freedom+t0*t0));
}

float student_outside_t0(float t0, float degrees_of_freedom) {
   return 1-student_within_t0(t0,degrees_of_freedom);
}

float student_up_to_t0(float t0, float degrees_of_freedom) {
   if   (t0 >=0) return 1.0-student_outside_t0(t0,degrees_of_freedom)/2;
   else          return student_outside_t0(t0,degrees_of_freedom)/2;
}

float student_from_t0(float t0, float degrees_of_freedom) {
   return 1-student_up_to_t0(t0,degrees_of_freedom);
}

float student_outside_probb(float p, float degrees_of_freedom) {
   // Make a Bolzano search for the right value
   float p1, p2, pm, t1, t2, tm;
   t1=0;
   t2=100;
   do {
      tm=(t1+t2)/2;
      p1=student_outside_t0(t1,degrees_of_freedom);
      p2=student_outside_t0(t2,degrees_of_freedom);
      pm=student_outside_t0(tm,degrees_of_freedom);
      if (pm>p) t1=tm; else t2=tm;
   } while (ABS(pm-p)/p>0.005);
   return tm;
}

// Log uniform distribution ................................................
float rnd_log(float a, float b)  {if (a==b) return a;
                                  else return exp(rnd_unif(log(a),log(b)));}


/* Exception handling ------------------------------------------------------ */
void _Xmipp_error (const int nerr, const string &what) {
   cout << nerr << ": " << what << endl;
   exit(nerr);
}
#ifndef _NO_EXCEPTION
   // Object Constructor
   Xmipp_error::Xmipp_error(const int nerr, const string &what) {
      __errno=nerr;
      msg=what;
   }

   // Show message 
   ostream& operator << (ostream& o, Xmipp_error &XE) {
      o << XE.__errno << ":" << XE.msg << endl;
   return o;
   }
#endif


/* Handling with filenames ------------------------------------------------- */
int exists (const FileName &fn) {
   FILE *aux;
   if ((aux = fopen (fn.c_str(), "r")) == NULL) return 0;
   fclose (aux);
   return 1;
}

/* Create empty file ------------------------------------------------------- */
void create_empty_file(const FileName &fn, size_t size,
   size_t block_size) _THROW {
   unsigned char * buffer = (unsigned char*) calloc(sizeof(unsigned char),
      block_size);
   if (buffer==NULL) 
      REPORT_ERROR(1,"create_empty_file: No memory left");
   FILE * fd = fopen(fn.c_str(), "w");
   if (fd==NULL)
      REPORT_ERROR(1,(string)"create_empty_file: Cannot open file" + fn);
   for (size_t i=0; i< size/block_size; i++)
       fwrite( buffer, sizeof(unsigned char), block_size, fd);
   fwrite(buffer, sizeof(unsigned char), size%block_size, fd);
   fclose(fd);
}

// Constructor with root, number and extension .............................
void FileName::compose(const string &str, int no, const string &ext) {
   *this=(FileName) str;
   if (no!=-1) {
      char aux_str[5];
      sprintf(aux_str,"%05d",no);
      *this += aux_str;
   }

   if (ext!="") *this += "."+ext;
}

// Get the root name of a filename .........................................
FileName FileName::get_root() const {
   int skip_directories=find_last_of("/")+1;
   int point=find_first_of(".",skip_directories);
   if (point==-1) point=length();
   int root_end=find_last_not_of("0123456789",point-1);
   if (root_end+1 != point)
      if (point-root_end>5) root_end=point-5-1;
   return (FileName) substr(0,root_end+1);
}

// Get the base name of a filename .........................................
string FileName::get_baseName() const {
    string basename = "";
    string myname = *this;
    int myindex = 0;
    for (int p = myname.size()-1; p >=0; p--) {
       if (myname[p] == '/') {
         myindex = p+1;
	 break;
       }
    }
    for (int p = myindex; p <myname.size(); p++) {
       if (myname[p] != '.')
         basename += myname[p];
       else break;
    }
   return basename;
}

// Get number from file ....................................................
int FileName::get_number() const {
   int skip_directories=find_last_of("/")+1;
   int point=find_first_of(".",skip_directories);
   if (point==-1) point=length();
   int root_end=find_last_not_of("0123456789",point-1);
   if (root_end+1 != point) {
      if (point-root_end>5) root_end=point-5-1;
      string aux=substr(root_end+1,point-root_end+1);
     
      return atoi(aux.c_str());
   } else return -1;
}

// Get the extension of a filename .........................................
string FileName::get_extension() const {
   int skip_directories=find_last_of("/")+1;
   int first_point=find_first_of(".",skip_directories);
   if (first_point==-1) return "";
   else  return substr(first_point+1);
}

// Add at beginning ........................................................
FileName FileName::add_prefix(const string &prefix) const {
   FileName retval=*this;
   int skip_directories=find_last_of("/")+1;
   return retval.insert(skip_directories,prefix);
}

// Add at the end ..........................................................
FileName FileName::add_extension(const string &ext) const {
   if (ext=="") return *this;
   else {FileName retval=*this; retval=retval.append((string)"."+ext); return retval;}
}

// Remove last extension ...................................................
FileName FileName::without_extension() const {
   FileName retval=*this;
   return retval.substr(0,rfind("."));
}

// Remove root .............................................................
FileName FileName::without_root() const {return without(get_root());}

// Insert before extension .................................................
FileName FileName::insert_before_extension(const string &str) const {
   int point=-1;
   bool done=FALSE;
   do {
      point=find(".",point+1);
      if (point==-1) {point=length(); done=TRUE;}
      else if (point==length()-1) done=TRUE;
      else if ((*this)[point+1]=='.' || (*this)[point+1]=='/') done=FALSE;
   } while (!done);
   FileName retval=*this;
   return retval.insert(point,str);
}

// Remove an extension wherever it is ......................................
FileName FileName::remove_extension(const string &ext) const {
   int first=find((string)"."+ext);
   if (first==-1) return *this;
   else {FileName retval=*this; return retval.erase(first,1+ext.length());}
}

// Remove all extensions....................................................
FileName FileName::remove_all_extensions() const {
   int first=find("/");
   first=find(".",first+1);
   if (first==-1) return *this;
   else return substr(0,first);
}

// Substitute one extension by other .......................................
FileName FileName::substitute_extension(const string &ext1,
   const string &ext2) const {
   int first=find((string)"."+ext1);
   if (first==-1) return *this;
   else {
      FileName retval=*this;
      return retval.replace(first,1+ext1.length(),(string)"."+ext2);
   }
}

// Remove a substring ......................................................
FileName FileName::without(const string &str) const {
   if (str.length()==0) return *this;
   int pos=find(str);
   if (pos==-1) return *this;
   else {FileName retval=*this; return retval.erase(pos,str.length());}
}

// Remove until prefix .....................................................
FileName FileName::remove_until_prefix(const string &str) const {
   if (str.length()==0) return *this;
   int pos=find(str);
   if (pos==-1) return *this;
   else {FileName retval=*this; return retval.erase(0,pos+str.length());}
}

// Remove directories ......................................................
FileName FileName::remove_directories() const {
   FileName retval=*this;
   return retval.substr(0,rfind("/"));
}

/* Time managing ----------------------------------------------------------- */
#ifdef _NO_TIME
   void time_config() {}
   void annotate_time(TimeStamp *time) {}
   void print_elapsed_time(TimeStamp &time) {}
   float elapsed_time(TimeStamp &time) {}
   float time_to_go(TimeStamp &time, float fraction_done) {}
   void TimeMessage(string message) {}
   void progress_bar(long rlen) {}
#else
// A global ................................................................
int XmippTICKS;

// Time configuration ......................................................
// The clock frequency for each machine must be known
void time_config() {XmippTICKS = sysconf(_SC_CLK_TCK);}

// Annotate actual time ....................................................
void annotate_time(TimeStamp *time) {times(time);}

// Show elapsed time since last annotation .................................
void print_elapsed_time(TimeStamp &time, int _IN_SECS) {
   TimeStamp now; times(&now);
   float userTime=now.tms_utime - time.tms_utime;
   float sysTime=now.tms_stime - time.tms_stime;
   if (_IN_SECS) {userTime /= XmippTICKS; sysTime /=XmippTICKS;}
   cout << "Elapsed time: User(" << userTime << ") System(" << sysTime
        << ")\n";
}

// Calculate elapsed time since last annotation .............................
float elapsed_time(TimeStamp &time, int _IN_SECS) {
   TimeStamp now; times(&now);
   float userTime=now.tms_utime - time.tms_utime;
   float sysTime=now.tms_stime - time.tms_stime;
   if (_IN_SECS) {userTime /= XmippTICKS; sysTime /=XmippTICKS;}
   return userTime+sysTime;
}

// Compute the predicted time left .........................................
float time_to_go(TimeStamp &time, float fraction_done) {
   TimeStamp now; times(&now);
   float totalTime=(now.tms_utime - time.tms_utime +
      now.tms_stime - time.tms_stime)/XmippTICKS;
   return totalTime*(1-fraction_done)/fraction_done;
}

// Show a message with the time it is produced .............................
void TimeMessage(string message) {
   struct tm *T;
   time_t     seconds;

   if( time(&seconds) <0) seconds=0;
   T=localtime(&seconds);

   printf("%2d:%2d:%2d (day=%2d) =>%s ",T->tm_hour,
      T->tm_min,T->tm_sec,T->tm_mday,message.c_str());
}

// Show a bar with the progress in time ....................................
// When the input is negative then we are setting the progress bar, this
// will be the total of elements to process. Afterwards the call to this
// routine must be in ascending order, ie, 0, 1, 2, ... No. elements
void progress_bar(long rlen) {
   static time_t startt, prevt;
   time_t currt;
   static long totlen;
   long t1, t2;
   int min, i, hour;
   float h1, h2, m1, m2;

   if (rlen==0) return;
   currt=time(NULL);

   if (rlen<0) {
      totlen=-rlen;
      prevt=startt=currt;
      fprintf(stderr,"0000/???? sec. ");
      for (i=0; i<10; i++) fprintf(stderr,"------");
      fflush(stderr);
   } else if (totlen>0) {
      t1=currt-startt;                   // Elapsed time
      t2=(long) (t1*(float)totlen/rlen); // Total time

      hour = 0; min = 0;
      if (t2>60) {
         m1 = (float)t1/60.0;
         m2 = (float)t2/60.0;
         min=1;
      	 if (m2>60) {
          h1 = (float)m1/60.0;
          h2 = (float)m2/60.0;
          hour=1;
	  min = 0;
         } else hour=0;	 
      } else min = 0; 	  
      
      if (hour)
         fprintf(stderr,"\r%3.2f/%3.2f %s ", h1, h2, "hours");
      else if (min)
         fprintf(stderr,"\r%3.2f/%3.2f %s ", m1, m2, "min");
      else 	 
         fprintf(stderr,"\r%4u/%4u %4s ", (int)t1, (int)t2, "sec.");

      i=(int) (60*(1-(float)(totlen-rlen)/totlen));
      while (i--) fprintf(stderr,".");
      if (rlen==totlen) {fprintf(stderr,"\n"); totlen=0;}
      fflush(stderr);
      prevt=currt;
   }
}


// Initialize progress bar.

void xmippTextualListener::OnInitOperation(unsigned long _est_it) {
	progress_bar(-(_est_it));
}
	
// Show a bar with the progress in time ....................................
// When the input is negative then we are setting the progress bar, this
// will be the total of elements to process. Afterwards the call to this
// routine must be in ascending order, ie, 0, 1, 2, ... No. elements
// Almost identical to previous progress bar function, in fact, we call
// those functions inside.

void xmippTextualListener::OnProgress(unsigned long _it) {
	progress_bar(_it);
}	

// Shows a message indicating the operation in progress.
void xmippTextualListener::OnReportOperation(const string& _rsOp) {
	fprintf(stderr, _rsOp.c_str());//	cout << _rsOp;
}


#endif

/* Little/big endian ------------------------------------------------------- */
// Read in reverse/normal order --------------------------------------------
size_t FREAD(void *dest, size_t size, size_t nitems, FILE * &fp, bool reverse) {
   size_t retval;
   if (!reverse)
      retval=fread(dest,size,nitems,fp);
   else {
      char *ptr=(char *)dest;
      bool end=FALSE;
      retval=0;
      for (int n=0; n<nitems; n++) {
         for (int i=size-1; i>=0; i--) {
             if (fread(ptr+i,1,1,fp)!=1) {end=TRUE; break;}
         }
         if (end) break; else retval++;
         ptr +=size;
      }
   }
   return retval;
}

// Read in reverse/normal order --------------------------------------------
size_t FWRITE(const void *src, size_t size, size_t nitems, FILE * &fp,
   bool reverse) {
   size_t retval;
   if (!reverse)
      retval=fwrite(src,size,nitems,fp);
   else {
      char *ptr=(char *)src;
      bool end=FALSE;
      retval=0;
      for (int n=0; n<nitems; n++) {
         for (int i=size-1; i>=0; i--) {
            if (fwrite(ptr+i,1,1,fp)!=1) {end=TRUE; break;}
         }
         if (end) break; else retval++;
         ptr +=size;
      }
   }
   return retval;
}

// Managing memory ---------------------------------------------------------
template <class T> void ask_Tvector(T* &v, int nl, int nh) _THROW {
   if (nh-nl+1>1) {
      v=(T *)malloc((unsigned) (nh-nl+1)*sizeof(T));
      if (!v) REPORT_ERROR(1,"allocation failure in vector()");
      v-=nl;
   } else v=NULL;
}

template <class T> void free_Tvector(T* &v, int nl, int nh) {
   if (v!=NULL) {
      free((char*) (v+nl));
      v=NULL;
   }
}

template <class T> void ask_Tmatrix(T ** &m, int nrl, int nrh,
   int ncl, int nch) _THROW {
   if (nrh-nrl+1>1 && nch-ncl+1>1) {
      m=(T **) malloc((unsigned) (nrh-nrl+1)*sizeof(T*));
      if (!m) REPORT_ERROR(1,"allocation failure 1 in matrix()");
      m -= nrl;

      for(int i=nrl;i<=nrh;i++) {
              m[i]=(T *) malloc((unsigned) (nch-ncl+1)*sizeof(T));
              if (!m[i]) REPORT_ERROR(1,"allocation failure 2 in matrix()");
              m[i] -= ncl;
      }
   } else m=NULL;
}

template <class T> void free_Tmatrix(T ** &m, int nrl, int nrh,
   int ncl, int nch) {
   if (m!=NULL) {
      for(int i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
      free((char*) (m+nrl));
      m=NULL;
   }
}

template <class T> void ask_Tvolume(T *** &m, int nsl, int nsh, int nrl,
   int nrh, int ncl, int nch) _THROW {
   if (nsh-nsl+1>1 && nrh-nrl+1>1 && nch-ncl+1>1) {
      m=(T ***) malloc((unsigned) (nsh-nsl+1)*sizeof(T**));
      if (!m) REPORT_ERROR(1,"allocation failure 1 in matrix()");
      m -= nsl;

      for (int k=nsl;k<=nsh;k++) {
          m[k]=(T **) malloc((unsigned) (nrh-nrl+1)*sizeof(T*));
          if (!m[k]) REPORT_ERROR(1,"allocation failure 2 in matrix()");
          m[k] -= nrl;

          for (int i=nrl;i<=nrh;i++) {
              m[k][i]=(T *) malloc((unsigned) (nch-ncl+1)*sizeof(T));
              if (!m[k][i]) REPORT_ERROR(1,"allocation failure 2 in matrix()");
              m[k][i] -= ncl;
          }
      }
   } else m=NULL;
}

template <class T> void free_Tvolume(T *** &m, int nsl, int nsh,
   int nrl, int nrh, int ncl, int nch) {
   if (m!=NULL) {
      for(int k=nsh;k>=nsl;k--) {
	 for(int i=nrh;i>=nrl;i--) free((char*) (m[k][i]+ncl));
         free((char*) (m[k]+nrl));
      }
      free((char*) (m+nsl));
      m=NULL;
   }
}

/* Marsaglia, fast random number generator --------------------------------- */
template <class T>
  void Marsaglia<T>::Init(const FileName &fn_in, int No_Numbers) _THROW
  {
   int Type_size;            //sizeof(type)
   
   pointer_in_memory=0;
   Number_of_Numbers=No_Numbers; // initializa class variable
   Type_size=sizeof(T);

   ifstream in(fn_in.c_str());
   in.seekg(0, ios::end);         // End of file
   streampos sp = in.tellg();     // Size of file
   if(sp < Number_of_Numbers*Type_size)
       {
       REPORT_ERROR(1,(string)"Marsaglia::Init: File "+fn_in+"is too small");
       }
   else{
   //get a random number to set the file pointer at a random position
   randomize_random_generator();  // seed the random generator

   random_vector =  new char[(Number_of_Numbers*Type_size)];
   T_random_vector = (T *) random_vector;
   in.seekg((streampos) FLOOR ( rnd_unif(0.f,(float) (sp - 
                                Number_of_Numbers*Type_size) ) ),ios::beg);
   in.read(random_vector, (Number_of_Numbers*Type_size));

    in.close();
   
   }
   if( typeid(float) == typeid(T))
      {
      Verify_float(); 
      }
      
}
//--------------------------------------------------------------
// produce valid random floats instead of random bits
template <class T>
  void Marsaglia<T>::Verify_float()
  {
  unsigned int * int_random_vector;
  long long MaxInteger;
  if(sizeof(float)!= sizeof(int))
       {
       cout << "\nMarsaglia: I do not know how to make the float correction\n";
       exit(0);
       } 
  MaxInteger = (long long)pow(2.0,sizeof(unsigned int)*8.0);
  int_random_vector = (unsigned int *) random_vector;
     for(int hh=0; hh< Number_of_Numbers; hh++)
        T_random_vector[hh]= (T)((double)int_random_vector[hh]/
                                                         (double)MaxInteger);
  }
//--------------------------------------------------------------
// calculate the logarithm. if zero make it small
template <class T>
  void Marsaglia<T>::Marsaglia_log(void)
  {
   if(typeid(float)!=typeid(T) && typeid(double)!=typeid(T))
       {
       cout << "\nMarsaglia: I do not know how to calculate integer logs\n";
       exit(0);
       } 
   
  for(int hh=0; hh< Number_of_Numbers; hh++)
     if(T_random_vector[hh]==0.) T_random_vector[hh]= -1e+20f;
     else T_random_vector[hh]=log(T_random_vector[hh]);
  }
//--------------------------------------------------------------
// Multiply by a constant
template <class T>
  void Marsaglia<T>::mul(T mul_cte)
  {
   
  for(int hh=0; hh< Number_of_Numbers; hh++)
     T_random_vector[hh] *= mul_cte;
  }
//--------------------------------------------------------------
// Multiply by a constant
template <class T>
  void Marsaglia<T>::operator &=(T mod_cte)
  {
   
  for(int hh=0; hh< Number_of_Numbers; hh++)
     T_random_vector[hh] &= mod_cte;
  }
//--------------------------------------------------------------
// Add a constant
template <class T>
  void Marsaglia<T>::add(T add_cte)
  {
   
  for(int hh=0; hh< Number_of_Numbers; hh++)
     T_random_vector[hh] += add_cte;
  }
// Fix the maximum value only valid for integers
template <class T>
  void Marsaglia<T>::M_max(const FileName &fn_in, T m_max)
  {
   int Type_size;                 //sizeof(type)
   Type_size=sizeof(T);

   ifstream in(fn_in.c_str());
   in.seekg(0, ios::end);         // End of file
   streampos sp = in.tellg();     // Size of file
   T power_of_2 =(T)NEXT_POWER_OF_2(m_max);
   if (power_of_2==m_max)
      power_of_2=(T)NEXT_POWER_OF_2(m_max+1);
   T mask=power_of_2-1;  
   T aux_number;
   m_max;
   //get a random number to set the file pointer at a random position
   in.seekg((streampos) FLOOR ( rnd_unif(0.f,(float) (sp - 
                                Number_of_Numbers*Type_size) ) ),ios::beg);
    for(int ii=0; ii<Number_of_Numbers;)
       {  
          
          aux_number  = T_random_vector[ii];
          aux_number &= mask;
          if(aux_number > m_max ||/* 
             aux_number < -m_max  ||*/
             (T_random_vector[ii] <= 0) && (aux_number==0) )
             {
             if(in.eof())
                 in.seekg((streampos) FLOOR ( rnd_unif(0.f,(float) (sp - 
                                Number_of_Numbers*Type_size) ) ),ios::beg);
        
             in.read((char*)&(T_random_vector[ii]),Type_size);
             }
          else
             {
             T_random_vector[ii] = aux_number*(T)SGN(T_random_vector[ii]); 
             ii++;
             }

       }//end for   
    
   in.close();
         
  }
  
/* Instantiate ------------------------------------------------------------- */
void instantiante_xmippFuncs() {
   char           *h, **hh, ***hhh;
   short int      *s, **ss, ***sss;
   int            *i, **ii, ***iii;
   float          *f, **ff, ***fff;
   double         *d, **dd, ***ddd;
   double_complex *c, **cc, ***ccc;

   ask_Tvector(h,  1,1);         free_Tvector(h,  1,1);
   ask_Tmatrix(hh, 1,1,1,1);     free_Tmatrix(hh, 1,1,1,1);
   ask_Tvolume(hhh,1,1,1,1,1,1); free_Tvolume(hhh,1,1,1,1,1,1);
   ask_Tvector(s,  1,1);         free_Tvector(s,  1,1);
   ask_Tmatrix(ss, 1,1,1,1);     free_Tmatrix(ss, 1,1,1,1);
   ask_Tvolume(sss,1,1,1,1,1,1); free_Tvolume(sss,1,1,1,1,1,1);
   ask_Tvector(i,  1,1);         free_Tvector(i,  1,1);
   ask_Tmatrix(ii, 1,1,1,1);     free_Tmatrix(ii, 1,1,1,1);
   ask_Tvolume(iii,1,1,1,1,1,1); free_Tvolume(iii,1,1,1,1,1,1);
   ask_Tvector(f,  1,1);         free_Tvector(f,  1,1);
   ask_Tmatrix(ff, 1,1,1,1);     free_Tmatrix(ff, 1,1,1,1);
   ask_Tvolume(fff,1,1,1,1,1,1); free_Tvolume(fff,1,1,1,1,1,1);
   ask_Tvector(d,  1,1);         free_Tvector(d,  1,1);
   ask_Tmatrix(dd, 1,1,1,1);     free_Tmatrix(dd, 1,1,1,1);
   ask_Tvolume(ddd,1,1,1,1,1,1); free_Tvolume(ddd,1,1,1,1,1,1);
   ask_Tvector(c,  1,1);         free_Tvector(c,  1,1);
   ask_Tmatrix(cc, 1,1,1,1);     free_Tmatrix(cc, 1,1,1,1);
   ask_Tvolume(ccc,1,1,1,1,1,1); free_Tvolume(ccc,1,1,1,1,1,1);
   
   Marsaglia<float> Random_pool_f("m",1);
   Random_pool_f.Get_One_Number();
   Random_pool_f.Marsaglia_log();//omly floats
   Random_pool_f.mul((float)1.);
   Random_pool_f.add((float)1.);

   Marsaglia<int> Random_pool_i("m",1);
   Random_pool_i.Get_One_Number();
   Random_pool_i.mul((int)1);
   Random_pool_i.add((int)1);
   Random_pool_i&=(int)1;//only integers
   Random_pool_i.M_max("m",(int) 1);//only integers

   Marsaglia<unsigned int> Random_pool_ui("m",1);
   Random_pool_ui.Get_One_Number();
   Random_pool_ui.mul((unsigned int)1);
   Random_pool_ui.add((unsigned int)1);
   Random_pool_ui&=(unsigned int)1;//only integers
   Random_pool_ui.M_max("m",(unsigned int) 1);//only integers

   Marsaglia<short int> Random_pool_s("m",1);
   Random_pool_s.Get_One_Number();
   Random_pool_s.mul((short int)1);
   Random_pool_s.add((short int)1);
   Random_pool_s&=(short int) 1;//only integers
   Random_pool_s.M_max("m",(short int)1);//only integers

   Marsaglia<unsigned short int> Random_pool_us("m",1);
   Random_pool_us.Get_One_Number();
   Random_pool_us.mul((unsigned short int)1);
   Random_pool_us.add((unsigned short int)1);
   Random_pool_us&=(unsigned short int)1;//only integers
   Random_pool_us.M_max("m",(unsigned short int)1);//only integers

   printb(cout, (char)0);
   printb(cout, (unsigned char)0);
   printb(cout, (int)0);
   printb(cout, (unsigned int)0);
   printb(cout, (long)0);
   printb(cout, (unsigned long)0);
   printb(cout, (long long)0);
   printb(cout, (unsigned long long)0);
}
