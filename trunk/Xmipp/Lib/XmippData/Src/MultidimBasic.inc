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
 * MERCHANTABILITY or FITNESS FOR Aouble_ PARTICULAR PURPOSE.  See the       
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

/* ------------------------------------------------------------------------- */
/* MULTIDIM BASIC                                                            */
/* ------------------------------------------------------------------------- */

#include <complex>
#define mi MULTIDIM_ELEM(*this,i)
#define msize MULTIDIM_SIZE(*this)
#ifndef double_complex
   typedef complex<double> double_complex;
#endif

/* Print stats ------------------------------------------------------------- */
template <class T>
   void maT::print_stats(ostream &out) const {
      T min_val, max_val;
      double avg_val, dev_val;
      
      compute_stats(avg_val, dev_val, min_val, max_val);

      out.setf(ios::showpoint); int old_prec=out.precision(7); 
      out << " min= "; out.width(9); out << min_val;
      out << " max= "; out.width(9); out << max_val;
      out << " avg= "; out.width(9); out << avg_val;
      out << " dev= "; out.width(9); out << dev_val;
      out.precision(old_prec);
   }

/* Compute max ------------------------------------------------------------- */
template <class T>
   T maT::compute_max() const {
      if (__dim<=0) return (T)0;
      T max_val=MULTIDIM_ELEM(*this,0);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         if (mi>max_val) max_val=mi;
      return max_val;
   }

/* Compute min ------------------------------------------------------------- */
template <class T>
   T maT::compute_min() const {
      if (__dim<=0) return (T)0;
      T min_val=MULTIDIM_ELEM(*this,0);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         if (mi<min_val) min_val=mi;
      return min_val;
   }

/* Compute min-max --------------------------------------------------------- */
template <class T>
   void maT::compute_double_minmax(double &min, double &max) const {
      min=max=0.0f;
      if (__dim<=0) return;
      min=max=(double) MULTIDIM_ELEM(*this,0);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) {
         if (mi<min) min=(double)mi;
         if (mi>max) max=(double)mi;
      }
   }

template <class T>
   void maT::compute_double_minmax(double &min, double &max,
      const matrix1D<int> &corner1, const matrix1D<int> &corner2) const {
      matrix1D<double> dcorner1, dcorner2;
      type_cast(corner1,dcorner1);
      type_cast(corner2,dcorner2);
      compute_double_minmax(min,max,dcorner1,dcorner2);
   }

template <class T>
   void maT::compute_stats(double &avg, double &stddev, T &min_val, T &max_val,
      const matrix1D<int> &corner1, const matrix1D<int> &corner2) const {
      matrix1D<double> dcorner1, dcorner2;
      type_cast(corner1,dcorner1);
      type_cast(corner2,dcorner2);
      compute_stats(avg,stddev,min_val,max_val,dcorner1,dcorner2);
   }

/* Compute avg ------------------------------------------------------------- */
template <class T>
   double maT::compute_avg() const {
      if (__dim<=0) return 0;
      double sum=0;
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         sum += (double) mi;
      return sum/msize;
   }

/* Compute stddev ---------------------------------------------------------- */
template <class T>
   double maT::compute_stddev() const {
      if (msize<=0) return 0;
      double avg=0, stddev=0;
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) {
         avg    += (double) mi;
         stddev += (double) mi * (double) mi;
      }
      if (msize>1) {
         avg/=msize;
         stddev = stddev/msize - avg*avg;
         stddev*= msize/(msize-1);
         stddev = sqrt((double)(ABS(stddev))); // Foreseeing numerical
                                               // instabilities
      } else stddev=0;
      return stddev;
   }

/* Compute stats ---------------------------------------------------------- */
template <class T>
   void maT::compute_stats(double &avg, double &stddev,
      T &min, T &max) const {
      avg=0; stddev=0; min=(T)0; max=(T)0;
      if (msize<=0) return;

      min=max=MULTIDIM_ELEM(*this,0);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) {
         avg    += (double) mi;
         stddev += (double) mi * (double) mi ;
         if (mi>max) max=mi;
         if (mi<min) min=mi;
      }
      avg /=msize;
      if (msize>1) {
         stddev = stddev/msize - avg*avg;
         stddev*= msize/(msize-1);
         stddev = sqrt((double)(ABS(stddev))); // Foreseeing numerical
                                               // instabilities
      } else stddev=0;
   }

/* Range adjust ------------------------------------------------------------ */
// This function takes a vector with range between min0 ... max0 and
// linearly transforms it to minF ... maxF.
// If the input vector is a constant then it is adjusted to minF
template <class T>
   void maT::range_adjust(T minF, T maxF) {
   if (msize==0) return;
   T min0=compute_min();
   T max0=compute_max();
   
   // If max0==min0, it means that the vector is a constant one, so the
   // only possible transformation is to a fixed minF
   double slope;
   if (max0!=min0) slope=(double)(maxF-minF)/(double)(max0-min0);
   else            slope=0;
   
   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
      mi=minF+(T) (slope * (double) (mi-min0));
}

/* Statistics Adjust ------------------------------------------------------- */
// This function takes a vector with an average avg0 and standard deviation
// dev0, and transforms it linearly into a vector with avgF and devF
// If input vector has got dev0==0, ie, is a constant vector then only the
// average is adjusted, and the deviation remains being 0.
template <class T>
   void maT::statistics_adjust(double avgF, double stddevF) {
   double avg0,stddev0;
   double a,b;

   if (msize==0) return;
   T min, max;
   compute_stats(avg0, stddev0, min, max);
   
   if (stddev0!=0) a=(double)stddevF/(double)stddev0;
   else         a=0;
   b=(double)avgF - a*(double)avg0;
   
   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
      mi = (T) (a*(double)mi+b);
}

/* Init random ------------------------------------------------------------- */
template <class T>
   void maT::init_random(double op1, double op2, const string &mode) _THROW {
      if (mode=="uniform")
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
            mi=(T) rnd_unif(op1,op2);
      else if (mode=="gaussian")
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
            mi=(T) rnd_gaus(op1,op2);
      else
         REPORT_ERROR(1005,(string)"Init_random: Mode not supported ("+mode+")");
   }

/* Add Noise --------------------------------------------------------------- */
// Supported random distributions: 
//    "uniform", between op1 and op2
//    "gaussian", with avg=op1 and var=op2
// It is not an error that the vector is empty
template <class T>
   void maT::add_noise(double op1, double op2, const string &mode) const _THROW {
      if (mode=="uniform")
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
            mi +=(T) rnd_unif(op1,op2);
      else if (mode=="gaussian")
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
            mi +=(T) rnd_gaus(op1,op2);
      else
         REPORT_ERROR(1005,(string)"Add_noise: Mode not supported ("+mode+")");
   }

/* Threshold --------------------------------------------------------------- */
// Substitute component values by other according to the type of threshold
// to apply
// abs_above, abs_below, above, below, range
template <class T>
   void maT::threshold(const string &type, T a, T b) {
      int mode;

      if      (type == "abs_above") mode=1;
      else if (type == "abs_below") mode=2;
      else if (type == "above")     mode=3;
      else if (type == "below")     mode=4;
      else if (type == "range")     mode=5;
      else
         REPORT_ERROR(1005, (string)"Threshold: mode not supported ("+type+")");

      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         switch (mode) {
            case 1: if (ABS(mi)>a) mi=SGN(mi)*b; break;
            case 2: if (ABS(mi)<a) mi=SGN(mi)*b; break;
            case 3: if (mi>a)      mi=b; break;
            case 4: if (mi<a)      mi=b; break;
            case 5: if      (mi<a) mi=a;
                    else if (mi>b) mi=b; break;
         }
   }

/* Count with threshold ---------------------------------------------------- */
template <class T>
   long maT::count_threshold(const string &type, T a, T b) {
      int mode;

      if      (type == "abs_above") mode=1;
      else if (type == "abs_below") mode=2;
      else if (type == "above")     mode=3;
      else if (type == "below")     mode=4;
      else if (type == "range")     mode=5;
      else
         REPORT_ERROR(1005, (string)"Threshold: mode not supported ("+type+")");

      long retval=0;
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         switch (mode) {
            case 1: if (ABS(mi)>a) retval++; break;
            case 2: if (ABS(mi)<a) retval++; break;
            case 3: if (mi>a)      retval++; break;
            case 4: if (mi<a)      retval++; break;
            case 5: if (mi>=a && mi<=b) retval++; break;
         }
      return retval;
   }

/* Equality for normal data types ------------------------------------------ */
template <class T>
   bool operator == (const maT &op1, const maT &op2) {return op1.equal(op2);}

/* Equality for complex numbers -------------------------------------------- */
#ifdef ENABLE_COMPLEX_MATRICES
bool operator == (const ma< complex<double> > &op1, const ma< complex<double> > &op2) {
   if (!op1.same_shape(op2)) return FALSE;
   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(op1)
      if (ABS(MULTIDIM_ELEM(op1,i).real()-MULTIDIM_ELEM(op2,i).real())
             >XMIPP_EQUAL_ACCURACY
        ||ABS(MULTIDIM_ELEM(op1,i).imag()-MULTIDIM_ELEM(op2,i).imag())
             >XMIPP_EQUAL_ACCURACY)
          return FALSE;
   return TRUE;
}
#endif

/* Core array by scalar ---------------------------------------------------- */
template <class T>
void core_array_by_scalar(const maT &op1, const T &op2,
   maT &result, char operation) _THROW {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) + op2; break;
            case '-':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) - op2; break;
            case '*':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) * op2; break;
            case '/':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) / op2; break;
            case '^':
               MULTIDIM_ELEM(result,i)=
                  (T) pow((double)MULTIDIM_ELEM(op1,i),(double)op2); break;
         }
   }

template <>
void core_array_by_scalar< complex<double> >(const maTC &op1,
   const complex<double> &op2, maTC &result, char operation) _THROW {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) + op2; break;
            case '-':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) - op2; break;
            case '*':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) * op2; break;
            case '/':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) / op2; break;
            case '^':
               MULTIDIM_ELEM(result,i)=
                  pow(MULTIDIM_ELEM(op1,i),op2); break;
         }
   }

/* Scalar by array --------------------------------------------------------- */
template <class T>
void core_scalar_by_array(const T &op1, const maT &op2,
   maT &result, char operation) _THROW {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=op1 + MULTIDIM_ELEM(op2,i); break;
            case '-':
               MULTIDIM_ELEM(result,i)=op1 - MULTIDIM_ELEM(op2,i); break;
            case '*':
               MULTIDIM_ELEM(result,i)=op1 * MULTIDIM_ELEM(op2,i); break;
            case '/':
               MULTIDIM_ELEM(result,i)=op1 / MULTIDIM_ELEM(op2,i); break;
            case '^':
               MULTIDIM_ELEM(result,i)=(T)
                  pow((double)op1,(double)MULTIDIM_ELEM(op2,i)); break;
         }
   }

template <>
void core_scalar_by_array< complex<double> >(const complex<double> &op1,
   const maTC &op2, maTC &result, char operation) _THROW {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=op1 + MULTIDIM_ELEM(op2,i); break;
            case '-':
               MULTIDIM_ELEM(result,i)=op1 - MULTIDIM_ELEM(op2,i); break;
            case '*':
               MULTIDIM_ELEM(result,i)=op1 * MULTIDIM_ELEM(op2,i); break;
            case '/':
               MULTIDIM_ELEM(result,i)=op1 / MULTIDIM_ELEM(op2,i); break;
            case '^':
               MULTIDIM_ELEM(result,i)=pow(op1,MULTIDIM_ELEM(op2,i)); break;
         }
   }

/* Array by array ---------------------------------------------------------- */
template <class T>
void core_array_by_array(const maT &op1, const maT &op2,
   maT &result, char operation) _THROW {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) +
                                       MULTIDIM_ELEM(op2,i); break;
            case '-':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) -
                                       MULTIDIM_ELEM(op2,i); break;
            case '*':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) *
                                       MULTIDIM_ELEM(op2,i); break;
            case '/':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) /
                                       MULTIDIM_ELEM(op2,i); break;
            case '^':
               MULTIDIM_ELEM(result,i)=
                  (T) pow((double)MULTIDIM_ELEM(op1,i),
                          (double)MULTIDIM_ELEM(op2,i)); break;
         }
   }

template <>
void core_array_by_array< complex<double> >(const maTC &op1, const maTC &op2,
   maTC &result, char operation) _THROW {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) +
                                       MULTIDIM_ELEM(op2,i); break;
            case '-':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) -
                                       MULTIDIM_ELEM(op2,i); break;
            case '*':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) *
                                       MULTIDIM_ELEM(op2,i); break;
            case '/':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) /
                                       MULTIDIM_ELEM(op2,i); break;
            case '^':
               MULTIDIM_ELEM(result,i)=pow(MULTIDIM_ELEM(op1,i),
                                       MULTIDIM_ELEM(op2,i)); break;
         }
   }


