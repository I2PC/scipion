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
/* ************************************************************************* */
/* CLASS DEFINITION AND PROTOTYPES                                           */
/* ************************************************************************* */
   // This file contains several definitions for statistics and arithmetic
   // operations. To use it you have to redirect the internal type maT
   // (multidimensional array<T>) to the type you need. For example, for
   // a vector library would be
   //
   //#define maT vector<T>
   //
   // These definitions are outside because in this way we can reuse the
   // module for other libraries (matrices1D, matrices2D, ...). This way
   // a further advantage is that all functions are named the same in
   // all libraries

/* Structure =============================================================== */
public:
   T          *__m;              // Array
   int         __dim;            // In matrix1D dim=xdim;
                                 // In matrix2D dim=xdim*ydim;
                                 // In matrix3D dim=xdim*ydim*zdim;
   int         __spcdim;         // In matrix1D spcdim=1;
                                 // In matrix2D spcdim=2;
                                 // In matrix3D spcdim=3;

/* Methods ================================================================= */
public:
/**@name Statistics functions */
//@{
   /** Print statistics in current line.
       No end of line character is written after this print out.
       \\Ex: a.compute_stats();
       cout << "Statistics of variable a "; a.print_stats();
       cout << endl;*/
   void print_stats(ostream &out=cout) const;

   /** Maximum of the values in the array.
       The returned value is of the same type as the type of the array.
       \\Ex: double maxval=double_array.compute_max();*/
   T compute_max() const;

   /** Minimum of the values in the array.
       The returned value is of the same type as the type of the array.
       \\Ex: double minval=double_array.compute_min();*/
   T compute_min() const;

   /** Minimum and maximum of the values in the array.*/
   void compute_double_minmax(double &min, double &max) const;

   /** Minimum and maximum within region. The region is specified
       by two corners. */
   void compute_double_minmax(double &min, double &max,
      const matrix1D<double> &corner1, const matrix1D<double> &corner2) const;

   /** Minimum and maximum within region. The region is specified
       by two corners. */
   void compute_double_minmax(double &min, double &max,
      const matrix1D<int> &corner1, const matrix1D<int> &corner2) const;

   /** Average of the values in the array.
       The returned value is always double, independently of the type
       of the array.
       \\Ex: double avgval=int_array.compute_avg();*/
   double compute_avg() const;

   /** Standard deviation of the values in the array.
       Be careful that the standard deviation and NOT the variance is
       returned.
       The returned value is always double, independently of the type
       of the array.
       \\Ex: double devval=int_array.compute_stddev();*/
   double compute_stddev() const;

   /** Compute statistics.
       The average, standard deviation, minimum and maximum value are
       returned. */
   void compute_stats(double &avg, double &stddev, T &min, T &max) const;
   
   /** Compute statistics within region. The region is specified
       by two corners. */
   void compute_stats(double &avg, double &stddev, T &min, T &max,
      const matrix1D<int> &corner1, const matrix1D<int> &corner2) const;

   /** Compute statistics within region. The region is specified
       by two corners. */
   void compute_stats(double &avg, double &stddev, T &min, T &max,
      const matrix1D<double> &corner1, const matrix1D<double> &corner2) const;

   /** Adjust the range of the array to a given one.
       A linear operation is performed on the values of the array such
       that after it, the values of the array are comprissed between
       the two values set. The actual array is modified itself
       \\Ex: v.range_adjust(0,1);
       \\The array is now ranging from 0 to 1 */
   // This function must be explictly implemented outside
   void range_adjust(T minF, T maxF);
   
   /** Adjust the average and stddev of the array to given values.
       A linear operation is performed on the values of the array such
       that after it, the average and standard deviation of the array
       are the two values set. The actual array is modified itself
       \\Ex: v.statistics_adjust(0,1);
       \\The array has got now 0 mean and stddev=1 */
   // This function must be explictly implemented outside
   void statistics_adjust(double avgF, double stddevF);
//@}

/**@name Arithmethic operations */
//@{
   #define mi MULTIDIM_ELEM(*this,i)
   /* Array "by" array ----------------------------------------------------- */
   /**@name Array "by" array operations
      These are operations that are performed between 2 arrays of the
      SAME type (two integer vectors, two double matrices, ...). If they
      are not of the same type you can convert one of the arrays to the
      desired type using the function type_cast. The result must have been
      defined to be of the same type as the operands.

      In this kind of operations each element of array 1 is operated with its
      homologous in array 2, it is very important that both have got the
      same size and starting origins. The result has also got the same
      shape as the two operated arrays and its former content is lost.

      @see type_cast
      */
   //@{
   /** Core array by array operation.
       It assumes that the result is already resized. */
   friend void core_array_by_array<>(const maT &op1, const maT &op2,
      maT &result, char operation);

   /// v3=v1+v2
   maT operator  + (const maT &op1) const
      {maT temp; array_by_array(*this,op1,temp,'+'); return temp;}
   /// v3=v1-v2.
   maT operator  - (const maT &op1) const
      {maT temp; array_by_array(*this,op1,temp,'-'); return temp;}
   /// v3=v1*v2.
   maT operator  * (const maT &op1) const
      {maT temp; array_by_array(*this,op1,temp,'*'); return temp;}
   /// v3=v1/v2.
   maT operator  / (const maT &op1) const
      {maT temp; array_by_array(*this,op1,temp,'/'); return temp;}
   /// v3=v1^v2.
   maT operator  ^ (const maT &op1) const
      {maT temp; array_by_array(*this,op1,temp,'^'); return temp;}

   /// v3+=v2
   void operator  += (const maT &op1)
      {array_by_array(*this,op1,*this,'+');}
   /// v3-=v2.
   void operator  -= (const maT &op1)
      {array_by_array(*this,op1,*this,'-');}
   /// v3*=v2.
   void operator  *= (const maT &op1)
      {array_by_array(*this,op1,*this,'x');}
   /// v3/=v2.
   void operator  /= (const maT &op1)
      {array_by_array(*this,op1,*this,'/');}
   /// v3^=v2.
   void operator  ^= (const maT &op1)
      {array_by_array(*this,op1,*this,'^');}
   //@}

   /* Array "by" scalar ---------------------------------------------------- */
   /**@name Array "by" scalar operations
      These operations are between an array and a scalar (of the same type as
      the array). The result must
      have been defined to be of the same type as the operands.
      
      In this kind of operations each element of array 1 is operated with the
      given constant. The result has also got the same
      shape as the input array and its former content is lost
      */
      
   //@{
   /** This function must take one vector and a constant, and operate
      element by element according to the operation required. This is the
      function which really implements the operations. Simple calls to it
      perform much faster than calls to the corresponding operators.
      Although it is supposed to be a hidden function not useable by
      normal programmers.*/
   friend void array_by_scalar(const maT &op1, T op2, maT &result,
      char operation) {
         result.resize(op1);
         core_array_by_scalar(op1, op2, result, operation);
      }

   /** Core array by scalar operation.
       It assumes that the result is already resized. */
   friend void core_array_by_scalar<>(const maT &op1, const T &op2,
      maT &result, char operation);

   /// v3=v1+k
   maT operator  + (T op1) const 
      {maT temp; array_by_scalar(*this,op1,temp,'+'); return temp;}
   /// v3=v1-k
   maT operator  - (T op1) const
      {maT temp; array_by_scalar(*this,op1,temp,'-'); return temp;}
   /// v3=v1*k
   maT operator  * (T op1) const
      {maT temp; array_by_scalar(*this,op1,temp,'*'); return temp;}
   /// v3=v1/k
   maT operator  / (T op1) const
      {maT temp; array_by_scalar(*this,op1,temp,'/'); return temp;}
   /// v3=v1^k
   maT operator  ^ (T op1) const
      {maT temp; array_by_scalar(*this,op1,temp,'^'); return temp;}

   /// v3+=k
   void operator  += (const T &op1)
      {array_by_scalar(*this,op1,*this,'+');}
   /// v3-=k
   void operator  -= (const T &op1)
      {array_by_scalar(*this,op1,*this,'-');}
   /// v3*=k
   void operator  *= (const T &op1)
      {array_by_scalar(*this,op1,*this,'*');}
   /// v3/=k
   void operator  /= (const T &op1)
      {array_by_scalar(*this,op1,*this,'/');}
   /// v3^=k
   void operator  ^= (const T &op1)
      {array_by_scalar(*this,op1,*this,'^');}
   //@}

   /* Scalar "by" array ---------------------------------------------------- */
   /**@name Scalar "by" array operations
      These operations are between a scalar (of the same type as the array)
      and an array. The result must
      have been defined to be of the same type as the operand. The
      former content of the result array is lost after the operation.
      
      In this kind of operations the constant is operated with each element
      of array 2. The result has also got the same
      shape as the input array and its former content is lost
      */
   //@{
   /** This function must take one scalar and a vector, and operate
      element by element according to the operation required. This is the
      function which really implements the operations. Simple calls to it
      perform much faster than calls to the corresponding operators.
      Although it is supposed to be a hidden function not useable by
      normal programmers.*/
   friend void scalar_by_array(T op1, const maT &op2, maT &result,
      char operation) {
         result.resize(op2);
         core_scalar_by_array(op1, op2, result, operation);
   }

   /** Core array by scalar operation.
       It assumes that the result is already resized. */
   friend void core_scalar_by_array<>(const T &op1, const maT &op2,
      maT &result, char operation);

  /// v3=k+v2
   friend maT operator  + (T op1, const maT &op2)
      {maT temp; scalar_by_array(op1,op2,temp,'+'); return temp;}
   /// v3=k-v2
   friend maT operator  - (T op1, const maT &op2)
      {maT temp; scalar_by_array(op1,op2,temp,'-'); return temp;}
   /// v3=k*v2
   friend maT operator  * (T op1, const maT &op2)
      {maT temp; scalar_by_array(op1,op2,temp,'*'); return temp;}
   /// v3=k/v2
   friend maT operator  / (T op1, const maT &op2)
      {maT temp; scalar_by_array(op1,op2,temp,'/'); return temp;}
   /// v3=k^v2
   friend maT operator  ^ (T op1, const maT &op2)
      {maT temp; scalar_by_array(op1,op2,temp,'^'); return temp;}
   //@}
//@}

   /* Memory related ------------------------------------------------------- */
/**@name Size management */
//@{
   /// Core initialize. __m=NULL, __dim=0;
   void core_init() {__m=NULL; __dim=0;}

   /// Allocate vector memory
   void core_allocate(int _dim) {
      if (_dim<=0) {clear(); return;}
      __m = new T [_dim];
      if (__m != NULL) __dim = _dim;
      else REPORT_ERROR(1001,"Allocate: No space left");
   }

   /// Deallocate
   void core_deallocate() {
      if (__m!=NULL) delete [] __m;
      __m=NULL; __dim=0;
   }

   /** Resize according to a pattern.
       This function resize the actual array to the same size and origin
       as the input pattern. If the actual array is larger than the pattern
       then the trailing values are lost, if it is smaller then 0's are
       added at the end
       \\Ex: v2.resize(v1);
       \\--> v2 has got now the same structure as v1 */
   template <class T1>
      void resize(const maT1 &op1) {copy_shape(op1);}

   /** Get size.
       Returns the size of the object in a 3D vector. If the object
       is a matrix or a vector, then the higher order dimensions
       will be set to 1, ie, (Xdim, 1, 1) or (Xdim, Ydim, 1). */
   void get_size(int *size) const;

   /** Print shape of multidimensional array.
       This function shows the size, starting and finishing indexes of the
       given array. No end of line is printed neither at the beginning nor
       the end.
       \\Ex: v.print_shape();
       \\Ex: ofstream fh; ...; v.print_shape(fh); */
   void print_shape(ostream &out=cout) const;

   /** Makes an array empty.
       An empty array is the one with size 0, origin at 0, no statistics, ...
       The memory associated to the previous array is deallocated.
       \\Ex: v.clear(); */
   void clear() {core_deallocate(); init_shape();}

   /** Index outside logical definition region.
       This function returns TRUE if the position defined by v is
       outside the logical definition region of this array. It throws
       an exception if there are not enough components in v for this
       array, that is, 1 for vectors, 2 for matrices and 3 for volumes.*/
   // This function must be explictly implemented outside
   bool outside(const matrix1D<double> &v) const;

   /** True if this object intersects logically the argument array. */
   // This function must be explictly implemented outside
   bool intersects(const maT &m) const;

   /** True if this object intersects logically the argument rectangle.
       corner1 is the top-left corner and corner2 the right-bottom corner
       of teh rectangle.*/
   // This function must be explictly implemented outside
   bool intersects(const matrix1D<double> &corner1,
      const matrix1D<double> &corner2) const;

   /** True if the given index is in one corner.
       It is understood that one corner in logical coordinates. It throws
       an exception if there are not enough components in v for this
       array, that is, 1 for vectors, 2 for matrices and 3 for volumes. */
   // This function must be explictly implemented outside
   bool isCorner(const matrix1D<double> &v);

   /** True if the given index is in a border
       It is understood that one border in logical coordinates. It throws
       an exception if there are not enough components in v for this
       array, that is, 1 for vectors, 2 for matrices and 3 for volumes. */
   // This function must be explictly implemented outside
   bool isBorder(const matrix1D<int> &v);
   
   /** Insert a patch in array.
       This function substitutes the array given as argument in this object.
       It works with logical indexes, only the overlapping area is
       patched. You can modify the patch behaviour with the following
       operations '=', '+', '-', '*' or '/' */
   // This function must be explictly implemented outside
   void patch(const maT &patch_array, char operation='=');
//@}
   
   /* Initialisation ------------------------------------------------------- */
/**@name Initialisation */
//@{
   /** Same value in all components.
       The constant must be of a type compatible with the array type, ie,
       you cannot  assign a double to an integer array without a casting.
       It is not an error if the array is empty, then nothing is done.
       \\Ex: v.init_constant(3.14);*/
   void init_constant(T val)
       {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=val;}

   /** Initialize to zeros following a pattern.
       All values are set to 0, and the origin and size of the
       pattern are adopted.
       \\Ex: v2.init_zeros(v1); */
   void init_zeros(const maT &op) {resize(op); init_constant((T)0);}

   /** Initialize to zeros with current size.
       All values are set to 0. The current size and origin are kept. 
       It is not an error if the array is empty, then nothing is done.
       \\Ex: v.init_zeros();*/
   void init_zeros() {init_constant((T)0);}

   /** Initialize with random values.
       This function allows you to initialize the array with a set of
       random values picked from a uniform random distribution or a
       gaussian one. You must choose two parameters for each, for the
       uniform distribution they mean the range where to generate the
       random numbers, while in the gaussian case they are the mean and
       the standard deviation. By default the uniform distribution is
       selected. The size and origin of the array are not modified.
       \\Ex: v.init_random(0,1);
       \\--> uniform distribution between 0 and 1
       \\Ex: v.init_random(0,1,"uniform");
       \\--> the same
       \\Ex: v.init_random(0,1,"gaussian");
       \\--> gaussian distribution with 0 mean and stddev=1 */
   void init_random(double op1, double op2,
       const string &mode="uniform");

   /** Add noise to actual values.
       This function add some noise to the actual values of the array
       according to a certain random distribution. You must choose
       two parameters for each, for the uniform distribution they
       mean the range where to generate the random numbers, while
       in the gaussian case they are the mean and the standard
       deviation. By default the uniform distribution is
       selected. The size and origin of the array are not modified.
       The array itself is modified.
       \\Ex: v1.add_noise(0,1);
       \\--> uniform distribution between 0 and 1
       \\Ex: v1.add_noise(0,1,"uniform");
       \\--> the same
       \\Ex: v1.add_noise(0,1,"gaussian");
       \\--> gaussian distribution with 0 mean and stddev=1 */
   void add_noise(double op1, double op2,
       const string &mode="uniform") const;
//@}

   /* Some operators ------------------------------------------------------- */
/**@name Operators */
//@{
   /** Assignment.
       You can build as complex assignment expressions as you like. Multiple
       assignment is allowed.
       \\Ex: v1=v2+v3;
       \\Ex: v1=v2=v3; */
   maT& operator =  (const maT &op1) {
      if (&op1!=this) {
         resize(op1);
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=MULTIDIM_ELEM(op1,i);
      }
      return *this;
   }

   /** Unary plus.
       It is used to build arithmetic expressions.
       \\Ex: v1=+v2;*/
   maT operator +  () {return this;}

   /** Unary minus.
       It is used to build arithmetic expressions. You can make a minus
       of anything as long as it is correct semantically.
       \\Ex: v1=-v2;
       \\Ex: v1=-v2.transpose();*/
   maT  operator -  () const
      {maT temp(*this); temp.core_unary_minus(); return temp;}

   /** Core Unary minus.
       This function inverts (-x) each component of the core multidimensional
       array, the result is stored in this same array. */
   void core_unary_minus()
       {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=-mi;}

   /** Equality.
       Used to build logical expressions. It checks exact equality value by
       value (the difference mustn't be greater than
       \Ref{XMIPP_EQUAL_ACCURACY}), origins, and size.
       \\ Ex: if (v1==v2) cout << "v1 and v2 are the same";*/
   bool equal(const maT &op, double accuracy=XMIPP_EQUAL_ACCURACY) const {
      if (!same_shape(op)) return FALSE;
      else return core_equality(op,accuracy);
   }

   /** Operator ==.
       The XMIPP_EQUAL_ACCURACY is used as accuracy measure */
   bool friend operator ==<> (const maT &op1, const maT &op2);

   /** Operator !=.
       Not ==. */
   bool friend operator != (const maT &op1, const maT &op2)
     {return !(op1==op2);}

   /** Core equal.
       This is the equality of the core array. */
   bool core_equality(const maT &op2, double accuracy) const
      {if (MULTIDIM_SIZE(*this)!=MULTIDIM_SIZE(op2)) return FALSE;
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
          if (ABS(mi-MULTIDIM_ELEM(op2,i))>accuracy) return FALSE;
      return TRUE;}
   
   /** Input from input stream.
       Actual size of the array is used to know how many values must be read.
       \\Ex: v.resize(3); cin >> v;*/
   // This function must be explictly implemented outside
   friend istream& operator >> (istream& in, maT &v)
      {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v) in >> MULTIDIM_ELEM(v,i);
      return in;}
   
   /** Read from an ASCII file.
       The array must be previously resized to the correct size.*/
   void read(const FileName &fn)
      {ifstream  fh;
      fh.open(fn.c_str(), ios::in);
      if (!fh)
         REPORT_ERROR(1,(string)"MultidimArray::read: File "+fn+" not found");
      fh >> *this;
      fh.close();}

   /** Read from a binary file.
       The array must be previously resized to the correct size.*/
   void read_binary(const FileName &fn)
      {ifstream  fh;
      fh.open(fn.c_str(), ios::in | ios::binary);
      if (!fh)
         REPORT_ERROR(1,(string)"MultidimArray::read: File "+fn+" not found");
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         fh.read(mi, sizeof(T));
      fh.close();}

   /** Output to an output stream.
       A proper format is used to make all values fit in the same cell width.
       \\Ex: cout << v;*/
   // This function must be explictly implemented outside
   friend ostream& operator <<<>(ostream& out, const maT &v);

   /** Write to an ASCII file. */
   void write(const FileName &fn) const
      {ofstream  fh;
      fh.open(fn.c_str(), ios::out);
      if (!fh)
         REPORT_ERROR(1,(string)"MultidimArray::write: File "
            +fn+" cannot be openned for output");
      fh << *this;
      fh.close();}

   /** Write to a binary file. */
   void write_binary(const FileName &fn) const
      {ofstream  fh;
       
      fh.open(fn.c_str(), ios::out | ios::binary);
      if (!fh)
         REPORT_ERROR(1,(string)"MultidimArray::write: File "+fn+" not found");      
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this){	  	 	 
	 fh.write((void *) &mi, sizeof(T));		  	 	 
      }
      fh.close();}

    /** Edit with xmipp_editor.
        This function generates a random filename starting with PPP and
        edits it with xmipp_editor. After closing the editor the file is
        removed. */
    void edit() {
       FileName fn_random; fn_random.init_random(15);
       fn_random=(string)"PPP"+fn_random+".txt";
       write(fn_random);
       system(((string)"xmipp_edit -i "+fn_random+" -remove &").c_str());
    }
//@}
 
   /* Iterator ------------------------------------------------------------- */
/**@name Iterators
   Iterators are functions which could be applied to all elements, or
   all columns, or whatever in an array. There is only one iterator that
   can be applied to any kind of array. */
//@{
   /** Apply function to each element.
       You can define a function and apply it to each element of the array.
       A new array with the same shape as the input one is generated and the
       old one is not modified. The function applied must take an argument
       of type T and return a result of type T.
       \\Ex: v2=v1.for_all(&sin);
       \\--> If v1 and v2 are double arrays, then v2 is a version of v1 where
       all values of v1 have been substituted by their respective sines */
   maT for_all(T (*f)(T)) const
      {maT temp; temp.resize(*this);
       FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=(*f)(mi);
       return temp;}
      
   /** Apply a function to each element and store in this same object */
   void for_all(T (*f)(T))
      {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=(*f)(mi);}
//@}

   /* Utilities ------------------------------------------------------------ */
   /**@name Utilities
      Here you have several easy functions to manage the values of the array.*/
   //@{
   /** Dissimilarity.
       It returns the effective maximum absolute difference between 2
       arrays, ie, compute the absolute value of the difference of the
       arrays, then the histogram of the differences is computed and that
       value for which the number of higher values is less than the percentil
       out is returned. By default, only the 0.25% of the highest values are
       rejected as outliers, although this value could be increased as the
       number of samples in the array decreases.
       \\Ex: if (dissimilarity(v1,v2)<1e-6) cout << "v1 and v2 are more or less
              the same";
       \\--> the default 0.25% value is used
       \\Ex: if (dissimilarity(v1,v2,1)<1e-6) ...
       \\--> Now the 1% of the highest values are rejected to compute the
             effective range */
   friend double dissimilarity(maT &op1, maT &op2, double percentil_out)
      {maT diff; array_by_array(op1,op2,diff,'-'); diff.abs();
      //*** histogram1D hist; compute_hist(diff,hist,200);
      //return hist.percentil(100-percentil_out);
       return 1;}

   /** Several thresholding.
       Apply a threshold to the array, the object is modified itself. There
       are several kinds of thresholding and you must specify it, the values
       given in the fuction have different meanings according to the threshold
       applied.
       \\ abs_above: if |x|>a => x=b
       \\ abs_below: if |x|<a => x=b
       \\ above:     if  x >a => x=b
       \\ below:     if  x <a => x=b
       \\ range:     if  x <a => x=a   and    if x>b => x=b
       \\Ex: v.threshold("abs_above",10,10);
       \\--> any value whose absolute value is above 10 will be substituted by
             -10 (if it is negative) or 10 (if it is positive)
       \\Ex: v.threshold("abs_below",0.1,0);
       \\--> any value whose absolute value is below 0.1 will be substituted by
             -0 (if it is negative) or 0 (if it is positive)
       \\Ex: v.threshold("above",10,10);
       \\--> any value above 10 will be substituted by 10
       \\Ex: v.threshold("below",-10,-10);
       \\--> any value below -10 will be substituted by -10
       \\Ex: v.threshold("range",0,1);
       \\--> v is "saturated" by values 0 and 1, any value outside this range
             will be substituted by its nearest border*/
   void threshold(const string &type, T a, T b);

   /** Count with threshold.
       This function returns the number of elements meeting the threshold
       condition. */
   long count_threshold(const string &type, T a, T b);

   /** Substitute a value by another.
       Substitute an old value by a new one. The accuracy is used to
       say if the value in the array is equal to the old value. Set it
       to 0 for perfect accuracy. */
   void substitute(T old_value, T new_value,
      double accuracy=XMIPP_EQUAL_ACCURACY)
      {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         if (ABS(mi-old_value)<=accuracy) mi=new_value;}

   /** Binarize.
       This functions substitutes all values in a volume which are greater
       than val+accuracy by 1 and the rest are set to 0.
       Use threshold to get a very powerful binarization */
   void binarize(double val=0, double accuracy=XMIPP_EQUAL_ACCURACY)
      {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this)
         if (mi<=val+accuracy) mi=0; else mi=1;}

   /** ROUND n-dimensional.
       Applies a ROUND (look for the nearest integer) to each array element.*/
   void ROUNDnD()
     {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=ROUND(mi);}

   /** ROUND n-dimensional.
       The same as before but the result is returned. */
   friend maT ROUNDnD(const maT &a)
     {maT temp(a); temp.ROUNDnD(); return temp;}

   /** CEILING n-dimensional.
       Applies a CEILING (look for the nearest larger integer) to each
       array element.*/
   void CEILnD()
     {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=CEIL(mi);}

   /** CEILING n-dimensional.
       The same as before but the result is returned. */
   friend maT CEILnD(const maT &a)
     {maT temp(a); temp.CEILnD(); return temp;}

   /** FLOOR n-dimensional.
       Applies a FLOOR (look for the nearest larger integer) to each
       array element.*/
   void FLOORnD()
     {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=FLOOR(mi);}

   /** FLOOR n-dimensional.
       The same as before but the result is returned. */
   friend maT FLOORnD(const maT &a)
     {maT temp(a); temp.FLOORnD(); return temp;}

   /** ABS n-dimensional.
       Applies an ABS (absolute value) to each array element.*/
   void ABSnD()
     {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=ABS(mi);}

   /** ABS n-dimensional.
       The same as before but the result is returned. */
   friend maT ABSnD(const maT &a)
     {maT temp(a); temp.ABSnD(); return temp;}

   /** MAX n-dimensional.
       Each component of the result is the maximum of the correspoing
       components of the two input arrays. They must have the same shape, if
       not an exception is thrown */
   friend void MAXnD(const maT &v1, const maT &v2, maT &result)
      {if (!v1.same_shape(v2))
         REPORT_ERROR(1007,"MAX: arrays of different shape");
      result.resize(v1);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         MULTIDIM_ELEM(result,i)=MAX(MULTIDIM_ELEM(v1,i),MULTIDIM_ELEM(v2,i));}

   /** MAX n-dimensional.
       The same as before but the result is returned. */
   friend maT MAXnD(const maT &v1, const maT &v2)
      {maT temp; MAXnD(v1,v2,temp); return temp;}

   /** MIN n-dimensional.
       Each component of the result is the minimum of the correspoing
       components of the two input arrays. They must have the same shape, if
       not an exception is thrown */
   friend void MINnD(const maT &v1, const maT &v2, maT &result)
      {if (!v1.same_shape(v2))
         REPORT_ERROR(1007,"MIN: arrays of different shape");
      result.resize(v1);
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         MULTIDIM_ELEM(result,i)=MIN(MULTIDIM_ELEM(v1,i),MULTIDIM_ELEM(v2,i));}

   /** MIN n-dimensional.
       The same as before but the result is returned. */
   friend maT MINnD(const maT &v1, const maT &v2)
      {maT temp; MINnD(v1,v2,temp); return temp;}

   /** Sqrt.
       Each component of the result is the square root of the original
       components.*/
   void SQRTnD()
     {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=(T)sqrt((double)mi);}

   /** Sqrt n-dimensional.
       The same as before but the result is returned. */
   friend maT SQRTnD(const maT &a)
     {maT temp(a); temp.SQRTnD(); return temp;}

   /** Sum of matrix values.
       This function returns the sum of all internal values.
       \\Ex: double sum=m.sum();*/
   double sum() const
         {double accumulate=0;
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) accumulate += mi;
         return accumulate;}

   /** Sum of squared vector values.
       This function returns the sum of all internal values to the second
       power.
       \\Ex: double sum2=m.sum2();*/
   double sum2() const
        {double accumulate=0;
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) accumulate += mi*mi;
         return accumulate;}

   /** Log10.
       Each component of the result is the log10 of the original
       components.*/
   void self_log10()
     {FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(*this) mi=(T)log10((double)mi);}

   /** Log10 n-dimensional.
       The same as before but the result is returned. */
   friend maT log10(const maT &a)
     {maT temp(a); temp.self_log10(); return temp;}

   /** Compute center of mass.
       If a mask is provided it must be of the same dimension of the object 
       and of type int (i.e., matrix2D<int> *). Only those logical indexes
       within the object whose mask value is 1 (also in logical indexes)
       are taken into account.*/
   void center_of_mass(matrix1D<double> &center, void * mask=NULL);

   //@}
#undef mi
#undef msize
