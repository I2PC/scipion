/* Speed up ------------------------------------------------------------- */
/**@name Macros for all multidimensional arrays
   This macros are defined to allow high speed in critical parts of
   your program. They shouldn't be used systemultiTically as usually there
   is no checking on the correctness of the operation you are performing.
   Speed comes from three facts: first, they are macros and no function
   call is performed (although most of the critical functions are
   inline functions), there is no checking on the correctness of the
   operation (it could be wrong and you are not warned of it), and
   destination vectors are not returned saving time in the copy
   constructor and in the creation/destruction of temporary vectors.*/
//@{

/** Access to dimension (size).*/
#ifndef MULTIDIM_SIZE
   #define MULTIDIM_SIZE(v) ((v).__dim)
#endif

/** Dimension of the space in which it is defined (1D, 2D or 3D).*/
#ifndef SPACE_DIM
   #define SPACE_DIM(v) ((v).__spcdim)
#endif

/** For all elements in the array.
    This macro is used to generate loops for the array in an easy manner.
    It defines an internal index 'i' which ranges the array using its
    physical definition
    \\Ex:
    \begin{verbatim}
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v) 
       cout << DIRECT_MULTIDIM_ELEM(v,i) << " ";
    \end{verbatim}*/
#ifndef FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY
   #define FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v) \
       for (int i=0; i<(v).__dim; i++)
#endif

/** Vector element: Physical access.
    Be careful because this is physical access, usually vectors follow
    the C convention of starting index==0.
    \\ Ex: MULTIDIM_ELEM(v,0)=1;
    \\ Ex: val=MULTIDIM_ELEM(v,0);*/
#ifndef MULTIDIM_ELEM
   #define MULTIDIM_ELEM(v,i) ((v).__m[(i)])
#endif

/** Array access.
    This macro gives you access to the array (T *).
    \\ Ex: cout << "This is an int *" << MULTIDIM_ARRAY(v) << endl; */
#ifndef MULTIDIM_ARRAY
   #define MULTIDIM_ARRAY(v) ((v).__m)
#endif
//@}

/* Other useful functions -------------------------------------------------- */
/**@name Functions for all multidimensional arrays */
//@{
/** Conversion from one type to another.
    If we have an integer vector and we need a double one, we can use this
    function. The conversion is done through a type casting of each element.*/
template <class T, class T1>
   void type_cast(const maT &v1, maT1 &v2) {
         if (MULTIDIM_SIZE(v1)==0) {v2.clear(); return;}
         v2.resize(v1);
         FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v2)
            MULTIDIM_ELEM(v2,i)=(T1) MULTIDIM_ELEM(v1,i);
      }
//@}
