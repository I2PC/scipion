/**@defgroup FindRoot Find root
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** Find root by bracketing.
  * Bracketing of a root of (Function) inside the interval
  * [LowerBound, UpperBound]. The purpose is to search for a pair of
  * arguments for which 'Function' differs in sign. The search is conducted
  * within [LowerBound, UpperBound].
  * At most [(UpperBound - LowerBound) / Tolerance] function evaluations
  * are performed. Even-order roots cannot be bracketed.
  * The evaluation of Function at LowerBound (UpperBound) is returned in
  * LowerSample (UpperSample)
  *
  * success: return(!ERROR); failure: return(ERROR)
  *
  * the function 'Function' must be declared as follows:
  * @code
  *    extern int myFunction(double myArgument, void *AuxilliaryData, double *myResult);
  * @endcode
  * It must return ERROR upon failure, and !ERROR upon success
  *
  * It is evaluated for the value of the variable 'myArgument'. The result
  * of the function evaluation must be returned in 'myResult'. The
  * generic pointer 'AuxilliaryData' can be used to pass additional parameters.
  *
  * What follows is a developed example of the function
  *    f(x) = a * x^2 + b * x + c (a, b, and c are free parameters)
  *
  * @code
  *   struct myStruct {
  *      double a, b, c;
  *   };
  *
  *   extern int myFunction (
  *      double  myArgument,
  *      void    *AuxilliaryData,
  *      double  *myResult) {
  *      struct myStruct myData;
  *
  *      myData = *((struct myStruct *)AuxilliaryData);
  *      *myResult = myArgument * (myArgument * myData.a + myData.b) + myData.c;
  *      return(!ERROR);
  *   }
  *
  *   int main() {
  *      struct myStruct myData;
  *      double LowerBound = -100.0, UpperBound = 100.0;
  *      double LowerSample, UpperSample;
  *      double Tolerance = FLT_EPSILON;
  *      int    ValidBracket;
  *      int    Status;
  *
  *      myData.a = 1.0;
  *      myData.b = 5.0;
  *      myData.c = 4.0;
  *      RootBracket(*myFunction, (void *)&myData, &LowerBound, &UpperBound,
  *         &LowerSample, &UpperSample, Tolerance, &ValidBracket, &Status);
  *      return(0);
  *   }
  *   @endcode
*/
extern int  RootBracket
(
    int(*Function)(double, void *, double *),
    /* function to bracket */
    void *AuxilliaryData, /* parameters used by the function */
    double *LowerBound,  /* lower interval bound to be updated */
    double *UpperBound,  /* upper interval bound to be updated */
    double *LowerSample,  /* value of Function for argument LowerBound */
    double *UpperSample,  /* value of Function for argument UpperBound */
    double Tolerance,   /* admissible relative error */
    int  *ValidBracket,  /* whether or not a root could be bracketed */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Find a root by bisection.
  * Search for a root of (Function) inside the bracketing interval
  * [LowerBound, UpperBound]. The strategy proceeds by iteratively
  * cutting the interval in two equal parts. Even-order roots generally
  * cannot be found. Only one root is returned, even if there are several
  * ones. Tolerance is relative to the size of the bracketing interval
  *
  * success: return(!ERROR); failure: return(ERROR)
  *
  * The function 'Function' must be declared as follows:
  * @code
  *    extern int myFunction(double myArgument, void *AuxilliaryData, double *myResult);
  * @endcode
  * It must return ERROR upon failure, and !ERROR upon success. It is
  * evaluated for the value of the variable 'myArgument'. The result of the
  * function evaluation must be returned in 'myResult'. The generic pointer
  * 'AuxilliaryData' can be used to pass additional parameters.
  *
  * What follows is a developed example of the function
  * f(x) = a * x^2 + b * x + c (a, b, and c are free parameters)
  * @code
  *   struct myStruct {
  *      double a, b, c;
  *   };
  *
  *   extern int myFunction (
  *      double  myArgument,
  *      void    *AuxilliaryData,
  *      double  *myResult) {
  *      struct myStruct myData;
  *
  *      myData = *((struct myStruct *)AuxilliaryData);
  *      *myResult = myArgument * (myArgument * myData.a + myData.b) + myData.c;
  *      return(!ERROR);
  *   }
  *
  *   int main() {
  *      struct myStruct myData;
  *      double LowerBound = -100.0, UpperBound = 100.0;
  *      double Root;
  *      double Tolerance = FLT_EPSILON;
  *      int    Status;
  *
  *      myData.a = 1.0;
  *      myData.b = 5.0;
  *      myData.c = 4.0;
  *      RootBracket(*myFunction, (void *)&myData, &LowerBound, &UpperBound, Tolerance, &Status);
  *      RootFindBisection(*myFunction, (void *)&myData, &Root, LowerBound,
  *         UpperBound, Tolerance, &Status);
  *      return(0);
  *   }
  *   @endcode
*/
extern int  RootFindBisection
(
    int(*Function)(double, void *, double *),
    /* function, of which a root is sought */
    void *AuxilliaryData, /* parameters used by the function */
    double *Root,    /* returned root */
    double LowerBound,   /* lower bound of an interval containing a root */
    double UpperBound,   /* upper bound of an interval containing a root */
    double Tolerance,   /* admissible relative error */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Find root using Brent algorithm.
  * Search for a root of (Function) inside the bracketing interval
  * [LowerBound, UpperBound]. The strategy proceeds by fitting an inverse
  * quadratic function and by using the root that belong to the interval.
  * If any even-order roots generally cannot be found.
  * Only one root is returned, even if there are several ones.
  * Tolerance is relative to the size of the bracketing interval.
  *
  * success: return(!ERROR); failure: return(ERROR)
  *
  * The function 'Function' must be declared as follows:
  * @code
  * extern int myFunction(double myArgument, void *AuxilliaryData, double *myResult);
  * @endcode
  *
  * It must return ERROR upon failure, and !ERROR upon success. It is
  * evaluated for the value of the variable 'myArgument'. The result of the
  * function evaluation must be returned in 'myResult'. The generic pointer
  * 'AuxilliaryData' can be used to pass additional parameters.
  *
  * What follows is a developed example of the function
  * f(x) = a * x^2 + b * x + c (a, b, and c are free parameters)
  * @code
  *   struct myStruct {
  *      double a, b, c;
  *   };
  *
  *   extern int myFunction (
  *      double  myArgument,
  *      void    *AuxilliaryData,
  *      double  *myResult) {
  *      struct myStruct myData;
  *
  *      myData = *((struct myStruct *)AuxilliaryData);
  *      *myResult = myArgument * (myArgument * myData.a + myData.b) + myData.c;
  *      return(!ERROR);
  *   }
  *
  *   int main() {
  *      struct myStruct myData;
  *      double LowerBound = -100.0, UpperBound = 100.0;
  *      double Root;
  *      double Tolerance = FLT_EPSILON;
  *      int    Status;
  *
  *      myData.a = 1.0;
  *      myData.b = 5.0;
  *      myData.c = 4.0;
  *      RootBracket(*myFunction, (void *)&myData, &LowerBound, &UpperBound, Tolerance, &Status);
  *      RootFindBisection(*myFunction, (void *)&myData, &Root, LowerBound,
  *         UpperBound, Tolerance, &Status);
  *      return(0);
  *   }
  *   @endcode
*/
extern int  RootFindBrent
(
    int(*Function)(double, void *, double *),
    /* function, of which a root is sought */
    void *AuxilliaryData, /* parameters used by the function */
    double *Root,    /* returned root */
    double LowerBound,   /* lower bound of an interval containing a root */
    double UpperBound,   /* upper bound of an interval containing a root */
    double Tolerance,   /* admissible relative error */
    int  *Status    /* error management */
);
//@}
