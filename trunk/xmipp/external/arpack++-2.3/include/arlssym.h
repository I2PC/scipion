/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLSSym.h.
   Arpack++ class ARluSymStdEig definition
   (SuperLU version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Kristi Maschhoff
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLSSYM_H
#define ARLSSYM_H

#include <stddef.h>
#include "arch.h"
#include "arssym.h"
#include "arlsmat.h"


template<class ARFLOAT>
class ARluSymStdEig:
  public virtual ARSymStdEig<ARFLOAT, ARluSymMatrix<ARFLOAT> > {

 protected:

 // a) Protected function:

  virtual void Copy(const ARluSymStdEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // b) Public functions:

 // b.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // b.2) Constructors and destructor.

  ARluSymStdEig() { }
  // Short constructor.

  ARluSymStdEig(int nevp, ARluSymMatrix<ARFLOAT>& A,
                char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluSymStdEig(int nevp, ARluSymMatrix<ARFLOAT>& A,
                ARFLOAT sigma, char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluSymStdEig(const ARluSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluSymStdEig() { }
  // Destructor.

 // c) Operators.

  ARluSymStdEig& operator=(const ARluSymStdEig& other);
  // Assignment operator.

}; // class ARluSymStdEig.


// ------------------------------------------------------------------------ //
// ARluSymStdEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluSymStdEig<ARFLOAT>::Copy(const ARluSymStdEig<ARFLOAT>& other)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARluSymMatrix<ARFLOAT> >:: Copy(other);
  if (mode > 2) objOP->FactorAsI(sigmaR);

} // Copy.


template<class ARFLOAT>
inline void ARluSymStdEig<ARFLOAT>::ChangeShift(ARFLOAT sigmaRp)
{

  sigmaR    = sigmaRp;
  sigmaI    = 0.0;
  mode      = 3;
  iparam[7] = mode;

  objOP->FactorAsI(sigmaR);
  Restart();

} // ChangeShift.


template<class ARFLOAT>
inline void ARluSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARluSymMatrix<ARFLOAT> >::
    SetRegularMode(objOP, &ARluSymMatrix<ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARluSymMatrix<ARFLOAT> >::
    SetShiftInvertMode(sigmap, objOP, &ARluSymMatrix<ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARluSymStdEig<ARFLOAT>::
ARluSymStdEig(int nevp, ARluSymMatrix<ARFLOAT>& A,
              char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)
{

  NoShift();
  DefineParameters(A.ncols(), nevp, &A, &ARluSymMatrix<ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluSymStdEig<ARFLOAT>::
ARluSymStdEig(int nevp, ARluSymMatrix<ARFLOAT>& A,
              ARFLOAT sigmap, char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A, &ARluSymMatrix<ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARluSymStdEig<ARFLOAT>& ARluSymStdEig<ARFLOAT>::
operator=(const ARluSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLSSYM_H
