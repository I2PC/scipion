/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARError.h.
   Definition of ArpackError, a class that handles errors
   occurred during Arpack execution.

   There are three ways of handling an error:
   a) Declaring a variable of type ArpackError and calling
      function Set with the correct ErrorCode (see codes below).
   b) Calling the constructor ArpackError(ErrorCode) to define
      a variable.
   c) Calling ArpackError::Set(ErrorCode) directly.

   If an error occurs, a brief description of the error is
   displayed on the "cerr" stream, unless the variable
   ARPACK_SILENT_MODE is defined.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARERROR_H
#define ARERROR_H

#include <iostream>
#include <cstdlib>


#include "arch.h"

class ArpackError {

 public:

  enum ErrorCode {  // Listing all kinds of errors.

    // Innocuous error type.

    NO_ERRORS           =    0,

    // Errors in parameter definitions.

    PARAMETER_ERROR     = -101,
    N_SMALLER_THAN_2    = -102,
    NEV_OUT_OF_BOUNDS   = -103,
    WHICH_UNDEFINED     = -104,
    PART_UNDEFINED      = -105,
    INVMODE_UNDEFINED   = -106,
    RANGE_ERROR         = -107,

    // Errors in Aupp and Eupp functions.

    LAPACK_ERROR        = -201,
    START_RESID_ZERO    = -202,
    NOT_ACCURATE_EIG    = -203,
    REORDERING_ERROR    = -204,
    ARNOLDI_NOT_BUILD   = -205,
    AUPP_ERROR          = -291,
    EUPP_ERROR          = -292,

    // Errors in main functions.

    CANNOT_PREPARE      = -301,
    CANNOT_FIND_BASIS   = -302,
    CANNOT_FIND_VALUES  = -303,
    CANNOT_FIND_VECTORS = -304,
    CANNOT_FIND_SCHUR   = -305,
    SCHUR_UNDEFINED     = -306,

    // Errors due to incorrect function calling sequence.

    CANNOT_GET_VECTOR   = -401,
    CANNOT_GET_PROD     = -402,
    CANNOT_PUT_VECTOR   = -403,
    PREPARE_NOT_OK      = -404,
    BASIS_NOT_OK        = -405,
    VALUES_NOT_OK       = -406,
    VECTORS_NOT_OK      = -407,
    SCHUR_NOT_OK        = -408,
    RESID_NOT_OK        = -409,

    // Errors in classes that perform LU decompositions.

    MATRIX_IS_SINGULAR  = -501,
    DATA_UNDEFINED      = -502,
    INSUFICIENT_MEMORY  = -503,
    NOT_SQUARE_MATRIX   = -504,
    NOT_FACTORED_MATRIX = -505,
    INCOMPATIBLE_SIZES  = -506,
    DIFFERENT_TRIANGLES = -507,
    INCONSISTENT_DATA   = -508,
    CANNOT_READ_FILE    = -509,

    // Errors in matrix files.

    CANNOT_OPEN_FILE    = -551, 
    WRONG_MATRIX_TYPE   = -552,
    WRONG_DATA_TYPE     = -553,
    RHS_IGNORED         = -554,
    UNEXPECTED_EOF      = -555,

    // Other severe errors.

    NOT_IMPLEMENTED     = -901,
    MEMORY_OVERFLOW     = -902,
    GENERIC_SEVERE      = -999,

    // Warnings.

    NCV_OUT_OF_BOUNDS   =  101,
    MAXIT_NON_POSITIVE  =  102,
    MAX_ITERATIONS      =  201,
    NO_SHIFTS_APPLIED   =  202,
    CHANGING_AUTOSHIFT  =  301,
    DISCARDING_FACTORS  =  401,
    GENERIC_WARNING     =  999

  };

 private:

  static ErrorCode code;

  static void Print(const char* where, const char* message);
  // Writes error messages on cerr stream.

 public:

  static void Set(ErrorCode error, const char* where="AREigenProblem");
  // Set error code and write error messages.

  static int Status() { return (int) code; }
  // Returns current value of error code.

  ArpackError(ErrorCode error, char const* where="AREigenProblem") {
    Set(error,where);
  }
  // Constructor that set error code.

  ArpackError() { code = NO_ERRORS; };
  // Constructor that does nothing.

};

#endif // ARERROR_H
