// -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-
/*
  Implementation of ArpackError

  ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
  
  Copyright (C) 2003 Christophe Prud'homme (prudhomm@debian.org)

  Public domain software.
*/
#include <arerror.h>

// "code" initialization.
ArpackError::ErrorCode ArpackError::code = NO_ERRORS;

void
ArpackError::Print(const char* where, const char* message)
{

#ifndef ARPACK_SILENT_MODE
  std::cerr << "Arpack error in " << where << "." << std::endl;
  std::cerr << "-> " << message << "." << std::endl;
#endif

} // Print

void
ArpackError::Set(ErrorCode error, char const* where)
{

  code = error;
  switch (code) {
  case NO_ERRORS          :
    return;
  case NOT_IMPLEMENTED    :
    Print(where, "This function was not implemented yet");
    return;
  case MEMORY_OVERFLOW    :
    Print(where, "Memory overflow");
    return;
  case GENERIC_SEVERE     :
    Print(where, "Severe error");
    return;
  case PARAMETER_ERROR    :
    Print(where, "Some parameters were not correctly defined");
    return;
  case N_SMALLER_THAN_2   :
    Print(where, "'n' must be greater than one");
    return;
  case NEV_OUT_OF_BOUNDS  :
    Print(where, "'nev' is out of bounds");
    return;
  case WHICH_UNDEFINED    :
    Print(where, "'which' was not correctly defined");
    return;
  case PART_UNDEFINED     :
    Print(where, "'part' must be one of 'R' or 'I'");
    return;
  case INVMODE_UNDEFINED  :
    Print(where, "'InvertMode' must be one of 'S' or 'B'");
    return;
  case RANGE_ERROR  :
    Print(where, "Range error");
    return;
  case LAPACK_ERROR       :
    Print(where, "Could not perform LAPACK eigenvalue calculation");
    return;
  case START_RESID_ZERO   :
    Print(where, "Starting vector is zero");
    return;
  case NOT_ACCURATE_EIG   :
    Print(where, "Could not find any eigenvalue to sufficient accuracy");
    return;
  case REORDERING_ERROR   :
    Print(where, "Reordering of Schur form was not possible");
    return;
  case ARNOLDI_NOT_BUILD  :
    Print(where, "Could not build an Arnoldi factorization");
    return;
  case AUPP_ERROR         :
    Print(where, "Error in ARPACK Aupd fortran code");
    return;
  case EUPP_ERROR         :
    Print(where, "Error in ARPACK Eupd fortran code");
    return;
  case CANNOT_PREPARE     :
    Print(where, "Could not correctly define internal variables");
    return;
  case CANNOT_FIND_BASIS  :
    Print(where, "Could not find the Arnoldi basis vectors");
    return;
  case CANNOT_FIND_VALUES :
    Print(where, "Could not find any eigenvalue");
    return;
  case CANNOT_FIND_VECTORS:
    Print(where, "Could not find any eigenvector");
    return;
  case CANNOT_FIND_SCHUR  :
    Print(where, "Could not find any Schur vector");
    return;
  case SCHUR_UNDEFINED    :
    Print(where, "FindEigenvectors must be used instead of FindSchurVectors");
    return;
  case CANNOT_GET_VECTOR  :
    Print(where, "Vector is not already available");
    return;
  case CANNOT_GET_PROD    :
    Print(where, "Matrix vector product is not already available");
    return;
  case CANNOT_PUT_VECTOR  :
    Print(where, "Could not store vector");
    return;
  case PREPARE_NOT_OK     :
    Print(where, "DefineParameters must be called prior to this function");
    return;
  case BASIS_NOT_OK       :
    Print(where, "Arnoldi basis is not available");
    return;
  case VALUES_NOT_OK      :
    Print(where, "Eigenvalues are not available");
    return;
  case VECTORS_NOT_OK     :
    Print(where, "Eigenvectors are not available");
    return;
  case SCHUR_NOT_OK       :
    Print(where, "Schur vectors are not available");
    return;
  case RESID_NOT_OK       :
    Print(where, "Residual vector is not available");
    return;
  case MATRIX_IS_SINGULAR :
    Print(where, "Matrix is singular and could not be factored");
    return;
  case DATA_UNDEFINED     :
    Print(where, "Matrix data was not defined");
    return;
  case INSUFICIENT_MEMORY :
    Print(where, "fill-in factor must be increased");
    return;
  case NOT_SQUARE_MATRIX  :
    Print(where, "Matrix must be square to be factored");
    return;
  case NOT_FACTORED_MATRIX:
    Print(where, "Matrix must be factored before solving a system");
    return;
  case INCOMPATIBLE_SIZES :
    Print(where, "Matrix dimensions must agree");
    return;
  case DIFFERENT_TRIANGLES:
    Print(where, "A.uplo and B.uplo must be equal");
    return;
  case INCONSISTENT_DATA  :
    Print(where, "Matrix data contain inconsistencies");
    return;
  case CANNOT_READ_FILE   :
    Print(where, "Data file could not be read");
    return;
  case CANNOT_OPEN_FILE   :
    Print(where, "Invalid path or filename");
    return;
  case WRONG_MATRIX_TYPE  :
    Print(where, "Wrong matrix type");
    return;
  case WRONG_DATA_TYPE    :
    Print(where, "Wrong data type");
    return;
  case RHS_IGNORED        :
    Print(where, "RHS vector will be ignored");
    return;
  case UNEXPECTED_EOF     :
    Print(where, "Unexpected end of file");
    return;
  case NCV_OUT_OF_BOUNDS  :
    Print(where, "'ncv' is out of bounds");
    return;
  case MAXIT_NON_POSITIVE :
    Print(where, "'maxit' must be greater than zero");
    return;
  case MAX_ITERATIONS     :
    Print(where, "Maximum number of iterations taken");
    return;
  case NO_SHIFTS_APPLIED  :
    Print(where, "No shifts could be applied during a cycle of IRAM iteration");
    return;
  case CHANGING_AUTOSHIFT :
    Print(where, "Turning to automatic selection of implicit shifts");
    return;
  case DISCARDING_FACTORS :
    Print(where, "Factors L and U were not copied. Matrix must be factored");
    return;
  case GENERIC_WARNING    :
  default: ;
    Print(where, "There is something wrong");
    return;
  }

} // Set.

