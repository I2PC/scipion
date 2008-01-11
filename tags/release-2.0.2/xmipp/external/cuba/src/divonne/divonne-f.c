/*
	divonne-f.c
		Fortran interface for Divonne
		this file is part of Divonne
		last modified 1 Mar 06 th
*/


#include "decl.h"

Extern void PREFIX(Divonne)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  creal border, creal maxchisq, creal mindeviation,
  cnumber ngiven, ccount ldxgiven, real *xgiven,
  cnumber nextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob);


Extern void USCORE(divonne)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cint *pkey1, cint *pkey2, cint *pkey3, ccount *pmaxpass,
  creal *pborder, creal *pmaxchisq, creal *pmindeviation,
  cnumber *pngiven, ccount *pldxgiven, real *xgiven,
  cnumber *pnextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  PREFIX(Divonne)(*pndim, *pncomp,
    integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey1, *pkey2, *pkey3, *pmaxpass,
    *pborder, *pmaxchisq, *pmindeviation,
    *pngiven, *pldxgiven, xgiven,
    *pnextra, peakfinder,
    pnregions, pneval, pfail,
    integral, error, prob);
}

