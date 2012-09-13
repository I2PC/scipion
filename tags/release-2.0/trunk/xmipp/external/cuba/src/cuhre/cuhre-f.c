/*
	cuhre-f.c
		Fortran interface for Cuhre
		this file is part of Cuhre
		last modified 1 Mar 06 th
*/


#include "decl.h"

Extern void PREFIX(Cuhre)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob);


Extern void USCORE(cuhre)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  ccount *pkey,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  PREFIX(Cuhre)(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey,
    pnregions, pneval, pfail,
    integral, error, prob);
}

