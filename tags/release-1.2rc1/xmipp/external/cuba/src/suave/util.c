/*
	util.c
		Utility functions
		this file is part of Suave
		last modified 9 Feb 05 th
*/


#include "decl.h"

static count ndim_, ncomp_, nregions_;
static number neval_;


#define RegionAlloc(p, n, nnew) \
  MemAlloc(p, sizeof(Region) + \
              (n)*(ndim_ + ncomp_ + 1)*sizeof(real) + \
              (nnew)*ndim_*sizeof(bin_t))


#ifdef DEBUG
#include "debug.c"
#endif

