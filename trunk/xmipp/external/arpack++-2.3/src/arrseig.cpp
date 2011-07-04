// -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-
/*
  Implementation of MemoryOverflow
  
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
#include <arrseig.h>


void
MemoryOverflow() throw()
{ 
    throw ArpackError(ArpackError::MEMORY_OVERFLOW); 
}

