// -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-
/*
  Implementation of arpack_debug
  
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
#include <arpackf.h>
#include <debug.h>


extern "C"
{
    arpack_debug F77NAME(debug);
}
