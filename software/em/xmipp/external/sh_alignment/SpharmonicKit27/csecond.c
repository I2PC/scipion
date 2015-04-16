/***************************************************************************
  **************************************************************************
  
                Spherical Harmonic Transform Kit 2.7
  
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
  
   Copyright 1997-2003  Sean Moore, Dennis Healy,
                        Dan Rockmore, Peter Kostelec
  
  
   Copyright 2004  Peter Kostelec, Dan Rockmore


     SpharmonicKit is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
  
     SpharmonicKit is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
  
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  
  
   Commercial use is absolutely prohibited.
  
   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


/*

  a function used for timing; provided by Mark Taylor

      real*8 second
      before = second()
      code_to_time
      after = second()
      cpu_time = after - before

 *----------------------------------------------------------------------
 * csecond.c
 *----------------------------------------------------------------------
 * Last written:
 * Time-stamp: "1995/11/10 16:52:56 thibaud@kether.cgd.ucar.edu"
 *----------------------------------------------------------------------
 */

#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>


#ifdef CLK_TCK
#define DIVIDER CLK_TCK
#else
#define DIVIDER HZ         /* for old BSD systems */
#endif




/* define this to return wallclock instead of cpu time */
/* #define WALLCLOCK */


#if (defined CRAY)
double CSECOND( )   /* should this really be double??? */
#elif (defined T3D)
double CSECOND( ) 
#elif (defined IBM )
double csecond( )
#else /* works for SUN, LINUX, DECALPHA, SGI */
double csecond( )
#endif
{
    struct tms buf;
    static struct tms buf0;			/* times structure */
    static int firstcall=1;
    clock_t  rv;
    static clock_t rv0;

    if (firstcall) {
       firstcall=0;
       rv0=times(&buf0);
    }

    rv=times(&buf);
    /*    if ( rv < 0 ) {fprintf( stderr, "csecond failed  %d \n",rv ); } */
#ifdef WALLCLOCK
    return( (double) (rv-rv0) / (double)DIVIDER );
#else
    return( (double)( (buf.tms_utime + buf.tms_stime + buf.tms_cutime + buf.tms_cstime) -
                      (buf0.tms_utime + buf0.tms_stime + buf0.tms_cutime + buf0.tms_cstime)
                     ) / (double)DIVIDER );
#endif
}











