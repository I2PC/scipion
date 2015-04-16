/*********************************************************************
*                          L I B _ S T D                             *
**********************************************************************
* Library is part of the Situs package (c) Willy Wriggers, 1998-2003 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    * 
* Auxiliary routines to read from stdin.                             *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/


#include "situs.h"    

int readlnint ()
/* returns first integer from input stream and checks for EOF or errors */ 

{ int i, ch;

  if (scanf("%d", &i) == EOF ) {
    fprintf(stderr, "lib_std> Error: EOF while reading input [e.c. 14010]\n");
    exit (14010);
  }  
  for ( ; ; )  {
    ch = getchar(); 
    if (ch == EOF ) {
      fprintf(stderr, "lib_std> Error: EOF while reading input [e.c. 14020]\n");
      exit (14020);
    }
    if (ch == '\n') break;
  } 
  return i;
}

double readlnmyflt ()
/* returns first floating point from input stream and checks for EOF or errors */ 

{ double ddx;
  int ch;

  if (scanf("%le", &ddx) == EOF ) {
    fprintf(stderr, "lib_std> Error: EOF while reading input [e.c. 14030]\n");
    exit (14030);
  }  
  for ( ; ; )  {
    ch = getchar(); 
    if (ch == EOF ) {
      fprintf(stderr, "lib_std> Error: EOF while reading input [e.c. 14040]\n");
      exit (14040);
    }
    if (ch == '\n') break;
  } 
  return ((double)ddx);
}

void removespaces(char *spaceyname, int fl) {
/* removes erroneously pasted spaces from file name string */

  int i;

  spaceyname[fl-1]='\0';
  
  /* remove leading space */
  for(;;) {
    if (spaceyname[0]==' ') {
      for (i=1;i<fl;++i) spaceyname[i-1]=spaceyname[i]; 
      --fl;
    }
    else break; 
  }

  /* remove trailing white space */
  for (i=0;i<fl;++i) if (spaceyname[i]=='\n' || spaceyname[i]==' ') { 
    spaceyname[i]='\0'; 
    break; 
  }
}


