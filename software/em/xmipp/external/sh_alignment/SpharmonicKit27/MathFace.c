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


/*************************************************************************/

#include <math.h>
#include <ctype.h>
#include <stdio.h>

#define LOG_10_2 0.3010299956639812  /* log base 10 of 2 */

int MMChopLevel = -15;  /* numbers smaller than 10^ChopLevel are treated as 0 */
int MMPrecision = 16;    /* number of digits after the radix */


int readNumber(FILE *f,
	       double *real, double *imag)
{
  int c, gotImag, foo, scale;
  double negate;

  int nread=fscanf(f, " ");

  gotImag = 0;
  if((c = getc(f)) == '-') {
    negate = -1.0;
    nread=fscanf(f, " ");
  } else {
    negate = 1.0;
    ungetc(c, f);
  }

  /* added by Sean */
  if ((c = getc(f)) == '(')
    nread=fscanf(f," ");
  else
    ungetc(c,f);
  /* end of added by Sean */

  if((c = getc(f)) == 'I') {
    *real = 0.0;
    *imag = negate;
    gotImag = 1;
    nread=fscanf(f, " ");
  } else if((c == '.') || isdigit(c)) {
    ungetc(c, f);
    foo = fscanf(f, " %lf ", real);
    *real *= negate;
    /*
     * scanf only reads scientific notation with small e's, so I have
     *  to catch capital E's
     */
    if(((c = getc(f)) != '*') && (c != 'E')) {
      ungetc(c, f);
    } else {
      nread=fscanf(f, " ");
      if(c == 'E') {
	if((c = getc(f)) != '+') {
	  ungetc(c, f);
	}
	nread=fscanf(f, "%d ", &scale);
	*real *= pow(10.0, (double) scale);	
      } else if((c = getc(f)) == '1') {
	nread=fscanf(f, "0 ^ ");
	if((c = getc(f)) != '+') {
	  ungetc(c, f);
	}
	nread=fscanf(f, "%d ", &scale);
	*real *= pow(10.0, (double) scale);

	/* added by Sean */
	if ((c = getc(f)) == ')')
	  nread=fscanf(f," ");
	else
	  ungetc(c,f);
	/* end of added by Sean */

      } else if(c == 'I') {
	*imag = *real;
	*real = 0.0;
	gotImag = 1;
	nread=fscanf(f, " ");
      } else {
	ungetc(c, f);
      }
    }
  } else {
    ungetc(c, f);
    return -1;
  }

  if((((c = getc(f)) == '+') || (c == '-')) && (gotImag == 0)) {
    if(c == '-') {
      negate = -1.0;
    } else {
      negate = 1.0;
    }
    nread=fscanf(f, " ");
    if((c = getc(f)) == 'I') {
      *imag = negate;
      nread=fscanf(f, " ");
    } else if((c == '.') || isdigit(c)) {
      ungetc(c, f);
      nread=fscanf(f, " %lf ", imag);
      *imag *= negate;
      /*
       * scanf only reads scientific notation with small e's, so I have
       *  to catch capital E's
       */
      if((c = getc(f)) == 'E') {
	nread=fscanf(f, " ");
	if((c = getc(f)) != '+') {
	  ungetc(c, f);
	}
	nread=fscanf(f, "%d * ", &scale);
	*imag *= pow(10.0, (double) scale);
      } else if(c != '*') {
	return -1;
      }
      nread=fscanf(f, " ");
      if((c = getc(f)) != '1') {
	ungetc(c, f);
      } else {
	nread=fscanf(f, "0 ^ ");
	if((c = getc(f)) != '+') {
	  ungetc(c, f);
	}
	nread=fscanf(f, "%d * ", &scale);
	*imag *= pow(10.0, (double) scale);
      }
      if((c = getc(f)) != 'I') {
	ungetc(c, f);
      }
      nread=fscanf(f, " ");
    } else {
      ungetc(c, f);
      return -1;
    }
  } else if((c == '*') && (gotImag == 0)) {
    nread=fscanf(f, " ");
    if((c = getc(f)) == 'I') {
      *imag = *real;
      *real = 0.0;
      nread=fscanf(f, " ");
    } else {
      ungetc(c, f);
    }
  } else {
    ungetc(c, f);
    if(gotImag == 0) {
      *imag = 0.0;
    }
  }

  return 0;
}


/* This may no longer work as I had to remove calls to isnand
   for this to work on trinidad */
void printNumber(FILE *f, double real, double imag)
{
  int needNum;
  double ipart;

  needNum = 1;
  ipart = floor(log10(fabs(real)) + 1.0);
  if(((real != 0.0) && (ipart >= MMChopLevel)) ||
     (imag == 0.0)) {
    fprintf(f, "%.*g*10^%d", MMPrecision,
	    real * pow(10.0, -ipart), (int) ipart);
     needNum = 0;
    }

  if(imag != 0.0) {
    ipart = floor(log10(fabs(imag)) + 1.0);
    if((ipart >= MMChopLevel) && (-4.0 < ipart) && (ipart < 4.0)) {
      if((imag == 1.0) && (needNum)) {
	fprintf(f, "I");
      } else if(imag == 1.0) {
	fprintf(f, "+I");
      } else if(imag == -1.0) {
	fprintf(f, "-I");
      } else if((imag < 0.0) || (needNum)) {
	fprintf(f, "%.*g*I", MMPrecision, imag);
      } else {
	fprintf(f, "+%.*g*I", MMPrecision, imag);
      }
    }

    else if(needNum) {
      /*
       * We didn't print the 0 for the real part, because the
       *  imaginary part was not 0, but the imaginary part fell
       *  below the MMChopLevel and we didn't print it, so we
       *  need to print something
       */
      fprintf(f, "0");
    }
  }
  fflush(f);
}

/* simpler print routine - written to avoid MM conflicts */
void seanprintNumber(FILE *f, double real, double imag)
{

    fprintf(f,"%17.15f + I * %17.15f",real,imag);
}



int readMMRealTable(FILE *f, double *table, int size,
		    int *width, int *height)
/*
 * On input, table points to a region of (size) doubles.  On
 *  return, table is filled in with data, and *width and *height are changed
 *  to indicate actual size.
 *
 * NEW: If the table is 1-D, height is set to 0.
 *
 * RETURNS:      1, if all goes well (table filled and width and height set)
 *               0, if there's a problem (all values undefined)
 */
{
  int i, j, foo, n, c, flat;
  double real, imag;


  n = 0;
  foo = fscanf(f, " { ");
  if((c = getc(f)) != '{') {
    ungetc(c, f);
    flat = 1;
  } else {
    int nread=fscanf(f, " ");
    flat = 0;
  }

  j = 0;
  while(((c = getc(f)) != '}') && (j < size)) {
    ungetc(c, f);
    readNumber(f, &real, &imag);
    if(imag != 0.0) {
      return 0;
    }
    table[n] = real;
    n++;
    j++;
    if((c = getc(f)) != ',') {
      ungetc(c, f);
    } else {
      int nread=fscanf(f, " ");
    }
  }
  if((j == size) && (c != '}')) { return 0; }
  *width = j;
  if(flat == 0) {
    int nread=fscanf(f, " ");
    i = 1;
    while(((c = getc(f)) != '}') && (i*(*width) < size)) {
      foo = fscanf(f, " { ");
      j = 0;
      while(((c = getc(f)) != '}') && (j < *width)) {
	ungetc(c, f);
	readNumber(f, &real, &imag);
	if(imag != 0.0) {
	  return 0;
	}
	table[n] = real;
	n++;
	j++;
	if((c = getc(f)) != ',') {
	  ungetc(c, f);
	} else {
	  nread=fscanf(f, " ");
	}
      }
      if(c != '}') {
	return 0;
      }
      nread=fscanf(f, " ");
      i++;
    }
    if((i*(*width) == size) && (c != '}')) { return 0; }
    *height = i;
  } else {
    *height = 0;
  }

  return 1;
}


void printMMRealTable(FILE *f, double *table, int width, int height)
/*
 * Prints out the table in Mathematica form.  If height is 0, prints it as a
 *  1-D list
 */
{
  int i, j, flat;


  if(height != 0) {
    fprintf(f, "{");
    flat = 0;
  } else {
    height = 1;
    flat = 1;
  }

  for(i = 0; i < height; i++) {
    if(i != 0) {
      fprintf(f, ",\n ");
    }
    fprintf(f, "{");
    for(j = 0; j < width; j++) {
      if(j != 0) {
	putc(',', f);
      }
      printNumber(f, table[width*i+j], 0.0);
    }
    fprintf(f, "}");
    fflush(f);
  }
  if(flat == 0) {
    fprintf(f, "}\n");
  } else {
    putc('\n', f);
  }
}



int readMMComplexTable(FILE *f, double *real, double *imag,
		       int size, int *width, int *height,
		       int *isReal)
/*
 * On input, real and imag point to regions of (size) doubles.  On
 *  return, they are filled in with data, and *width and *height are changed
 *  to indicate actual size.  If none of the values are actually complex,
 *  isReal is set to 1, otherwise it is set to 0.
 *
 * NEW: If the table is 1-D, height is set to 0.
 *
 * RETURNS:      1, if all goes well (real and imag filled and width and height set)
 *               0, if there's a problem (all values undefined)
 */
{
  int i, j, n, c, flat, foo;

  n = 0;
  foo = fscanf(f, " { ");
  if((c = getc(f)) != '{') {
    ungetc(c, f);
    flat = 1;
  } else {
    int nread=fscanf(f, " ");
    flat = 0;
  }

  *isReal = 1;

  j = 0;
  while(((c = getc(f)) != '}') && (j < size)) {
    ungetc(c, f);
    readNumber(f, &(real[n]), &(imag[n]));
    if(imag[n] != 0.0) {
      *isReal = 0;
    }
    n++;
    j++;
    if((c = getc(f)) != ',') {
      ungetc(c, f);
    } else {
      int nread=fscanf(f, " ");
    }
  }

  if((j == size) && (c != '}')) {
    return 0;
  }

  *width = j;
  if(flat == 0) {
    int nread=fscanf(f, " ");
    i = 1;
    while(((c = getc(f)) != '}') && (i*(*width) < size)) {
      foo = fscanf(f, " { ");
      j = 0;
      while(((c = getc(f)) != '}') && (j < *width)) {
	ungetc(c, f);
	readNumber(f, &(real[n]), &(imag[n]));
	if(imag[n] != 0.0) {
	  *isReal = 0;
	}
	n++;
	j++;
	if((c = getc(f)) != ',') {
	  ungetc(c, f);
	} else {
	  nread=fscanf(f, " ");
	}
      }
      if(c != '}') {
	return 0;
      }
      nread=fscanf(f, " ");
      i++;
    }
    if((i*(*width) == size) && (c != '}')) {
      return 0;
    }
    *height = i;
  } else {
    *height = 0;
  }

  return 1;
}


void printMMComplexTable(FILE *f, double *real, double *imag,
			 int width, int height)
/*
 * Prints out the table in Mathematica form.  If height is 0, prints it as a
 *  1-D list
 */
{
  int i, j, flat;

  if(imag == (double *) NULL) {
    printMMRealTable(f, real, width, height);
    return;
  }

  if(height != 0) {
    fprintf(f, "{");
    flat = 0;
  } else {
    height = 1;
    flat = 1;
  }

  for(i = 0; i < height; i++) {
    if(i != 0) {
      fprintf(f, ",\n ");
    }
    fprintf(f, "{");
    for(j = 0; j < width; j++) {
      if(j != 0) {
	putc(',', f);
      }
      printNumber(f, real[width*i+j], imag[width*i+j]);
    }
    fprintf(f, "}");
    fflush(f);
  }

  if(flat == 0) {
    fprintf(f, "}\n");
  } else {
    putc('\n', f);
  }
}

void seanprintMMComplexTable(FILE *f, double *real, double *imag,
			     int width, int height)
/*
 * Prints out the table in Mathematica form.  If height is 0, prints it as a
 *  1-D list
 */
{
  int i, j, flat;

  if(imag == (double *) NULL) {
    printMMRealTable(f, real, width, height);
    return;
  }

  if(height != 0) {
    fprintf(f, "{");
    flat = 0;
  } else {
    height = 1;
    flat = 1;
  }

  for(i = 0; i < height; i++) {
    if(i != 0) {
      fprintf(f, ",\n ");
    }
    fprintf(f, "{");
    for(j = 0; j < width; j++) {
      if(j != 0) {
	putc(',', f);
      }
      seanprintNumber(f, real[width*i+j], imag[width*i+j]);
    }
    fprintf(f, "}");
    fflush(f);
  }

  if(flat == 0) {
    fprintf(f, "}\n");
  } else {
    putc('\n', f);
  }
}

void seanprintMMRealTable(FILE *f, double *table,
			  int width, int height)
/*
 * Prints out the table in Mathematica form.  If height is 0, prints it as a
 *  1-D list
 */
{
  int i, j, flat;


  if(height != 0) {
    fprintf(f, "{");
    flat = 0;
  } else {
    height = 1;
    flat = 1;
  }

  for(i = 0; i < height; i++) {
    if(i != 0) {
      fprintf(f, ",\n ");
    }
    fprintf(f, "{");
    for(j = 0; j < width; j++) {
      if(j != 0) {
	putc(',', f);
      }
      seanprintNumber(f, table[width*i+j], 0.0);
    }
    fprintf(f, "}");
    fflush(f);
  }
  if(flat == 0) {
    fprintf(f, "}\n");
  } else {
    putc('\n', f);
  }
}

