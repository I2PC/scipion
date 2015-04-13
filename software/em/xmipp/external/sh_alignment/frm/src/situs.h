/*********************************************************************
*                           s i t u s . h                            *
**********************************************************************
* Header file for the Situs package    URL: http://situs.scripps.edu *
* (c) Willy Wriggers and Pablo Chacon, 1998-2002                     *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>

#ifndef SITUS_H
#define SITUS_H

/********************************************************************/
/**** Warning: changes in the following affect the performance ******/
/**** of the programs. See documentation for more information. ******/
/********************************************************************/

#define NNMIN 2       /* minimum possible # of codebook vectors     */
#define NNMAX 20      /* maximum possible # of codebook vectors     */
#define NHMAX 200     /* max. pos # of codeb. vectors for helices   */
#define MAXCYCLE 8    /* max # of cycles in cluster analysis        */
#define MAXMON 20     /* max # of monomers for helices              */
#define SMAXS 100000  /* # of neural gas iteration steps            */
#define SMAXH 1000000 /* # of neural gas iteration steps f. helices */      
#define MAXPDB 100000 /* maximum number of lines in pdb file        */
#define BINS 50       /* number of voxel histogram bins             */
#define BARL 70       /* available space for histogram bars         */
#define FLENGTH 1000  /* file name length                           */

/********************************************************************/
/***** Warning: changes below this line are not recommended! ********/
/********************************************************************/

#define MAXDT (NNMAX*MAXCYCLE) 
#define MAXDH (NHMAX*MAXMON) 
#define PI 3.1415926535
#define ROT_CONV 0.017453293 
#define	TRUE 1
#define	FALSE 0
#define SWAPPING(_a,_b,_type) \
{\
  _type _tmp;\
  \
  _tmp = (_a);\
  (_a) = (_b);\
  (_b) = _tmp;\
}

typedef struct
{ char  recdName[7];    /*       1 -  6 */
  int   serial;         /*       7 - 11 */
  char  atomType[3];
  char  atomLoc[3];
  char  altLoc[2];      /*           17 */
  char  resName[5];     /*      18 - 21 */
  char  chainID[2];     /*           22 */
  int   resSeq;         /*      23 - 26 */
  char  iCode[2];       /*           27 */
  float x;              /*      31 - 38 */
  float y;              /*      39 - 46 */
  float z;              /*      47 - 54 */
  float occupancy;      /*      55 - 60 */
  float tempFactor;     /*      61 - 66 */
  int   ftNote;         /*      68 - 70 */
  char  segID[5];       /*      73 - 76 */
  char  element[3];     /*      77 - 78 */
  char  charge[3];      /*      79 - 80 */
  float weight;         /* mass of atom */
} PDB;

typedef struct {
  double score;
  double pos[3];
  double euler[3];
} FIT;    

typedef struct{
  unsigned eu;
  float score;
} SAV;

typedef struct{
  unsigned long ifft;
  unsigned long ireal;
  unsigned ix;
  unsigned iy;
  unsigned iz;
} POS;


typedef double Rmat3 [3][3];
typedef int Iseq3 [3];
typedef double Rseq3 [3];
typedef Rseq3 Rmat3NN [NNMAX];
typedef Rseq3 Rmat3NH [NHMAX];
typedef Rseq3 Rmat3DT [MAXDT];
typedef Rseq3 Rmat3DH [MAXDH];
typedef int IseqNN [NNMAX];
typedef int IseqNH [NHMAX];
typedef double RseqNN [NNMAX];
typedef double RseqNH [NHMAX];
typedef int ImatNNDT [NNMAX][MAXDT];
typedef int ImatNHDH [NHMAX][MAXDH];

typedef struct timeval the_time;

#endif








