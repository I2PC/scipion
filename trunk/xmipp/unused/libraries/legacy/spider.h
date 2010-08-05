/**************************************************************************/
/*IMAGE.DOC  2/23/83  FORMAT OF 2D AND 3D FILES IN THE SPIDER SYSTEM.     */
/*(LAST UPDATE 3/30/89 AL)                                                */
/*                                                                        */
/*EACH IMAGE IS CONSIDERED AS A SPECIAL CASE OF A VOLUME CONSISTING OF    */
/*A SINGLE SLICE (OR PLANE).                                              */
/*                                                                        */
/*A VOLUME OF NSAM X NROW X NSLICE VOXELS IS STORED IN THE                */
/*FOLLOWING WAY:                                                          */
/*                                                                        */
/* UNFORMATTED, DIRECT ACCESS FILE CONTAINING A                           */
/* TOTAL OF NSLICE X NROW +LABREC RECORDS, EACH RECORD                    */
/* CONTAINING NSAM 4-BYTE WORDS TO BE INTERPRETED                         */
/* AS FLOATING POINT NUMBERS.                                             */
/*                                                                        */
/*THERE ARE CURRENTLY TWO DIFFERENT FILE FORMATS IN USE.  THE OLD         */
/*SHORT LABEL FORMAT HAS A SINGLE FIRST RECORD FOR THE SPIDER LABEL.      */
/*THE NEW LONG LABEL FORMAT HAS AS MANY RECORDS IN THE LABEL AS ARE       */
/*NECESSARY TO HOLD 256 FLOATING POINT VARIABLES.  THE NEW FORMAT         */
/*ALLOWS MUCH GREATER INFORMATION STORAGE IN THE LABEL AND AVOIDS         */
/*TRUNCATION OF TITLE WHEN USING SMALL IMAGES (NSAM < 20)                 */
/*                                                                        */
/*THE SEQUENCE IN WHICH THE INFORMATION IS STORED IN THE                  */
/*FILE IS AS FOLLOWS:                                                     */
/*                                                                        */
/* RECORDS NO. 1 THRU LABREC -- SPIDER LABEL                              */
/*          FOR OLD FORMAT FILES LABREC = 1                               */
/*          FOR NEW FORMAT FILES LABREC = CIELING OF (256/NSAM)           */
/*                                                                        */
/* RECORDS NO. (LABREC+1) THRU NO. (NROW+LABREC) -- SLICE NO.1            */
/* RECORDS NO. (NROW+LABREC+1) THRU (2*NROW+LABREC) -- SLICE NO.2         */
/* .                                                                      */
/* .                                                                      */
/* .                                                                      */
/* RECORDS NO. ((ISL-1)*NROW+LABREC+1) THRU (ISL*NROW+LABREC) --          */
/*                SLICE NO. ISL                                           */
/* .                                                                      */
/* .                                                                      */
/* RECORDS NO. ((NSLICE-1)*NROW+LABREC+1) THRU (NSLICE*NROW+LABREC)       */
/*                 SLICE NO. NSLICE                                       */
/*                                                                        */
/*                                                                        */
/*LAYOUT OF THE SPIDER LABEL IS AS FOLLOWS:                               */
/*                                                                        */
/**************************************************************************/
#ifndef H_SPIDER
#define H_SPIDER
typedef struct CABECERO
{
    /*1 */      float fNslice;/* = NUMBER OF SLICES (PLANES) IN VOLUME          */
    /*  */                  /*   (=1 FOR AN IMAGE)  FOR NEW LONG LABEL        */
    /*  */                /*   FORMAT THE VALUE OF NSLICE STORED IN         */
    /*  */                /*   THE FILE IS NEGATIVE.                        */
    /*2 */
    float fNrow;   /* = NUMBER OF ROWS PER SLICE                     */
    /*3 */
    float fNrec;   /* = TOTAL NUMBER OF RECORDS (SEE NOTE #3).       */
    /*4 */
    float fNlabel; /* = AUXILIARY NUMBER TO COMPUTE TOTAL NUMBER     */
    /*  */                /*   OF RECS.                                     */
    /*5 */
    float fIform;  /* = FILE TYPE SPECIFIER.                         */
    /*  */                  /* = +3 FOR A 3-D FILE                          */
    /*  */              /* = +1   FOR A 2-D IMAGE                         */
    /*  */        /* = -1 FOR A 2-D FOURIER TRANSFORM             */
    /*  */        /* = -3 FOR A 3-D FOURIER TRANSFORM             */
    /*  */            /* = -5 FOR A NEW 2-D FOURIER TRANSFORM         */
    /*  */            /* = -7 FOR A NEW 3-D FOURIER TRANSFORM         */
    /*  */            /* = +8 FOR A 2-D EIGHT BIT IMAGE FILE          */
    /*  */                    /* = 11   FOR A 2-D EIGHT BIT COLOR IMAGE         */
    /*  */              /*        FILE                                    */
    /*6 */
    float fImami;  /* = MAXIMUM/MINIMUM FLAG. IS SET AT 0 WHEN THE   */
    /*  */              /*  FILE IS CREATED, AND AT 1 WHEN THE            */
    /*  */            /*  MAXIMUM AND MINIMUM HAVE BEEN COMPUTED,       */
    /*  */            /*  AND HAVE BEEN STORED INTO THIS LABEL          */
    /*  */              /*  RECORD (SEE FOLLOWING WORDS)                  */
    /*7 */
    float fFmax;   /* = MAXIMUM VALUE                                */
    /*8 */
    float fFmin;   /* = MINIMUM VALUE                                */
    /*9 */
    float fAv;     /* = AVERAGE VALUE                                */
    /*10*/
    float fSig;    /* = STANDARD DEVIATION. A VALUE OF -1.           */
    /*  */                /*  INDICATES THAT SIG HAS NOT BEEN COMPUTED      */
    /*  */                /*  PREVIOUSLY.                                   */
    /*11*/
    float fIhist;  /* = FLAG INDICATING IF THE HISTOGRAM HAS BE      */
    /*  */                /*  COMPUTED. NOT USED IN 3D FILES!               */
    /*  */                    /*                                                */
    /*12*/
    float fNsam;   /* = NUMBER OF PIXELS PER LINE                    */
    /*13*/
    float fLabrec; /* = NUMBER OF LABEL RECORDS IN FILE HEADER       */
    /*14*/
    float fIangle; /* = FLAG THAT TILT ANGLES HAVE BEEN FILLED       */
    /*15*/
    float fPhi;    /* = TILT ANGLE                                   */
    /*16*/
    float fTheta;  /* = TILT ANGLE                                   */
    /*17*/
    float fGamma;  /* = PSI  = TILT ANGLE                            */
    /*18*/
    float fXoff;   /* = X TRANSLATION                                */
    /*19*/
    float fYoff;   /* = Y TRANSLATION                                */
    /*20*/
    float fZoff;   /* = Z TRANSLATION                                */
    /*21*/
    float fScale;  /* = SCALE                                        */
    /*22*/
    float fLabbyt; /* = TOTAL NUMBER OF BYTES IN LABEL               */
    /*23*/
    float fLenbyt; /* = RECORD LENGTH IN BYTES                       */
    char  fNada[24];/*this is a spider incongruence*/
    /*30*/
    float fFlag;   /*  THAT ANGLES ARE SET. 1 = ONE ADDITIONAL       */
    /*  */                    /*  ROTATION IS PRESENT, 2 = ADDITIONAL ROTATION  */
    /*  */                    /*  THAT PRECEEDS THE ROTATION THAT WAS STORED IN */
    /*  */                    /*  15 FOR DETAILS SEE MANUAL CHAPTER VOCEUL.MAN  */
    /*31*/
    float fPhi1;   /*                                                */
    /*32*/
    float fTheta1; /*                                                */
    /*33*/
    float fPsi1;   /*                                                */
    /*34*/
    float fPhi2;   /*                                                */
    /*35*/
    float fTheta2; /*                                                */
    /*36*/
    float fPsi2;   /*    */
    /*char  fNada2[700];*/

    double fGeo_matrix[3][3];  /* 8x9 = 72 bytes: Geometric info     */
    float fAngle1;             /*                 angle info         */

    float fr1;
    float fr2; /** lift up cosine mask parameters **/

    /** Fraga 23/05/97  For Radon transforms **/
    float RTflag; /** 1=RT, 2=FFT(RT) **/
    float Astart;
    float Aend;
    float Ainc;
    float Rsigma;  /** 4*7 = 28 bytes **/
    float Tstart;
    float Tend;
    float Tinc;   /** 4*3 = 12, 12+28 = 40B **/

    char  fNada2[584]; /* empty     700-76-40=624-40= 584 bytes */

    /*212-214*/
    char szIDat[12];/*= LOGICAL * 1 ARRAY DIMENSIONED 10,             */
    /*  */                /*  CONTAINING THE DATE OF CREATION               */
    /*  */                /*  (10 CHARS)                                    */
    /*215-216*/
    char szITim[8]; /*= LOGICAL * 1 ARRAY DIMENSIONED 8,              */
    /*  */                /* CONTAINING THE TIME OF CREATION                */
    /*  */                /* (8 CHARS)                                      */
    /*217-256*/
    char szITit[160];/*= LOGICAL * 1 ARRAY DIMENSIONED 160,            */
    /* CONTAINING TITLE (160 CHARS)                   */
}
CABECERO;
/**************************************************************************/
/*NOTE#1 : IN SHORT LABEL FORMAT FILES ALL INTEGERS OR LOGICAL*2 ARRAYS   */
/*         ARE RETRIEVED FROM THE FLOATING POINT BUFFER ARRAY CONTAINING  */
/*  THE LABEL RECORD(S) BY ARITHMETIC ASSIGNMENTS; I.E., NSLICE =         */
/*  BUF(1), ETC.                                                          */
/*                                                                        */
/*NOTE#2 : IN LONG LABEL FORMAT FILES ALL INTEGERS ARE RETRIEVED FROM THE */
/*         FLOATING POINT BUFFER ARRAY BY ARITHMETIC ASSIGNMENT BUT ALL   */
/*   LOGICAL ARRAYS ARE RETRIEVED FROM THE FLOATING POINT BUFFER        */
/*  ARRAY CONTAINING THE LABEL RECORD(S) BY EQUIVALENCE ASSIGNMENTS.      */
/*  THUS LOGICAL ARRAYS ARE STORED IN THE LABEL WITHOUT ANY               */
/*  CONVERSION.                                                           */
/*                                                                        */
/*NOTE#3 : NREC IS THE MAXIMUM NUMBER OF RECORDS USED.  FOR A 3-D FILE,   */
/*  SPACE IS RESERVED FOR THE LABEL RECORD(S), FOR THE 3-D DATA,          */
/*  AND FOR AN EXTRA FOURIER PLANE.  FOR A 2-D FILE, SPACE IS             */
/*  RESERVED FOR THE LABEL RECORD(S), THE IMAGE DATA, AN EXTRA            */
/*  FOURIER LINE, AND FOR HISTOGRAM INFORMATION.  HISTOGRAM AND           */
/*         EXTRA FOURIER RECORDS ARE NOT OCCUPIED UNLESS HISTOGRAM AND    */
/*  FOURIER TRANSFORM HAVE BEEN COMPUTED.                                 */
/**************************************************************************/


typedef struct FICHERO
{
    char szNombre[16];    /*File Name*/
    char szNombre2[16];    /*File Name*/
    char szNombre3[16];    /*File Name*/
    /*long lTamano;    */     /*File size*/
}
FICHERO;



/******************************************************************/
/******************************************************************/

typedef struct
{

    float fImami;
    float fFmax;
    float fFmin;
    float fAv;
    float fSig;
    double fGeo_matrix[3][3];
    float fAngle1;

}
GEO_INFO;
#endif
