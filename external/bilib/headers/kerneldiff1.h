/*--------------------------------------------------------------------------*/
extern int  BlipDiff1
    (
        long Degree,    /* degree */
        double Argument,   /* input */
        double *Result    /* output */
    );

/*--------------------------------------------------------------------------*/
extern double Blip00Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Blip01Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Blip03Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern int  BsplineDiff1
    (
        long Degree,    /* degree */
        double Argument,   /* input */
        double *Result    /* output */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline00Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline01Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline02Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline03Diff1
    (
        double Argument   /* input */
    );

/** Bspline03Diff1 as a macro */
#define BSPLINE03DIFF1(y,x) \
{\
	double a = fabs(x); \
	if (a < 1.0) \
	{ \
		a *= a * 1.5 - 2.0; \
		y=(x>0.0) ? (a) : (-a); \
	} \
	else if (a < 2.0) { \
		a = 2.0 - a; \
		a *= a * -0.5; \
		y=(x>0.0) ? (a) : (-a); \
	} \
	else \
		y = 0.0; \
}

/*--------------------------------------------------------------------------*/
extern double Bspline04Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline05Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline06Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline07Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline08Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline09Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline10Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Bspline11Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double DodgsonDiff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double German04Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double KeysOptimalDiff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern int  OmomsDiff1
    (
        long Degree,    /* degree */
        double Argument,   /* input */
        double *Result    /* output */
    );

/*--------------------------------------------------------------------------*/
extern double Omoms00Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Omoms01Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Omoms02Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double Omoms03Diff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double PositiveDiff1
    (
        double Argument   /* input */
    );

/*--------------------------------------------------------------------------*/
extern double SincDiff1
    (
        double Argument   /* input */
    );

