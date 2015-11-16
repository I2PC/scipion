#ifndef DEFINES_H
#define DEFINES_H

//#define DEBUG
//#define DEBUG_CLOSEST
//#define DEBUG_TWO_CLOSEST
//#define WINDOWS

// Trigonometric defines.
#ifndef PI
#define PI          3.1415926535
#endif

#define FAILURE     0
#define SUCCESS     1

// Boolean defines.
#define FALSE 	    0
#define TRUE	    1

// Validity defines.
#define VALID       1
#define INVALID     -1

#define DELTA_DIFF  0.05

// Imaginary points.
#define     P_MINUS_1           -1
#define     P_MINUS_2           -2

// Window size
#define MAX_X_COORD	10000
#define MAX_Y_COORD	10000

// File name max length.
#ifndef FILENAME_MAX
#define FILENAME_MAX			200
#endif

#define TYPE		float

// Point coordinate type.
#define POINT_T		TYPE

#endif
