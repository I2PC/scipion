/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<stdio.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"
#include	"ttraceline.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"../headers/convert.h"
#include	"messagedisplay.h"
#include	"traceline.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<float.h>
#include	<math.h>
#include	<string.h>
#endif

/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
/* None */

/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/* None */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		FirstIndex3DLine6Connected
				(
					struct TTraceLine
							*LineData			/* line description */
				)

/* get precomputed elements for tracing a discrete 3D 6-connected line */
/* must be called prior to calling NextIndex3DLine6Connected() */
/* k[] returns the first coordinates of the discrete line */
/* Quadrant returns a code for the general orientation (sign) of q[] */
/* this code is as follows: */
/* -1 -1 -1 -1 -1 -1 -1 -1 -1 */
/* -1 -1 -1  0  0  0  1  1  1 */
/* -1  0  1 -1  0  1 -1  0  1 */
/* == == == == == == == == == */
/*  0  9 18  3 12 21  6 15 24 */
/*                            */
/*  0  0  0  0     0  0  0  0 */
/* -1 -1 -1  0     0  1  1  1 */
/* -1  0  1 -1     1 -1  0  1 */
/* == == == ==    == == == == */
/*  1 10 19  4    22  7 16 25 */
/*                            */
/*  1  1  1  1  1  1  1  1  1 */
/* -1 -1 -1  0  0  0  1  1  1 */
/* -1  0  1 -1  0  1 -1  0  1 */
/* == == == == == == == == == */
/*  2 11 20  5 14 23  8 17 26 */
/* ModificationCode returns the direction along which a step will be taken next */
/* the code is as follows: */
/* -3 <-> z--; -2 <-> y--; -1 <-> x--; 1 <-> x++; 2 <-> y++; 3 <-> z++ */
/* (P + Entry * q) belongs to the bounding box of the voxel surrounding k */
/* (P + Exit * q) belongs to the bounding box of the voxel surrounding k */
/* (Entry < Exit) is satisfied */
/* LineData->P[] is a 3D coordinate that belongs to the line */
/* this coordinate is modified after calling FirstIndex3DLine6Connected */
/* q[] is a 3D unit vector that gives the direction in which to draw the line */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin FirstIndex3DLine6Connected */

	double	t0, t1, t2;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FirstIndex3DLine6Connected, LineData, Status,
/**/	"Missing LineData")
/**/DEBUG_RETURN_ON_ERROR(FirstIndex3DLine6Connected, Status)
/**/DEBUG_WRITE_ENTERING(FirstIndex3DLine6Connected,
/**/	"About to precompute elements for drawing a 3D discrete line")

	LineData->k[0] = ConvertDoubleToLong(LineData->P[0]);
	LineData->k[1] = ConvertDoubleToLong(LineData->P[1]);
	LineData->k[2] = ConvertDoubleToLong(LineData->P[2]);
	if (LineData->q[0] < 0.0) {
		LineData->Quadrant = 0;
	}
	else {
		if (LineData->q[0] == 0.0) {
			LineData->Quadrant = 1;
		}
		else {
			LineData->Quadrant = 2;
		}
	}
	if (LineData->q[1] == 0.0) {
		LineData->Quadrant += 3;
	}
	else if (0.0 < LineData->q[1]) {
		LineData->Quadrant += 6;
	}
	if (LineData->q[2] == 0.0) {
		LineData->Quadrant += 9;
	}
	else if (0.0 < LineData->q[2]) {
		LineData->Quadrant += 18;
	}
	if (LineData->Quadrant == 13) {
		Status = ERROR;
		WRITE_ERROR(FirstIndex3DLine6Connected,
			"Invalid vector q (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(FirstIndex3DLine6Connected, "Done")
		return(Status);
	}
	switch (LineData->Quadrant) {
		case 0:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] += 0.5;
			LineData->P[1] += 0.5;
			LineData->P[2] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 1:
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t1 < t2) ? (t2) : (t1);
			LineData->P[1] += 0.5;
			LineData->P[2] += 0.5;
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 2:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] -= 0.5;
			LineData->P[1] += 0.5;
			LineData->P[2] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 3:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t2) ? (t2) : (t0);
			LineData->P[0] += 0.5;
			LineData->P[2] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 4:
			LineData->Entry = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->P[2] += 0.5;
			LineData->Exit = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			LineData->ModificationCode = -3;
			break;
		case 5:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t2) ? (t2) : (t0);
			LineData->P[0] -= 0.5;
			LineData->P[2] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 6:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] += 0.5;
			LineData->P[1] -= 0.5;
			LineData->P[2] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 7:
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t1 < t2) ? (t2) : (t1);
			LineData->P[1] -= 0.5;
			LineData->P[2] += 0.5;
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 8:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] + 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] -= 0.5;
			LineData->P[1] -= 0.5;
			LineData->P[2] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 9:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			LineData->Entry = (t0 < t1) ? (t1) : (t0);
			LineData->P[0] += 0.5;
			LineData->P[1] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			break;
		case 10:
			LineData->Entry = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			LineData->P[1] += 0.5;
			LineData->Exit = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			LineData->ModificationCode = -2;
			break;
		case 11:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			LineData->Entry = (t0 < t1) ? (t1) : (t0);
			LineData->P[0] -= 0.5;
			LineData->P[1] += 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			break;
		case 12:
			LineData->Entry = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			LineData->P[0] += 0.5;
			LineData->Exit = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			LineData->ModificationCode = -1;
			break;
		case 14:
			LineData->Entry = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			LineData->P[0] -= 0.5;
			LineData->Exit = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			LineData->ModificationCode = 1;
			break;
		case 15:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			LineData->Entry = (t0 < t1) ? (t1) : (t0);
			LineData->P[0] += 0.5;
			LineData->P[1] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			break;
		case 16:
			LineData->Entry = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			LineData->P[1] -= 0.5;
			LineData->Exit = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			LineData->ModificationCode = 2;
			break;
		case 17:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			LineData->Entry = (t0 < t1) ? (t1) : (t0);
			LineData->P[0] -= 0.5;
			LineData->P[1] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			break;
		case 18:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] += 0.5;
			LineData->P[1] += 0.5;
			LineData->P[2] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		case 19:
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t1 < t2) ? (t2) : (t1);
			LineData->P[1] += 0.5;
			LineData->P[2] -= 0.5;
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 20:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] + 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] -= 0.5;
			LineData->P[1] += 0.5;
			LineData->P[2] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		case 21:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t2) ? (t2) : (t0);
			LineData->P[0] += 0.5;
			LineData->P[2] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 22:
			LineData->Entry = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->P[2] -= 0.5;
			LineData->Exit = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			LineData->ModificationCode = 3;
			break;
		case 23:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t2) ? (t2) : (t0);
			LineData->P[0] -= 0.5;
			LineData->P[2] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 24:
			t0 = ((double)LineData->k[0] + 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] += 0.5;
			LineData->P[1] -= 0.5;
			LineData->P[2] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		case 25:
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t1 < t2) ? (t2) : (t1);
			LineData->P[1] -= 0.5;
			LineData->P[2] -= 0.5;
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 26:
			t0 = ((double)LineData->k[0] - 0.5 - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - 0.5 - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - 0.5 - LineData->P[2]) / LineData->q[2];
			LineData->Entry = (t0 < t1) ? ((t1 < t2) ? (t2) : (t1)) : ((t0 < t2) ? (t2) : (t0));
			LineData->P[0] -= 0.5;
			LineData->P[1] -= 0.5;
			LineData->P[2] -= 0.5;
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(FirstIndex3DLine6Connected,
				"Invalid Quadrant (shoud belong to [0..26] - {13})")
/**/		DEBUG_WRITE_LEAVING(FirstIndex3DLine6Connected, "Done")
			return(Status);
	}
/**/DEBUG_WRITE_LEAVING(FirstIndex3DLine6Connected, "Done")
	return(Status);
} /* end FirstIndex3DLine6Connected */

/*--------------------------------------------------------------------------*/
extern int		NextIndex3DLine6Connected
				(
					struct TTraceLine
							*LineData			/* line description */
				)

/* updates the 3D index k[] to run along a discrete 6-connected 3D line */
/* see other usage notes at FirstIndex3DLine6Connected() */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin NextIndex3DLine6Connected */

	static double
			t0, t1, t2;
	static int
			Status;

/**/DEBUG_CHECK_NULL_POINTER(NextIndex3DLine6Connected, LineData, Status,
/**/	"Missing LineData")
/**/DEBUG_RETURN_ON_ERROR(NextIndex3DLine6Connected, Status)
/**/DEBUG_WRITE_ENTERING(NextIndex3DLine6Connected,
/**/	"About to precompute elements for drawing a 3D discrete line")

	Status = !ERROR;
	LineData->Entry = LineData->Exit;
	if (LineData->ModificationCode < 0) {
		switch (LineData->ModificationCode) {
			case -3:
				(LineData->k[2])--;
				break;
			case -2:
				(LineData->k[1])--;
				break;
			case -1:
				(LineData->k[0])--;
				break;
			default:
				Status = ERROR;
				WRITE_ERROR(NextIndex3DLine6Connected,
					"Invalid ModificationCode (shoud belong to [-3..3] - {0})")
/**/			DEBUG_WRITE_LEAVING(NextIndex3DLine6Connected, "Done")
				return(Status);
		}
	}
	else {
		switch (LineData->ModificationCode) {
			case 1:
				(LineData->k[0])++;
				break;
			case 2:
				(LineData->k[1])++;
				break;
			case 3:
				(LineData->k[2])++;
				break;
			default:
				Status = ERROR;
				WRITE_ERROR(NextIndex3DLine6Connected,
					"Invalid ModificationCode (shoud belong to [-3..3] - {0})")
/**/			DEBUG_WRITE_LEAVING(NextIndex3DLine6Connected, "Done")
				return(Status);
		}
	}
	switch (LineData->Quadrant) {
		case 0:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 1:
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 2:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 3:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 4:
			LineData->Exit = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			break;
		case 5:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 6:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 7:
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = -3;
			}
			break;
		case 8:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = -3;
				}
			}
			break;
		case 9:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			break;
		case 10:
			LineData->Exit = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			break;
		case 11:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			break;
		case 12:
			LineData->Exit = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			break;
		case 14:
			LineData->Exit = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			break;
		case 15:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			break;
		case 16:
			LineData->Exit = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			break;
		case 17:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			if (t0 < t1) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			break;
		case 18:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		case 19:
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = -2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 20:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = -2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		case 21:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = -1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 22:
			LineData->Exit = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			break;
		case 23:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t2) {
				LineData->Exit = t0;
				LineData->ModificationCode = 1;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 24:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = -1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		case 25:
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t1 < t2) {
				LineData->Exit = t1;
				LineData->ModificationCode = 2;
			}
			else {
				LineData->Exit = t2;
				LineData->ModificationCode = 3;
			}
			break;
		case 26:
			t0 = ((double)LineData->k[0] - LineData->P[0]) / LineData->q[0];
			t1 = ((double)LineData->k[1] - LineData->P[1]) / LineData->q[1];
			t2 = ((double)LineData->k[2] - LineData->P[2]) / LineData->q[2];
			if (t0 < t1) {
				if (t0 < t2) {
					LineData->Exit = t0;
					LineData->ModificationCode = 1;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			else {
				if (t1 < t2) {
					LineData->Exit = t1;
					LineData->ModificationCode = 2;
				}
				else {
					LineData->Exit = t2;
					LineData->ModificationCode = 3;
				}
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(NextIndex3DLine6Connected,
				"Invalid Quadrant (shoud belong to [0..26] - {13})")
/**/		DEBUG_WRITE_LEAVING(NextIndex3DLine6Connected, "Done")
			return(Status);
	}
/**/DEBUG_WRITE_LEAVING(NextIndex3DLine6Connected, "Done")
	return(Status);
} /* end NextIndex3DLine6Connected */

