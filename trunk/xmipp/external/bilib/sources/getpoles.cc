/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stdio.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"getpoles.h"
#include	"kernel.h"
#include	"linearalgebra.h"
#include	"messagedisplay.h"
#include	"polynomial.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<float.h>
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
extern int		GetBsplinePoles
				(
					double	RealPoles[],		/* returned array of poles */
					long	Degree,				/* spline degree */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* fill an array with the values of spline poles */
/* the number of returned poles is (Degree / 2L) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetBsplinePoles */

	double	*p;
	double	*Polynomial = (double *)NULL;
	double	*AllPolesReal = (double *)NULL;
	long	i, PoleNumber;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetBsplinePoles, RealPoles, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_DOUBLE(GetBsplinePoles, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(GetBsplinePoles, *Status)
/**/DEBUG_WRITE_ENTERING(GetBsplinePoles,
/**/	"About to fill an array of spline poles")

	switch (Degree) {
		case 0L:
		case 1L:
			*Status = ERROR;
			WRITE_ERROR(GetBsplinePoles, "Invalid spline degree")
			break;
		case 2L:
			RealPoles[0] = sqrt(8.0) - 3.0;
			break;
		case 3L:
			RealPoles[0] = sqrt(3.0) - 2.0;
			break;
		case 4L:
			RealPoles[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			RealPoles[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			break;
		case 5L:
			RealPoles[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			RealPoles[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			break;
		case 6L:
			RealPoles[0] = -0.488294589303044755130118038883789062112279161239377608394;
			RealPoles[1] = -0.081679271076237512597937765737059080653379610398148178525368;
			RealPoles[2] = -0.00141415180832581775108724397655859252786416905534669851652709;
			break;
		case 7L:
			RealPoles[0] = -0.5352804307964381655424037816816460718339231523426924148812;
			RealPoles[1] = -0.122554615192326690515272264359357343605486549427295558490763;
			RealPoles[2] = -0.0091486948096082769285930216516478534156925639545994482648003;
			break;
		case 8L:
			RealPoles[0] = -0.57468690924876543053013930412874542429066157804125211200188;
			RealPoles[1] = -0.163035269297280935240551896860737052234768145508298548489731;
			RealPoles[2] = -0.0236322946948448500234039192963613206126659208546294370651457;
			RealPoles[3] = -0.000153821310641690911739352530184021607629640540700430019629940;
			break;
		case 9L:
			RealPoles[0] = -0.60799738916862577900772082395428976943963471853990829550220;
			RealPoles[1] = -0.201750520193153238796064685055970434680898865757470649318867;
			RealPoles[2] = -0.043222608540481752133321142979429688265852380231497069381435;
			RealPoles[3] = -0.00212130690318081842030489655784862342205485609886239466341517;
			break;
		case 10L:
			RealPoles[0] = -0.63655066396942385875799205491349773313787959010128860432339;
			RealPoles[1] = -0.238182798377573284887456162200161978666543494059728787251924;
			RealPoles[2] = -0.065727033228308551538201803949684252205121694392255863103034;
			RealPoles[3] = -0.0075281946755486906437698340318148831650567567441314086093636;
			RealPoles[4] = -0.0000169827628232746642307274679399688786114400132341362095006930;
			break;
		case 11L:
			RealPoles[0] = -0.66126606890073470691013126292248166961816286716800880802421;
			RealPoles[1] = -0.272180349294785885686295280258287768151235259565335176244192;
			RealPoles[2] = -0.089759599793713309944142676556141542547561966017018544406214;
			RealPoles[3] = -0.0166696273662346560965858360898150837154727205519335156053610;
			RealPoles[4] = -0.00051055753444650205713591952840749392417989252534014106289610;
			break;
		default:
			AllocateVector(&Polynomial, Degree + 1L, Status);
			if (*Status == ERROR) {
				FreeVector(&Polynomial);
/**/			DEBUG_WRITE_LEAVING(GetBsplinePoles, "Done")
				return(*Status);
			}
			p = Polynomial;
			for (i = -Degree / 2L; (i <= (Degree / 2L)); i++) {
				*Status = Bspline(Degree, (double)i, p++);
				if (*Status == ERROR) {
					FreeVector(&Polynomial);
/**/				DEBUG_WRITE_LEAVING(GetBsplinePoles, "Done")
					return(*Status);
				}
			}
			AllocateVector(&AllPolesReal, Degree, Status);
			if (*Status == ERROR) {
				FreeVector(&Polynomial);
/**/			DEBUG_WRITE_LEAVING(GetBsplinePoles, "Done")
				return(*Status);
			}
			PolynomialRealRoots(Polynomial, Degree, AllPolesReal, &PoleNumber,
				Tolerance, Status);
			if (*Status == ERROR) {
				FreeVector(&AllPolesReal);
				FreeVector(&Polynomial);
/**/			DEBUG_WRITE_LEAVING(GetBsplinePoles, "Done")
				return(*Status);
			}
			*Status = FreeVector(&Polynomial);
			if (*Status == ERROR) {
				FreeVector(&AllPolesReal);
/**/			DEBUG_WRITE_LEAVING(GetBsplinePoles, "Done")
				return(*Status);
			}
			for (i = -Degree / 2L; (i <= (Degree / 2L)); i++) {
				if (*AllPolesReal < 0.0) {
					*RealPoles++ = *AllPolesReal;
				}
				AllPolesReal++;
			}
			*Status = FreeVector(&AllPolesReal);
			break;
	}
/**/DEBUG_WRITE_LEAVING(GetBsplinePoles, "Done")
	return(*Status);
} /* end GetBsplinePoles */

/*--------------------------------------------------------------------------*/
extern int		GetOmomsPoles
				(
					double	RealPoles[],		/* returned array of poles */
					long	Degree,				/* oMoms degree */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* fill an array with the values of oMoms poles */
/* the number of returned poles is (Degree / 2L) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetOmomsPoles */

	double	*p;
	double	*Polynomial = (double *)NULL;
	double	*AllPolesReal = (double *)NULL;
	long	i, PoleNumber;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetOmomsPoles, RealPoles, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_DOUBLE(GetOmomsPoles, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(GetOmomsPoles, *Status)
/**/DEBUG_WRITE_ENTERING(GetOmomsPoles,
/**/	"About to fill an array of oMoms poles")

	switch (Degree) {
		case 0L:
		case 1L:
			*Status = ERROR;
			WRITE_ERROR(GetOmomsPoles, "Invalid oMoms degree")
			break;
		case 2L:
			RealPoles[0] = sqrt(1560.0 / 289.0) - 43.0 / 17.0;
			break;
		case 3L:
			RealPoles[0] = sqrt(105.0 / 64.0) - 13.0 / 8.0;
			break;
		case 4L:
			RealPoles[0] = sqrt(68880840.0 / 552049.0 - sqrt(4666909808998080.0
				/ 304758098401.0)) + sqrt(28511280.0 / 552049.0) - 6397.0 / 743.0;
			RealPoles[1] = sqrt(68880840.0 / 552049.0 + sqrt(4666909808998080.0
				/ 304758098401.0)) - sqrt(28511280.0 / 552049.0) - 6397.0 / 743.0;
			break;
		case 5L:
			RealPoles[0] = sqrt(285420.0 / 11449.0 - sqrt(77202800640.0 / 131079601.0))
				+ sqrt(96165.0 / 11449.0) - 448.0 / 107.0;
			RealPoles[1] = sqrt(285420.0 / 11449.0 + sqrt(77202800640.0 / 131079601.0))
				- sqrt(96165.0 / 11449.0) - 448.0 / 107.0;
			break;
		case 6L:
			RealPoles[0] = -0.52667681509090929598702625804052518367299232611182543044727;
			RealPoles[1] = -0.113602213794490707459277514309782377945743188648838906585927;
			RealPoles[2] = -0.0062184195886762433029497449362169006376506731708848957234005;
			break;
		case 7L:
			RealPoles[0] = -0.56853761800229298164787483927617716585997958614632832393829;
			RealPoles[1] = -0.155700774677357760841565001104679442932629669048969307010564;
			RealPoles[2] = -0.0197684253838613956123727080799270645154754590429447748430377;
			break;
		case 8L:
			RealPoles[0] = -0.60339620680495477357624485934607992088475202356468905280774;
			RealPoles[1] = -0.195760866811992318305809476612370514400584120988635886262742;
			RealPoles[2] = -0.039301266108611678563535490932330519287035489696802783105375;
			RealPoles[3] = -0.00131723647225732476410835675436309476443645309809984514964350;
			break;
		case 9L:
			RealPoles[0] = -0.63302645078114260206493405039483007841576319677177837016686;
			RealPoles[1] = -0.233245076299699204682329507753694890665166631856015048882768;
			RealPoles[2] = -0.062016051554856745489712385333682728202357092758651943221208;
			RealPoles[3] = -0.0061077992746160882560927708014738068983721892883552093101361;
			break;
		case 10L:
			RealPoles[0] = -0.65850232025797824268236835412538454392471957294049823149805;
			RealPoles[1] = -0.268073111692664696090244736364513359311223180893644790376100;
			RealPoles[2] = -0.086377390153589325877599225694789952508471024430869300205191;
			RealPoles[3] = -0.0148944273821827884756957892122847609801573110427489810826383;
			RealPoles[4] = -0.000287827212515962754748623276177592144223736764032594845958710;
			break;
		case 11L:
			RealPoles[0] = -0.68065728168114552507405373070786114528043918248200630786792;
			RealPoles[1] = -0.300332056088588604876865222180030612313181015580081525329237;
			RealPoles[2] = -0.111312514706830149369412472075399334332922094176016717237651;
			RealPoles[3] = -0.0269288910734721128228921792213836378705084951700165010332794;
			RealPoles[4] = -0.00197961213746736842587514349996322223537646370088176083984645;
			break;
		default:
			AllocateVector(&Polynomial, Degree + 1L, Status);
			if (*Status == ERROR) {
				FreeVector(&Polynomial);
/**/			DEBUG_WRITE_LEAVING(GetOmomsPoles, "Done")
				return(*Status);
			}
			p = Polynomial;
			for (i = -Degree / 2L; (i <= (Degree / 2L)); i++) {
				*Status = Omoms(Degree, (double)i, p++);
				if (*Status == ERROR) {
					FreeVector(&Polynomial);
/**/				DEBUG_WRITE_LEAVING(GetOmomsPoles, "Done")
					return(*Status);
				}
			}
			AllocateVector(&AllPolesReal, Degree, Status);
			if (*Status == ERROR) {
				FreeVector(&Polynomial);
/**/			DEBUG_WRITE_LEAVING(GetOmomsPoles, "Done")
				return(*Status);
			}
			PolynomialRealRoots(Polynomial, Degree, AllPolesReal, &PoleNumber,
				Tolerance, Status);
			if (*Status == ERROR) {
				FreeVector(&AllPolesReal);
				FreeVector(&Polynomial);
/**/			DEBUG_WRITE_LEAVING(GetOmomsPoles, "Done")
				return(*Status);
			}
			*Status = FreeVector(&Polynomial);
			if (*Status == ERROR) {
				FreeVector(&AllPolesReal);
/**/			DEBUG_WRITE_LEAVING(GetOmomsPoles, "Done")
				return(*Status);
			}
			for (i = -Degree / 2L; (i <= (Degree / 2L)); i++) {
				if (*AllPolesReal < 0.0) {
					*RealPoles++ = *AllPolesReal;
				}
				AllPolesReal++;
			}
			*Status = FreeVector(&AllPolesReal);
			break;
	}
/**/DEBUG_WRITE_LEAVING(GetOmomsPoles, "Done")
	return(*Status);
} /* end GetOmomsPoles */

