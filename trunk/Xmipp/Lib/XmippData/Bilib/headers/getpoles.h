/**@name Poles of B-spline and Omom */
//@{
/*--------------------------------------------------------------------------*/
/** Get poles of B-spline.
    Fill an array with the values of spline poles.
    The number of returned poles is (Degree / 2L).
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		GetBsplinePoles
				(
					double	RealPoles[],		/* returned array of poles */
					long	Degree,				/* spline degree */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Get poles of Omom.
    Fill an array with the values of oMoms poles.
    The number of returned poles is (Degree / 2L).
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		GetOmomsPoles
				(
					double	RealPoles[],		/* returned array of poles */
					long	Degree,				/* oMoms degree */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				);
//@}