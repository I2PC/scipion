/*............................................................................
	Filename:	debug.h

	Project:	Biomedical Imaging Library

	Author:		Philippe Thevenaz
				Swiss Federal Institute of Technology--Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA
				BM-Ecublens
				CH-1015 Lausanne
				Switzerland

	Date:		February 4, 1999

	Purpose:	Declaration of macros that deal with debug messages
............................................................................*/



/****************************************************************************/
/*	Defines																	*/
/****************************************************************************/

/*==========================================================================*/
#ifdef		DEBUG

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_NULL_POINTER
#define		DEBUG_CHECK_NULL_POINTER(FunctionName, Pointer, Status, String) \
				{ \
					if ((void *)Pointer == (void *)NULL) { \
						Status = ERROR; \
						if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
							|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
							&& ((strstr(DEBUG_LEVEL, "|RANGE_CHECK|") != (char *)NULL) \
							|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
							char	ErrorMessage[256]; \
							void	*isValidIdentifier; \
							isValidIdentifier = (void *)(*FunctionName); \
							if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
								__FILE__ ": NULL POINTER---" #Pointer, (long)__LINE__) == EOF) \
								MessageDisplay("\a"); \
							else { \
								MessageDisplay(ErrorMessage); \
								MessageDisplay("	" String "\n"); \
							} \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_CHAR
#define		DEBUG_CHECK_RANGE_CHAR(FunctionName, CharVariable, Min, Max, Status, String) \
				{ \
					if (((char)CharVariable < (char)Min) \
						|| ((char)Max < (char)CharVariable)) { \
						Status = ERROR; \
						if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
							|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
							&& ((strstr(DEBUG_LEVEL, "|RANGE_CHECK|") != (char *)NULL) \
							|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
							char	ErrorMessage[256]; \
							void	*isValidIdentifier; \
							isValidIdentifier = (void *)(*FunctionName); \
							if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
								__FILE__ ": INVALID CHAR---" #CharVariable " = %c", \
								(long)__LINE__, (char)CharVariable) == EOF) \
								MessageDisplay("\a"); \
							else { \
								MessageDisplay(ErrorMessage); \
								MessageDisplay("	" String "\n"); \
							} \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_DOUBLE
#define		DEBUG_CHECK_RANGE_DOUBLE(FunctionName, DoubleVariable, Min, Max, Status, String) \
				{ \
					if (((double)DoubleVariable < (double)Min) \
						|| ((double)Max < (double)DoubleVariable)) { \
						Status = ERROR; \
						if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
							|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
							&& ((strstr(DEBUG_LEVEL, "|RANGE_CHECK|") != (char *)NULL) \
							|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
							char	ErrorMessage[256]; \
							void	*isValidIdentifier; \
							isValidIdentifier = (void *)(*FunctionName); \
							if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
								__FILE__ ": INVALID DOUBLE---" #DoubleVariable " = %E", \
								(long)__LINE__, (double)DoubleVariable) == EOF) \
								MessageDisplay("\a"); \
							else { \
								MessageDisplay(ErrorMessage); \
								MessageDisplay("	" String "\n"); \
							} \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_FLOAT
#define		DEBUG_CHECK_RANGE_FLOAT(FunctionName, FloatVariable, Min, Max, Status, String) \
				{ \
					if (((float)FloatVariable < (float)Min) \
						|| ((float)Max < (float)FloatVariable)) { \
						Status = ERROR; \
						if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
							|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
							&& ((strstr(DEBUG_LEVEL, "|RANGE_CHECK|") != (char *)NULL) \
							|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
							char	ErrorMessage[256]; \
							void	*isValidIdentifier; \
							isValidIdentifier = (void *)(*FunctionName); \
							if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
								__FILE__ ": INVALID FLOAT---" #FloatVariable " = %E", \
								(long)__LINE__, (float)FloatVariable) == EOF) \
								MessageDisplay("\a"); \
							else { \
								MessageDisplay(ErrorMessage); \
								MessageDisplay("	" String "\n"); \
							} \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_INT
#define		DEBUG_CHECK_RANGE_INT(FunctionName, IntVariable, Min, Max, Status, String) \
				{ \
					if (((int)IntVariable < (int)Min) \
						|| ((int)Max < (int)IntVariable)) { \
						Status = ERROR; \
						if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
							|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
							&& ((strstr(DEBUG_LEVEL, "|RANGE_CHECK|") != (char *)NULL) \
							|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
							char	ErrorMessage[256]; \
							void	*isValidIdentifier; \
							isValidIdentifier = (void *)(*FunctionName); \
							if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
								__FILE__ ": INVALID INT---" #IntVariable " = %d", \
								(long)__LINE__, (int)IntVariable) == EOF) \
								MessageDisplay("\a"); \
							else { \
								MessageDisplay(ErrorMessage); \
								MessageDisplay("	" String "\n"); \
							} \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_LONG
#define		DEBUG_CHECK_RANGE_LONG(FunctionName, LongVariable, Min, Max, Status, String) \
				{ \
					if (((long)LongVariable < (long)Min) \
						|| ((long)Max < (long)LongVariable)) { \
						Status = ERROR; \
						if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
							|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
							&& ((strstr(DEBUG_LEVEL, "|RANGE_CHECK|") != (char *)NULL) \
							|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
							char	ErrorMessage[256]; \
							void	*isValidIdentifier; \
							isValidIdentifier = (void *)(*FunctionName); \
							if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
								__FILE__ ": INVALID LONG---" #LongVariable " = %ld", \
								(long)__LINE__, (long)LongVariable) == EOF) \
								MessageDisplay("\a"); \
							else { \
								MessageDisplay(ErrorMessage); \
								MessageDisplay("	" String "\n"); \
							} \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_SHORT
#define		DEBUG_CHECK_RANGE_SHORT(FunctionName, ShortVariable, Min, Max, Status, String) \
				{ \
					if (((short)ShortVariable < (short)Min) \
						|| ((short)Max < (short)ShortVariable)) { \
						Status = ERROR; \
						if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
							|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
							&& ((strstr(DEBUG_LEVEL, "|RANGE_CHECK|") != (char *)NULL) \
							|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
							char	ErrorMessage[256]; \
							void	*isValidIdentifier; \
							isValidIdentifier = (void *)(*FunctionName); \
							if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
								__FILE__ ": INVALID SHORT---" #ShortVariable " = %hd", \
								(long)__LINE__, (short)ShortVariable) == EOF) \
								MessageDisplay("\a"); \
							else { \
								MessageDisplay(ErrorMessage); \
								MessageDisplay("	" String "\n"); \
							} \
						} \
					} \
				}

/*----------------------------------------------------------------------------*/
#undef		DEBUG_RETURN_ON_ERROR
#define		DEBUG_RETURN_ON_ERROR(FunctionName, Status) \
				{ \
					if (Status == ERROR) \
						return(ERROR); \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_ENTERING
#define		DEBUG_WRITE_ENTERING(FunctionName, String) \
				{ \
					if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
						|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
						&& ((strstr(DEBUG_LEVEL, "|ENTERING_FUNCTION|") != (char *)NULL) \
						|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
						char	ErrorMessage[256]; \
						void	*isValidIdentifier; \
						isValidIdentifier = (void *)(*FunctionName); \
						if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
							__FILE__ ": ENTERING---", (long)__LINE__) == EOF) \
							MessageDisplay("\a"); \
						else { \
							MessageDisplay(ErrorMessage); \
							MessageDisplay("	" String "\n"); \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_INFO
#define		DEBUG_WRITE_INFO(FunctionName, String) \
				{ \
					if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
						|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
						&& ((strstr(DEBUG_LEVEL, "|INFO|") != (char *)NULL) \
						|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
						char	ErrorMessage[256]; \
						void	*isValidIdentifier; \
						isValidIdentifier = (void *)(*FunctionName); \
						if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
							__FILE__ ": INFO---", (long)__LINE__) == EOF) \
							MessageDisplay("\a"); \
						else { \
							MessageDisplay(ErrorMessage); \
							MessageDisplay("	" String "\n"); \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_LEAVING
#define		DEBUG_WRITE_LEAVING(FunctionName, String) \
				{ \
					if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
						|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
						&& ((strstr(DEBUG_LEVEL, "|LEAVING_FUNCTION|") != (char *)NULL) \
						|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
						char	ErrorMessage[256]; \
						void	*isValidIdentifier; \
						isValidIdentifier = (void *)(*FunctionName); \
						if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
							__FILE__ ": LEAVING---", (long)__LINE__) == EOF) \
							MessageDisplay("\a"); \
						else { \
							MessageDisplay(ErrorMessage); \
							MessageDisplay("	" String "\n"); \
						} \
					} \
				}

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_PERSONAL
#define		DEBUG_WRITE_PERSONAL(FunctionName, String) \
				{ \
					if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
						|| (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
						&& ((strstr(DEBUG_LEVEL, "|PERSONAL|") != (char *)NULL) \
						|| (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
						char	ErrorMessage[256]; \
						void	*isValidIdentifier; \
						isValidIdentifier = (void *)(*FunctionName); \
						if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
							__FILE__ ": PERSONAL---", (long)__LINE__) == EOF) \
							MessageDisplay("\a"); \
						else { \
							MessageDisplay(ErrorMessage); \
							MessageDisplay("	" String "\n"); \
						} \
					} \
				}


#else

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_NULL_POINTER
#define		DEBUG_CHECK_NULL_POINTER(FunctionName, Pointer, Status, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_CHAR
#define		DEBUG_CHECK_RANGE_CHAR(FunctionName, CharVariable, Min, Max, Status, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_DOUBLE
#define		DEBUG_CHECK_RANGE_DOUBLE(FunctionName, DoubleVariable, Min, Max, Status, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_FLOAT
#define		DEBUG_CHECK_RANGE_FLOAT(FunctionName, FloatVariable, Min, Max, Status, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_INT
#define		DEBUG_CHECK_RANGE_INT(FunctionName, IntVariable, Min, Max, Status, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_LONG
#define		DEBUG_CHECK_RANGE_LONG(FunctionName, LongVariable, Min, Max, Status, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_CHECK_RANGE_SHORT
#define		DEBUG_CHECK_RANGE_SHORT(FunctionName, ShortVariable, Min, Max, Status, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_RETURN_ON_ERROR
#define		DEBUG_RETURN_ON_ERROR(FunctionName, Status)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_ENTERING
#define		DEBUG_WRITE_ENTERING(FunctionName, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_INFO
#define		DEBUG_WRITE_INFO(FunctionName, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_LEAVING
#define		DEBUG_WRITE_LEAVING(FunctionName, String)

/*--------------------------------------------------------------------------*/
#undef		DEBUG_WRITE_PERSONAL
#define		DEBUG_WRITE_PERSONAL(FunctionName, String)


#endif
