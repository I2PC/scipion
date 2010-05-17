/*****************************************************************************
 *  The definition that the user may change are in between @ signs
 *  To modify a definition, permute the two lines that define or undefine things
 ****************************************************************************/


/*--- Defines ----------------------------------------------------------------*/
/*****************************************************************************
 *  Undefine user-modifiable names
 ****************************************************************************/
#undef  CODEWARRIOR
#undef  CC
#undef  GCC
#undef  DEBUG

/*****************************************************************************
 *  Define some constants
 ****************************************************************************/
#undef  FALSE
#define  FALSE    ((int)(0 != 0))
#undef  TRUE
#define  TRUE    (!FALSE)
#undef  ERROR
#define  ERROR    TRUE
#undef  PI
#define  PI     ((double)3.14159265358979323846264338327950288419716939937510)


/*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
------- User-modifiable names start here ---------------------------------------
  ________________________________           */
/*    @    @    @    @     */
/*    @    @    @    @     */
/*    @    @    @    @     */
/*     @ @ @     @ @ @     @ @ @     @ @ @     */
/*      @@@      @@@      @@@      @@@     */
/*    @    @    @    @     */


/*****************************************************************************
 *  Selection of the compiler
 ****************************************************************************/
#define  CODEWARRIOR

/*****************************************************************************
 *  Selection of the debugging mode (defined or undefined)
 ****************************************************************************/
#define  DEBUG
#undef  DEBUG

/*****************************************************************************
 *  Selection of the functions to debug
 ****************************************************************************/
/* Use "||"      for an empty list */
/* Use "|Function1|Function2|"  for specific functions */
/* Use "|*|"      for every function */
#undef  DEBUG_CONTEXT
#define  DEBUG_CONTEXT  "|*|"

/*****************************************************************************
 *  Selection of the debug level
 ****************************************************************************/
/* Use "||" for an empty list */
/* Recognized possibilities are
 "|ENTER_FUNCTION|INFO|LEAVE_FUNCTION|PERSONAL|RANGE_CHECK|WARNING|*|" */
#undef  DEBUG_LEVEL
#define  DEBUG_LEVEL   "|*|"


/*    @    @    @    @     */
/*      @@@      @@@      @@@      @@@     */
/*     @ @ @     @ @ @     @ @ @     @ @ @     */
/*    @    @    @    @     */
/*    @    @    @    @     */
/*    @    @    @    @     */
/*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
------- User-modifiable names end here -----------------------------------------
  ______________________________            */


/*=== Define private switches ================================================*/
#ifdef CODEWARRIOR
#pragma  ANSI_strict   on
#pragma  only_std_keywords on
#pragma  warning_errors  on
#pragma  warn_emptydecl  on
#pragma  warn_extracomma  on
#pragma  warn_implicitconv on
#pragma  warn_illpragma  on
#pragma  warn_possunwant  on
#pragma  warn_unusedarg  on
#pragma  warn_unusedvar  on
#endif

