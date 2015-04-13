/*............................................................................
 Filename: error.h

 Project: Biomedical Imaging Library

 Author:  Philippe Thevenaz
    Swiss Federal Institute of Technology--Lausanne
    Biomedical Imaging Group
    EPFL/DMT/IOA
    BM-Ecublens
    CH-1015 Lausanne
    Switzerland

 Date:  February 24, 1999

 Purpose: Declaration of macro for error message
............................................................................*/



/****************************************************************************/
/* Defines                 */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
#undef  WRITE_ERROR
#define  WRITE_ERROR(FunctionName, String) \
    { \
        char ErrorMessage[256]; \
        if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
                    __FILE__ ": ERROR---", (long)__LINE__) == EOF) \
            MessageDisplay("\a"); \
        else { \
            MessageDisplay(ErrorMessage); \
            MessageDisplay(" " String "\n"); \
        } \
    }

/*--------------------------------------------------------------------------*/
#undef  WRITE_WARNING
#ifdef  DEBUG
#define  WRITE_WARNING(FunctionName, String) \
    { \
        char ErrorMessage[256]; \
        void *isValidIdentifier; \
        isValidIdentifier = (void *)(*FunctionName); \
        if (((strstr(DEBUG_CONTEXT, "|" #FunctionName "|") != (char *)NULL) \
             || (strstr(DEBUG_CONTEXT, "|*|") != (char *)NULL)) \
            && ((strstr(DEBUG_LEVEL, "|WARNING|") != (char *)NULL) \
                || (strstr(DEBUG_LEVEL, "|*|") != (char *)NULL))) { \
            if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
                        __FILE__ ": WARNING---", (long)__LINE__) == EOF) \
                MessageDisplay("\a"); \
            else { \
                MessageDisplay(ErrorMessage); \
                MessageDisplay(" " String "\n"); \
            } \
        } \
    }
#else
#define  WRITE_WARNING(FunctionName, String) \
    { \
        char ErrorMessage[256]; \
        if (sprintf(ErrorMessage, #FunctionName " at line %ld in " \
                    __FILE__ ": WARNING---", (long)__LINE__) == EOF) \
            MessageDisplay("\a"); \
        else { \
            MessageDisplay(ErrorMessage); \
            MessageDisplay(" " String "\n"); \
        } \
    }
#endif

