#include "memory.h"
char*  askMemory(size_t memsize)
{
    char*		ptr = NULL;

    if ( memsize == 0 ) {
        REPORT_ERROR(ERR_MEM_BADREQUEST, "Error in askMemory: Memory allocation size requested is zero!");
        return(NULL);
    }

    if ( ( ptr = (char *) calloc(1,memsize*sizeof(char)) ) == NULL ) {
        std::cerr<<"Memory allocation of %ld bytes failed, memsize= "<< memsize<<std::endl;
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Error in askMemory");
        return(NULL);
    }

    //memset(ptr, 0, memsize);

    return(ptr);
}

int  freeMemory(void* ptr, unsigned long memsize)
{
    if ( ptr == NULL )
        return(0);

    if ( memsize < 1 )
    {
        return(-1);
    }

    free(ptr);
    ptr = NULL;
    return(0);
}

