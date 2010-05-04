#include "memory.h"
char*  askMemory(unsigned long memsize) 
{ 
    char*		ptr = NULL;
    
    if ( memsize == 0 ) {
        REPORT_ERROR(1, "Error in askMemory: Memory allocation size requested is zero!");
        return(NULL);
    }
	
    if ( ( ptr = (char *) malloc(memsize*sizeof(char)) ) == NULL ) { 
        std::cerr<<"Memory allocation of %ld bytes failed, memsize= "<< memsize<<std::endl;
        REPORT_ERROR(2, "Error in askMemory");
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

