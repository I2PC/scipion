#include "image.h"

template<>
void Image< std::complex< double > >::castPage2T(char * page, 
                                                 std::complex<double> * ptrDest, 
                                                 size_t pageSize)
{

    switch (datatype)
    {
    case ComplexShort:
    {
        std::complex<short> * ptr = (std::complex<short> *) page;
        for(int i=0; i<pageSize;i++) ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
        break;
    }
    case ComplexInt:
    {
        std::complex<int> * ptr = (std::complex<int> *) page;
        for(int i=0; i<pageSize;i++) ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
        break;
    }
    case ComplexFloat:
    {
        std::complex<float> * ptr = (std::complex<float> *) page;
        for(int i=0; i<pageSize;i++) ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
        break;
    }
    case ComplexDouble:
    {
        memcpy(ptrDest, page, pageSize*sizeof(std::complex<double>));
        break;
    }
    default:
    {
        std::cerr<<"Datatype= "<<datatype<<std::endl;
        REPORT_ERROR(16," ERROR: cannot cast datatype to std::complex<double>");
        break;
    }
    }    

}

template<>
void Image< std::complex< double > >::castPage2Datatype(std::complex<double> * srcPtr, 
                                                        char * page, 
                                                        DataType outputDataType, 
                                                        size_t pageSize)
{
    switch (outputDataType)
    {
    case ComplexShort:
    {
        std::complex<short> * ptr = (std::complex<short> *) page;
        for(int i=0; i<pageSize;i++) ptr[i] = (std::complex<short>)srcPtr[i];
        break;
    }
    case ComplexInt:
    {
        std::complex<int> * ptr = (std::complex<int> *) page;
        for(int i=0; i<pageSize;i++) ptr[i] = (std::complex<int>)srcPtr[i];
        break;
    }
    case ComplexFloat:
    {
        std::complex<float> * ptr = (std::complex<float> *) page;
        for(int i=0; i<pageSize;i++) ptr[i] = (std::complex<float>)srcPtr[i];
        break;
    }
    case ComplexDouble:
    {
        memcpy(page, srcPtr, pageSize*sizeof(std::complex<double>));
        break;
    }
    default:
    {
        std::cerr<<"Datatype= "<<datatype<<std::endl;
        REPORT_ERROR(16," ERROR: cannot cast datatype to complex<double>");
        break;
    }
    }


}

// Get size of datatype
unsigned long  gettypesize(DataType type)
{
    unsigned long   size;

    switch ( type ) {
        case UChar: case SChar:  size = sizeof(char); break;
        case UShort: case Short: size = sizeof(short); break;
        case Int:                size = sizeof(int); break;
        case Float:              size = sizeof(float); break;
        case Double:             size = sizeof(double); break;
        case ComplexShort:       size = sizeof(std::complex<short>); break;
        case ComplexInt:         size = sizeof(std::complex<int>); break;
        case ComplexFloat:       size = sizeof(std::complex<float>); break;
        case ComplexDouble:      size = sizeof(std::complex<double>); break;
        default: size = 0;
    }

    return(size);
}
