#include "newimage.h"

template<>
void NewImage< std::complex< double > >::castPage2T(char * page, 
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

