%{
#include "../multidimensional_array.h"
%}

%ignore coreArrayByArray<>( const MultidimArray<T>& op1,
                            const MultidimArray<T>& op2,
                            MultidimArray<T>& result,
                            char operation);
%ignore coreArrayByScalar<>(const MultidimArray<T>& op1,
                            const T& op2,
                            MultidimArray<T>& result,
                            char operation);
%ignore coreScalarByArray<>(const T& op1,
                            const MultidimArray<T>& op2,
                            MultidimArray<T>& result,
                            char operation);
%ignore coreArrayByArray;
%ignore coreArrayByScalar;
%ignore coreScalarByArray;
%ignore arrayByArray;
%ignore DBL_EPSILON;

/* rename operator friends*/
%include "../multidimensional_array.h"

PRINT(MultidimArray)

%extend MultidimArray {
    %template(resize) resize<double>; 
    %template(resize) resize<int>; 
    %template(sameShape) sameShape<double>; 
    %template(sameShape) sameShape<int>; 
    %template(initZeros) initZeros<double>; 
    %template(initZeros) initZeros<int>; 
}

%template(MultidimArrayd) MultidimArray<double>;
%template(MultidimArrayi) MultidimArray<int>; 

%pythoncode
%{
    MultidimArray=MultidimArrayd
    MultidimArray=MultidimArrayi
%}
