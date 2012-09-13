// Forward declarations
template<typename T>
class Matrix1D;

template<typename T>
class Matrix2D;

template<typename T>
class Matrix3D;

// FIXME Remove this... bad practice
#define vT Matrix1D< T >
#define mT Matrix2D< T >
#define VT Matrix3D< T >

template<typename T>
std::ostream& operator<<(std::ostream& o, const maT& v);

template<typename T>
bool operator==(const maT& op1, const maT& op2);

template<typename T>
void coreArrayByScalar(const maT& op1, const T& op2, maT& result,
                          char operation);

template<typename T>
void coreScalarByArray(const T& op1, const maT& op2, maT& result,
                          char operation);

template<typename T>
void coreArrayByArray(const maT& op1, const maT& op2, maT& result,
                         char operation);

// Specializations
template<>
std::ostream& operator<<(std::ostream& ostrm, const maTC& m);

template<>
void coreArrayByScalar< complex< double > >(const maTC&op1,
                                               const complex< double >& op2,
                                               maTC& result,
                                               char operation);

template<>
void coreScalarByArray< complex< double > >(const complex< double >& op1,
                                               const maTC& op2,
                                               maTC& result,
                                               char operation);

template<>
void coreArrayByArray< complex< double > >(const maTC& op1, const maTC& op2,
                                              maTC& result,
                                              char operation);
