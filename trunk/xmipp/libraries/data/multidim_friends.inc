// Forward declarations
template<typename T>
class matrix1D;

template<typename T>
class matrix2D;

template<typename T>
class matrix3D;

// FIXME Remove this... bad practice
#define vT matrix1D< T >
#define mT matrix2D< T >
#define VT matrix3D< T >

template<typename T>
std::ostream& operator<<(std::ostream& o, const maT& v);

template<typename T>
bool operator==(const maT& op1, const maT& op2);

template<typename T>
void core_array_by_scalar(const maT& op1, const T& op2, maT& result,
                          char operation);

template<typename T>
void core_scalar_by_array(const T& op1, const maT& op2, maT& result,
                          char operation);

template<typename T>
void core_array_by_array(const maT& op1, const maT& op2, maT& result,
                         char operation);

// Specializations
template<>
std::ostream& operator<<(std::ostream& ostrm, const maTC& m);

template<>
void core_array_by_scalar< complex< double > >(const maTC&op1,
                                               const complex< double >& op2,
                                               maTC& result,
                                               char operation);

template<>
void core_scalar_by_array< complex< double > >(const complex< double >& op1,
                                               const maTC& op2,
                                               maTC& result,
                                               char operation);

template<>
void core_array_by_array< complex< double > >(const maTC& op1, const maTC& op2,
                                              maTC& result,
                                              char operation);
