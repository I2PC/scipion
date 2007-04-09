/** @defgroup MultidimMacros Macros for all multidimensional arrays
 *
 * This macros are defined to allow high speed in critical parts of your
 * program. They shouldn't be used systematically as usually there is no
 * checking on the correctness of the operation you are performing.
 *
 * Speed comes from three facts: first, they are macros and no function call is
 * performed (although most of the critical functions are inline functions),
 * there is no checking on the correctness of the operation (it could be wrong
 * and you are not warned of it), and destination vectors are not returned
 * saving time in the copy constructor and in the creation/destruction of
 * temporary vectors
 */

/** Access to dimension (size).
 * @ingroup MultidimMacros
 */
#ifndef MULTIDIM_SIZE
#define MULTIDIM_SIZE(v) ((v).__dim)
#endif

/** Dimension of the space in which it is defined (1D, 2D or 3D).
 * @ingroup MultidimMacros
 */
#ifndef SPACE_DIM
#define SPACE_DIM(v) ((v).__spcdim)
#endif

/** For all elements in the array.
 * @ingroup MultidimMacros
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * defines an internal index 'i' which ranges the array using its physical
 * definition
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
 *     std::cout << MULTIDIM_ELEM(v, i) << " ";
 * @endcode
 */
#ifndef FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY
#define FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v) \
    for (int i=0; i<(v).__dim; i++)
#endif

/** For all elements in the array (pointer version).
 * @ingroup MultidimMacros
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * uses an external pointer 'ptr' which ranges the array using its physical
 * definition
 *
 * @code
 * T* ptr;
 * FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY_ptr(v)
 *     std::cout << *ptr << " ";
 * @endcode
 */
#ifndef FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY_ptr
#define FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY_ptr(v) \
    for (int i=0, ptr=__m; i<(v).__dim; i++, ptr++)
#endif

/** Vector element: Physical access.
 * @ingroup MultidimMacros
 *
 * Be careful because this is physical access, usually vectors follow the C
 * convention of starting index==0
 *
 * @code
 * MULTIDIM_ELEM(v, 0) = 1;
 * val = MULTIDIM_ELEM(v, 0);
 * @endcode
 */
#ifndef MULTIDIM_ELEM
#define MULTIDIM_ELEM(v, i) ((v).__m[(i)])
#endif

/** Array access.
 * @ingroup MultidimMacros
 *
 * This macro gives you access to the array (T *)
 *
 * @code
 * std::cout << "This is an int *" << MULTIDIM_ARRAY(v) << std::endl;
 * @endcode
 */
#ifndef MULTIDIM_ARRAY
#define MULTIDIM_ARRAY(v) ((v).__m)
#endif

/// @defgroup MultidimFunctions Functions for all multidimensional arrays

/** Conversion from one type to another.
 * @ingroup MultidimFunctions
 *
 * If we have an integer vector and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 */
template<typename T, typename T1>
void type_cast(const maT& v1, maT1& v2)
{
    if (MULTIDIM_SIZE(v1) == 0)
    {
        v2.clear();
        return;
    }

    v2.resize(v1);

    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v2)
        MULTIDIM_ELEM(v2, i) = (T1) MULTIDIM_ELEM(v1, i);
}
