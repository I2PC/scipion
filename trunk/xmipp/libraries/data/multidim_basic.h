/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

// This file contains several definitions for statistics and arithmetic
// operations. To use it you have to redirect the internal type maT
// (multidimensional array<T>) to the type you need. For example, for
// a vector library would be
//
//#define maT vector<T>
//
// These definitions are outside because in this way we can reuse the
// module for other libraries (matrices1D, matrices2D, ...). This way
// a further advantage is that all functions are named the same in
// all libraries


public:
    // The array itself
    T* data;

    // Number of elements
    int size;

    // Dimension (1 for vectors, 2 for matrices..)
    int dimension;

    /// @defgroup Statistics Statistics functions.

    /** Print statistics in current line.
     * @ingroup Statistics
     *
     * No end of line character is written after this print out.
     *
     * @code
     * a.compute_stats();
     * std::cout << "Statistics of variable a ";
     * a.print_stats();
     * std::cout << std::endl;
     * @endcode
     */
    void print_stats(std::ostream& out=std::cout) const
    {
        T min, max;
        double avg, dev;

        compute_stats(avg, dev, min, max);

        out.setf(std::ios::showpoint);
        int old_prec = out.precision(7);

        out << " min= "; out.width(9); out << min;
        out << " max= "; out.width(9); out << max;
        out << " avg= "; out.width(9); out << avg;
        out << " dev= "; out.width(9); out << dev;

        out.precision(old_prec);
    }

    /** Maximum of the values in the array.
     * @ingroup Statistics
     *
     * The returned value is of the same type as the type of the array.
     *
     * @code
     * double max = a.compute_max();
     * @endcode
     */
    T compute_max() const
    {
        if (size <= 0)
            return static_cast< T >(0);

        T max = data[0];

        for (int i=0; i<size; i++)
            if (data[i] > max)
                max = data[i];

        return max;
    }

    /** Minimum of the values in the array.
     * @ingroup Statistics
     *
     * The returned value is of the same type as the type of the array.
     *
     * @code
     * double min = a.compute_min();
     * @endcode
     */
    T compute_min() const
    {
        if (size <= 0)
            return static_cast< T >(0);

        T min = data[0];

        for (int i=0; i<size; i++)
            if (data[i] < min)
                min = data[i];

        return min;
    }

    /** Minimum and maximum of the values in the array.
     * @ingroup Statistics
     *
     * As doubles.
     */
    void compute_double_minmax(double& min, double& max) const
    {
        if (size <= 0)
            return;

        min = max = static_cast< double >(data[0]);

        for (int i=0; i<size; i++)
        {
            if (data[i] < min)
                min = static_cast< double >(data[i]);

            if (data[i] > max)
                max = static_cast< double >(data[i]);
        }
    }

    /** Minimum and maximum within region.
     * @ingroup Statistics
     *
     * The region is specified by two corners.
     */
    void compute_double_minmax(double& min,
                               double& max,
                               const matrix1D< double >& corner1,
                               const matrix1D< double >& corner2) const;

    /** Minimum and maximum within region.
     * @ingroup Statistics
     *
     * The region is specified by two corners.
     */
    void compute_double_minmax(double& min,
                               double& max,
                               const matrix1D< int >& corner1,
                               const matrix1D< int >& corner2) const
    {
        matrix1D< double > dcorner1, dcorner2;
        type_cast(corner1, dcorner1);
        type_cast(corner2, dcorner2);

        compute_double_minmax(min, max, dcorner1, dcorner2);
    }

    /** Average of the values in the array.
     * @ingroup Statistics
     *
     * The returned value is always double, independently of the type of the
     * array.
     *
     * @code
     * double avg = a.compute_avg();
     * @endcode
     */
    double compute_avg() const
    {
        if (size <= 0)
            return 0;

        double sum = 0;

        for (int i=0; i<size; i++)
            sum += static_cast< double >(data[i]);

        return sum / size;
    }

    /** Standard deviation of the values in the array.
     * @ingroup Statistics
     *
     * Be careful that the standard deviation and NOT the variance is returned.
     * The returned value is always double, independently of the type of the
     * array.
     *
     * @code
     * ouble dev = a.compute_stddev();
     * @endcode
     */
    double compute_stddev() const
    {
        if (size <= 0)
            return 0;

        double avg = 0, stddev = 0;

        for (int i=0; i<size; i++)
        {
            avg += static_cast< double >(data[i]);
            stddev += static_cast< double >(data[i]) *
                static_cast< double >(data[i]);
        }

        if (size > 1)
        {
            avg /= size;
            stddev = stddev / size - avg * avg;
            stddev *= size / (size - 1);

            // Foreseeing numerical instabilities
            stddev = sqrt(static_cast< double>((ABS(stddev))));
        }
        else
            stddev = 0;

        return stddev;
    }


    /** Compute statistics.
     * @ingroup Statistics
     *
     * The average, standard deviation, minimum and maximum value are
     * returned.
     */
    void compute_stats(double& avg, double& stddev, T& min, T& max) const
    {
        if (size <= 0)
            return;

        avg = 0;
        stddev = 0;

        min = max = data[0];

        for (int i=0; i<size; i++)
        {
            avg += static_cast< double >(data[i]);
            stddev += static_cast< double >(data[i]) *
                static_cast< double >(data[i]);

            if (data[i] > max)
                max = data[i];

            if (data[i] < min)
                min = data[i];
        }

        avg /= size;

        if (size > 1)
        {
            stddev = stddev / size - avg * avg;
            stddev *= size / (size - 1);

            // Foreseeing numerical instabilities
            stddev = sqrt(static_cast< double >(ABS(stddev)));
        }
        else
            stddev = 0;
    }

    /** Compute statistics within region.
     * @ingroup Statistics
     *
     * The region is specified by two corners.
     */
    void compute_stats(double& avg,
                       double& stddev,
                       T& min,
                       T& max,
                       const matrix1D< int >& corner1,
                       const matrix1D< int >& corner2) const
    {
        matrix1D< double > dcorner1, dcorner2;
        type_cast(corner1, dcorner1);
        type_cast(corner2, dcorner2);

        compute_stats(avg, stddev, min, max, dcorner1, dcorner2);
    }

    /** Compute statistics within region.
     * @ingroup Statistics
     *
     * The region is specified by two corners.
     */
    void compute_stats(double& avg,
                       double& stddev,
                       T& min,
                       T& max,
                       const matrix1D< double >& corner1,
                       const matrix1D< double >& corner2) const;

    /** Adjust the range of the array to a given one.
     * @ingroup Statistics
     *
     * A linear operation is performed on the values of the array such that
     * after it, the values of the array are comprissed between the two values
     * set. The actual array is modified itself
     *
     * @code
     * v.range_adjust(0, 1);
     * // The array is now ranging from 0 to 1
     * @endcode
     */
    // This function must be explictly implemented outside
    void range_adjust(T minF, T maxF)
    {
        if (size == 0)
            return;

        T min0 = compute_min();
        T max0 = compute_max();

        // If max0==min0, it means that the vector is a constant one, so the
        // only possible transformation is to a fixed minF
        double slope;
        if (max0 != min0)
            slope = static_cast< double >(maxF - minF) /
                static_cast< double >(max0 - min0);
        else
            slope = 0;

        for (int i=0; i<size; i++)
            data[i] = minF + static_cast< T >(slope *
                static_cast< double >(data[i] - min0));
    }


    /** Adjust the average and stddev of the array to given values.
     * @ingroup Statistics
     *
     * A linear operation is performed on the values of the array such
     * that after it, the average and standard deviation of the array
     * are the two values set. The actual array is modified itself
     *
     * @code
     * v.statistics_adjust(0,1);
     *
     * // The array has got now 0 mean and stddev=1
     * @endcode
     */
    // This function must be explictly implemented outside.
    void statistics_adjust(double avgF, double stddevF)
    {
        double avg0, stddev0;
        double a, b;

        if (size == 0)
            return;

        T min, max;
        compute_stats(avg0, stddev0, min, max);

        if (stddev0 != 0)
            a = static_cast< double >(stddevF) / static_cast< double >(stddev0);
        else
            a = 0;

        b = static_cast< double >(avgF) - a * static_cast< double >(avg0);

        for (int i=0; i<size; i++)
            data[i] = static_cast< T >(a * static_cast< double >(data[i] + b));
    }

    /// @defgroup Arithmethic Arithmethic operations.

    /** @defgroup ArrayByArray Array "by" array operations.
     * @ingroup Arithmethic
     *
     * These are operations that are performed between 2 arrays of the
     * SAME type (two integer vectors, two double matrices, ...). If they
     * are not of the same type you can convert one of the arrays to the
     * desired type using the function type_cast. The result must have been
     * defined to be of the same type as the operands.
     *
     * In this kind of operations each element of array 1 is operated with its
     * homologous in array 2, it is very important that both have got the
     * same size and starting origins. The result has also got the same
     * shape as the two operated arrays and its former content is lost.
     */

    /** Core array by array operation.
     * @ingroup ArrayByArray
     *
     * It assumes that the result is already resized.
     */
    friend void core_array_by_array<>(const maT& op1,
                                      const maT& op2,
                                      maT& result,
                                      char operation);

    /** v3 = v1 + v2.
     * @ingroup ArrayByArray
     */
    maT operator+(const maT& op1) const
    {
        maT tmp;
        array_by_array(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - v2.
     * @ingroup ArrayByArray
     */
    maT operator-(const maT& op1) const
    {
        maT tmp;
        array_by_array(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * v2.
     * @ingroup ArrayByArray
     */
    maT operator*(const maT& op1) const
    {
        maT tmp;
        array_by_array(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / v2.
     * @ingroup ArrayByArray
     */
    maT operator/(const maT& op1) const
    {
        maT tmp;
        array_by_array(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 = v1 ^ v2.
     * @ingroup ArrayByArray
     */
    maT operator^(const maT& op1) const
    {
        maT tmp;
        array_by_array(*this, op1, tmp, '^');
        return tmp;
    }

    /** v3 += v2.
     * @ingroup ArrayByArray
     */
    void operator+=(const maT& op1)
    {
        array_by_array(*this, op1, *this, '+');
    }

    /** v3 -= v2.
     * @ingroup ArrayByArray
     */
    void operator-=(const maT& op1)
    {
        array_by_array(*this, op1, *this, '-');
    }

    /** v3 *= v2.
     * @ingroup ArrayByArray
     */
    void operator*=(const maT& op1)
    {
        array_by_array(*this, op1, *this, 'x');
    }

    /** v3 /= v2.
     * @ingroup ArrayByArray
     */
    void operator/=(const maT& op1)
    {
        array_by_array(*this, op1, *this, '/');
    }

    /** v3 ^= v2.
     * @ingroup ArrayByArray
     */
    void operator^=(const maT& op1)
    {
        array_by_array(*this, op1, *this, '^');
    }

    /** @defgroup ArrayByScalar Array "by" scalar operations.
     * @ingroup Arithmethic
     *
     * These operations are between an array and a scalar (of the same type as
     * the array). The result must have been defined to be of the same type as
     * the operands.
     *
     * In this kind of operations each element of array 1 is operated with the
     * given constant. The result has also got the same shape as the input
     * array and its former content is lost
     */

    /** Array by scalar.
     * @ingroup ArrayByScalar
     *
     * This function must take one vector and a constant, and operate element
     * by element according to the operation required. This is the function
     * which really implements the operations. Simple calls to it perform much
     * faster than calls to the corresponding operators. Although it is
     * supposed to be a hidden function not useable by normal programmers.
     */
    friend void array_by_scalar(const maT& op1,
                                T op2,
                                maT& result,
                                char operation)
    {
        result.resize(op1);
        core_array_by_scalar(op1, op2, result, operation);
    }

    /** Core array by scalar operation.
     * @ingroup ArrayByScalar
     *
     * It assumes that the result is already resized.
     */
    friend void core_array_by_scalar<>(const maT& op1,
                                       const T& op2,
                                       maT& result,
                                       char operation);

    /** v3 = v1 + k.
     * @ingroup ArrayByScalar
     */
    maT operator+(T op1) const
    {
        maT tmp;
        array_by_scalar(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - k.
     * @ingroup ArrayByScalar
     */
    maT operator-(T op1) const
    {
        maT tmp;
        array_by_scalar(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * k.
     * @ingroup ArrayByScalar
     */
    maT operator*(T op1) const
    {
        maT tmp;
        array_by_scalar(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / k.
     * @ingroup ArrayByScalar
     */
    maT operator/(T op1) const
    {
        maT tmp;
        array_by_scalar(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 = v1 ^ k.
     * @ingroup ArrayByScalar
     */
    maT operator^(T op1) const
    {
        maT tmp;
        array_by_scalar(*this, op1, tmp, '^');
        return tmp;
    }

    /** v3 += k.
     * @ingroup ArrayByScalar
     */
    void operator+=(const T& op1)
    {
        array_by_scalar(*this, op1, *this, '+');
    }

    /** v3 -= k.
     * @ingroup ArrayByScalar
     */
    void operator-=(const T& op1)
    {
        array_by_scalar(*this, op1, *this, '-');
    }

    /** v3 *= k.
     * @ingroup ArrayByScalar
     */
    void operator*=(const T& op1)
    {
        array_by_scalar(*this, op1, *this, '*');
    }

    /** v3 /= k.
     * @ingroup ArrayByScalar
     */
    void operator/=(const T& op1)
    {
        array_by_scalar(*this, op1, *this, '/');
    }

    /** v3 ^= k.
     * @ingroup ArrayByScalar
     */
    void operator^=(const T& op1)
    {
        array_by_scalar(*this, op1, *this, '^');
    }

    /** @defgroup ScalarByArray Scalar "by" array operations.
     * @ingroup Arithmethic
     *
     * These operations are between a scalar (of the same type as the array)
     * and an array. The result must have been defined to be of the same type
     * as the operand. The former content of the result array is lost after
     * the operation.
     *
     * In this kind of operations the constant is operated with each element
     * of array 2. The result has also got the same shape as the input array
     * and its former content is lost
     */

    /** Scalar by array.
     * @ingroup ScalarByArray
     *
     * This function must take one scalar and a vector, and operate element by
     * element according to the operation required. This is the function which
     * really implements the operations. Simple calls to it perform much faster
     * than calls to the corresponding operators. Although it is supposed to
     * be a hidden function not useable by normal programmers.
     */
    friend void scalar_by_array(T op1,
                                const maT& op2,
                                maT& result,
                                char operation)
    {
        result.resize(op2);
        core_scalar_by_array(op1, op2, result, operation);
    }

    /** Core array by scalar operation.
     * @ingroup ScalarByArray
     *
     * It assumes that the result is already resized.
     */
    friend void core_scalar_by_array<>(const T& op1,
                                       const maT& op2,
                                       maT& result,
                                       char operation);

    /** v3 = k + v2.
     * @ingroup ScalarByArray
     */
    friend maT operator+(T op1, const maT& op2)
    {
        maT tmp;
        scalar_by_array(op1, op2, tmp, '+');
        return tmp;
    }

    /** v3 = k - v2.
     * @ingroup ScalarByArray
     */
    friend maT operator-(T op1, const maT& op2)
    {
        maT tmp;
        scalar_by_array(op1, op2, tmp, '-');
        return tmp;
    }

    /** v3 = k * v2.
     * @ingroup ScalarByArray
     */
    friend maT operator*(T op1, const maT& op2)
    {
        maT tmp;
        scalar_by_array(op1, op2, tmp, '*');
        return tmp;
    }

    /** v3 = k / v2
     * @ingroup ScalarByArray
     */
    friend maT operator/(T op1, const maT& op2)
    {
        maT tmp;
        scalar_by_array(op1, op2, tmp, '/');
        return tmp;
    }

    /** v3 = k ^ v2
     * @ingroup ScalarByArray
     */
    friend maT operator^(T op1, const maT& op2)
    {
        maT tmp;
        scalar_by_array(op1, op2, tmp, '^');
        return tmp;
    }

    /// @defgroup Size Size management

    /** Core initialize.
     * @ingroup Size
     */
    void core_init()
    {
        data = NULL;
        size = 0;
    }

    /** Allocate array memory.
     * @ingroup Size
     */
    void core_allocate(long int n)
    {
        if (n <= 0)
        {
            clear();
            return;
        }

        data = new T[n];
        if (data != NULL)
            size = n;
        else
            REPORT_ERROR(1001, "Allocate: No space left");
    }

    /** Deallocate.
     * @ingroup Size
     */
    void core_deallocate()
    {
        if (data != NULL)
            delete[] data;

        data = NULL;
        size = 0;
    }

    /** Resize according to a pattern.
     * @ingroup Size
     *
     * This function resize the actual array to the same size and origin
     * as the input pattern. If the actual array is larger than the pattern
     * then the trailing values are lost, if it is smaller then 0's are
     * added at the end
     *
     * @code
     * v2.resize(v1);
     * // v2 has got now the same structure as v1
     * @endcode
     */
    template<typename T1>
    void resize(const maT1 &op1)
    {
        copy_shape(op1);
    }

    /** Get size.
     * @ingroup Size
     *
     * Returns the size of the object in a 3D vector. If the object is a matrix
     * or a vector, then the higher order dimensions will be set to 1, ie,
     * (Xdim, 1, 1) or (Xdim, Ydim, 1).
     */
    void get_size(int* size) const;

    /** Print shape of multidimensional array.
     * @ingroup Size
     *
     * This function shows the size, starting and finishing indexes of the
     * given array. No end of line is printed neither at the beginning nor
     * the end.
     *
     * @code
     * v.print_shape();
     *
     * std::ofstream fh;
     * ...;
     * v.print_shape(fh);
     * @endcode
     */
    void print_shape(std::ostream& out=std::cout) const;

    /** Makes an array empty.
     * @ingroup Size
     *
     * An empty array is the one with size 0, origin at 0, no statistics, ...
     * The memory associated to the previous array is deallocated.
     *
     * @code
     * v.clear();
     * @endcode
     */
    void clear()
    {
        core_deallocate();
        init_shape();
    }

    /** Index outside logical definition region.
     * @ingroup Size
     *
     * This function returns TRUE if the position defined by v is outside the
     * logical definition region of this array. It throws an exception if
     * there are not enough components in v for this array, that is, 1 for
     * vectors, 2 for matrices and 3 for volumes.
     */
    // This function must be explictly implemented outside
    bool outside(const matrix1D< double >& v) const;

    /** True if this object intersects logically the argument array.
     * @ingroup Size
     */
    // This function must be explictly implemented outside
    bool intersects(const maT& m) const;

    /** True if this object intersects logically the argument rectangle.
     * @ingroup Size
     *
     * corner1 is the top-left corner and corner2 the right-bottom corner
     * of the rectangle.
     */
    // This function must be explictly implemented outside
    bool intersects(const matrix1D< double >& corner1,
                    const matrix1D< double >& corner2) const;

    /** True if the given index is in one corner.
     * @ingroup Size
     *
     * It is understood that one corner in logical coordinates. It throws
     * an exception if there are not enough components in v for this array,
     * that is, 1 for vectors, 2 for matrices and 3 for volumes.
     */
    // This function must be explictly implemented outside
    bool isCorner(const matrix1D< double >& v);

    /** True if the given index is in a border.
     * @ingroup Size
     *
     * It is understood that one border in logical coordinates. It throws an
     * exception if there are not enough components in v for this array, that
     * is, 1 for vectors, 2 for matrices and 3 for volumes.
     */
    // This function must be explictly implemented outside
    bool isBorder(const matrix1D< int >& v);

    /** Insert a patch in array.
     * @ingroup Size
     *
     * This function substitutes the array given as argument in this object.
     * It works with logical indexes, only the overlapping area is patched.
     * You can modify the patch behaviour with the following operations '=',
     * '+', '-', '*' or '/'
     */
    // This function must be explictly implemented outside
    void patch(const maT& patch, char operation='=');

    /// @defgroup Initialization Initialization


    /** Same value in all components.
     * @ingroup Initialization
     *
     * The constant must be of a type compatible with the array type, ie,
     * you cannot  assign a double to an integer array without a casting.
     * It is not an error if the array is empty, then nothing is done.
     *
     * @code
     * v.init_constant(3.14);
     * @endcode
     */
    void init_constant(T val)
    {
        for (int i=0; i<size; i++)
            data[i] = val;
    }

    /** Initialize to zeros following a pattern.
     * @ingroup Initialization
     *
     * All values are set to 0, and the origin and size of the pattern are
     * adopted.
     *
     * @code
     * v2.init_zeros(v1);
     * @endcode
     */
    void init_zeros(const maT& op)
    {
        resize(op);
        init_constant(static_cast< T >(0));
    }

    /** Initialize to zeros with current size.
     * @ingroup Initialization
     *
     * All values are set to 0. The current size and origin are kept. It is not
     * an error if the array is empty, then nothing is done.
     *
     * @code
     * v.init_zeros();
     * @endcode
     */
    void init_zeros()
    {
        init_constant(static_cast< T >(0));
    }

    /** Initialize with random values.
     * @ingroup Initialization
     *
     * This function allows you to initialize the array with a set of random
     * values picked from a uniform random distribution or a gaussian one. You
     * must choose two parameters for each, for the uniform distribution they
     * mean the range where to generate the random numbers, while in the
     * gaussian case they are the mean and the standard deviation. By default
     * the uniform distribution is selected. The size and origin of the array
     * are not modified.
     *
     * @code
     * v.init_random(0, 1);
     * // uniform distribution between 0 and 1
     *
     * v.init_random(0, 1, "uniform");
     * // the same
     *
     * v.init_random(0, 1, "gaussian");
     * // gaussian distribution with 0 mean and stddev=1
     * @endcode
     */
    void init_random(double op1, double op2, const std::string& mode="uniform")
    {
        if (mode == "uniform")
            for (int i=0; i<size; i++)
                data[i] = static_cast< T >(rnd_unif(op1, op2));
        else if (mode == "gaussian")
            for (int i=0; i<size; i++)
                data[i] = static_cast< T >(rnd_gaus(op1, op2));
        else
            REPORT_ERROR(1005,
                static_cast< std::string >("Init_random: Mode not supported (" +
                mode + ")"));
    }

    /** Add noise to actual values.
     * @ingroup Initialization
     *
     * This function add some noise to the actual values of the array according
     * to a certain random distribution. You must choose two parameters for
     * each, for the uniform distribution they mean the range where to generate
     * the random numbers, while in the gaussian case they are the mean and the
     * standard deviation. By default the uniform distribution is selected. The
     * size and origin of the array are not modified. The array itself is
     * modified.
     *
     * @code
     * v1.add_noise(0, 1);
     * // uniform distribution between 0 and 1
     *
     * v1.add_noise(0, 1, "uniform");
     * // the same
     *
     * v1.add_noise(0, 1, "gaussian");
     * // gaussian distribution with 0 mean and stddev=1
     * @endcode
     */
    void add_noise(double op1,
                   double op2,
                   const std::string& mode="uniform") const
    {
        if (mode == "uniform")
            for (int i=0; i<size; i++)
                data[i] += static_cast< T >(rnd_unif(op1, op2));
        else if (mode == "gaussian")
            for (int i=0; i<size; i++)
                data[i] += static_cast< T >(rnd_gaus(op1, op2));
        else
            REPORT_ERROR(1005,
                static_cast< std::string >("Add_noise: Mode not supported (" +
                mode + ")"));

    }


    /// @defgroup Operators Operators

    /** Assignment.
     * @ingroup Operators
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    maT& operator=(const maT& op1)
    {
        if (&op1 != this)
        {
            resize(op1);

            for (int i=0; i<size; i++)
                data[i] = op1.data[i];
        }

        return *this;
    }

    /** Unary plus.
     * @ingroup Operators
     *
     * It is used to build arithmetic expressions.
     */
    maT operator+()
    {
        return *this;
    }

    /** Unary minus.
     * @ingroup Operators
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     *
     * @code
     * v1 = -v2;
     * v1 = -v2.transpose();
     * @endcode
     */
    maT operator-() const
    {
        maT tmp(*this);
        tmp.core_unary_minus();
        return tmp;
    }

    /** Core Unary minus.
     * @ingroup Operators
     *
     * This function inverts (-x) each component of the core multidimensional
     * array, the result is stored in this same array.
     */
    void core_unary_minus()
    {
        for (int i=0; i<size; i++)
            data[i] = -data[i];
    }

    /** Equality.
     * @ingroup Operators
     *
     * Used to build logical expressions. It checks exact equality value by
     * value (the difference mustn't be greater than EQUAL_ACCURACY, origins,
     * and size.
     *
     * @code
     * if (v1 == v2)
     *     std::cout << "v1 and v2 are the same";
     * @endcode
     */
    bool equal(const maT& op, double accuracy=XMIPP_EQUAL_ACCURACY) const
    {
        if (!same_shape(op))
            return false;
        else
            return core_equality(op, accuracy);
    }

    /** Operator ==.
     * @ingroup Operators
     *
     * The EQUAL_ACCURACY is used as accuracy measure
     */
    bool friend operator==<>(const maT& op1, const maT& op2);

    /** Operator !=.
     * @ingroup Operators
     *
     * Not ==.
     */
    bool friend operator!=(const maT& op1, const maT& op2)
    {
        return !(op1 == op2);
    }

    /** Core equal.
     * @ingroup Operators
     *
     * This is the equality of the core array.
     */
    bool core_equality(const maT& op2, double accuracy) const
    {
        if (size != op2.size)
            return false;

        for (int i=0; i<size; i++)
            if (ABS(data[i] - op2.data[i]) > accuracy)
                return false;

        return true;
    }

    /** Input from input stream.
     * @ingroup Operators
     *
     * Actual size of the array is used to know how many values must be read.
     *
     * @code
     * v.resize(3);
     * std::cin >> v;
     * @endcode
     */
    // This function must be explictly implemented outside
    friend std::istream& operator>>(std::istream& in, maT& v)
    {
        for (int i=0; i<v.size; i++)
             in >> v.data[i];

        return in;
    }

    /** Read from an ASCII file.
     * @ingroup Operators
     *
     * The array must be previously resized to the correct size.
     */
    void read(const FileName& fn)
    {
        std::ifstream in;
        in.open(fn.c_str(), std::ios::in);

        if (!in)
            REPORT_ERROR(1,
                static_cast< std::string >("MultidimArray::read: File " +
                fn + " not found"));

        in >> *this;
        in.close();
    }

    /** Read from a binary file.
     * @ingroup Operators
     *
     * The array must be previously resized to the correct size.
     */
    void read_binary(const FileName& fn)
    {
        std::ifstream in;
        in.open(fn.c_str(), std::ios::in | std::ios::binary);
        if (!in)
            REPORT_ERROR(1,
            static_cast< std::string>("MultidimArray::read: File " + fn
            + " not found"));

        for (int i=0; i<size; i++)
            in.read(static_cast< char* >(&data[i]), sizeof(T));

        in.close();
    }

    /** Output to an output stream.
     * @ingroup Operators
     *
     * A proper format is used to make all values fit in the same cell width.
     *
     * @code
     * std::cout << v;
     * @endcode
     */
    // This function must be explictly implemented outside
    friend std::ostream& operator<<<>(std::ostream& out, const maT& v);

    /** Write to an ASCII file.
     * @ingroup Operators
     */
    void write(const FileName& fn) const
    {
        std::ofstream out;
        out.open(fn.c_str(), std::ios::out);
        if (!out)
            REPORT_ERROR(1,
                static_cast< std::string >("MultidimArray::write: File " + fn
                + " cannot be opened for output"));

        out << *this;
        out.close();
    }

    /** Write to a binary file.
     * @ingroup Operators
     */
    void write_binary(const FileName& fn) const
    {
        std::ofstream out;

        out.open(fn.c_str(), std::ios::out | std::ios::binary);
        if (!out)
            REPORT_ERROR(1,
                static_cast< std::string >("MultidimArray::write: File " + fn
                + " not found"));

        for (int i=0; i<size; i++)
            out.write(static_cast< char* >(&data[i]), sizeof(T));

        out.close();
    }

    /** Edit with xmipp_editor.
     * @ingroup Operators
     *
     *
     * This function generates a random filename starting with PPP and
     * edits it with xmipp_editor. After closing the editor the file is
     * removed.
     */
    // FIXME This is not good practice.. (system)
    void edit()
    {
        FileName nam;
        nam.init_random(15);

        nam = static_cast< std::string >("PPP" + nam+ ".txt");
        write(nam);

        system((static_cast< std::string >("xmipp_edit -i " + nam +
            " -remove &").c_str()));
    }


    /** @defgroup Iterators Iterators.
     *
     * Iterators are functions which could be applied to all elements, or
     * all columns, or whatever in an array. There is only one iterator that
     * can be applied to any kind of array.
     */

    /** Apply function to each element.
     * @ingroup Iterators
     *
     * You can define a function and apply it to each element of the array.
     * A new array with the same shape as the input one is generated and the
     * old one is not modified. The function applied must take an argument
     * of type T and return a result of type T.
     *
     * @code
     * v2 = v1.for_all(&sin);
     * // If v1 and v2 are double arrays, then v2 is a version of v1 where
     * // all values of v1 have been substituted by their respective sines
     * @endcode
     */
    maT for_all(T (*f)(T)) const
    {
        maT tmp;
        tmp.resize(*this);

        for (int i=0; i<size; i++)
            data[i] = (*f)(data[i]);

        return tmp;
    }

    /** Apply a function to each element and store in this same object.
     * @ingroup Iterators
     */
    void for_all(T (*f)(T))
    {
        for (int i=0; i<size; i++)
            data[i] = (*f)(data[i]);
    }

    /** @defgroup Utilities Utilities.
     *
     * Here you have several easy functions to manage the values of
     * the array.
     */

    /** Dissimilarity.
     * @ingroup Utilities
     *
     * It returns the effective maximum absolute difference between 2 arrays,
     * ie, compute the absolute value of the difference of the arrays, then
     * the histogram of the differences is computed and that value for which
     * the number of higher values is less than the percentil out is returned.
     * By default, only the 0.25% of the highest values are rejected as
     * outliers, although this value could be increased as the number of
     * samples in the array decreases.
     *
     * @code
     * if (dissimilarity(v1, v2) < 1e-6)
     *     std::cout << "v1 and v2 are more or less the same";
     * // the default 0.25% value is used
     *
     * if (dissimilarity(v1, v2, 1) < 1e-6)
     *     ...
     * // Now the 1% of the highest values are rejected to compute the
     * // effective range
     * @endcode
     */
    friend double dissimilarity(maT& op1, maT& op2, double percentil_out)
    {
        maT diff;
        array_by_array(op1, op2, diff, '-');
        diff.ABSnD();

        // FIXME
        //*** histogram1D hist; compute_hist(diff,hist,200);
        // return hist.percentil(100-percentil_out);

        return 1;
    }

    /** Several thresholding.
     * @ingroup Utilities
     *
     * Apply a threshold to the array, the object is modified itself. There
     * are several kinds of thresholding and you must specify it, the values
     * given in the fuction have different meanings according to the threshold
     * applied.
     *
     * abs_above: if |x|>a => x=b
     * abs_below: if |x|<a => x=b
     * above:     if  x >a => x=b
     * below:     if  x <a => x=b
     * range:     if  x <a => x=a   and    if x>b => x=b
     *
     * @code
     * v.threshold("abs_above", 10, 10);
     * // any value whose absolute value is above 10 will be substituted by
     * // -10 (if it is negative) or 10 (if it is positive)
     *
     * v.threshold("abs_below", 0.1, 0);
     * // any value whose absolute value is below 0.1 will be substituted by
     * // -0 (if it is negative) or 0 (if it is positive)
     *
     * v.threshold("above", 10, 10);
     * // any value above 10 will be substituted by 10
     *
     * v.threshold("below", -10, -10);
     * // any value below -10 will be substituted by -10
     *
     * v.threshold("range", 0, 1);
     * // v is "saturated" by values 0 and 1, any value outside this range
     * // will be substituted by its nearest border
     * @endcode
     */
    void threshold(const std::string& type, T a, T b)
    {
        int mode;

        if (type == "abs_above")
            mode = 1;
        else if (type == "abs_below")
            mode = 2;
        else if (type == "above")
            mode = 3;
        else if (type == "below")
            mode = 4;
        else if (type == "range")
            mode = 5;
        else
            REPORT_ERROR(1005,
                static_cast< std::string >("Threshold: mode not supported (" +
                type + ")"));

        for (int i=0; i<size; i++)
            switch (mode)
            {
            case 1:
                if (ABS(data[i]) > a)
                    data[i] = SGN(data[i]) * b;
                break;
            case 2:
                if (ABS(data[i]) < a)
                    data[i] = SGN(data[i]) * b;
                break;
            case 3:
                if (data[i] > a)
                    data[i] = b;
                break;
            case 4:
                if (data[i] < a)
                    data[i] = b;
                break;
            case 5:
                if (data[i] < a)
                    data[i] = a;
                else if (data[i] > b)
                    data[i] = b;
                break;
            }
    }

    /** Count with threshold.
     * @ingroup Utilities
     *
     * This function returns the number of elements meeting the threshold
     * condition.
     */
    long count_threshold(const std::string& type, T a, T b)
    {
        int mode;

        if (type == "abs_above")
            mode = 1;
        else if (type == "abs_below")
            mode = 2;
        else if (type == "above")
            mode = 3;
        else if (type == "below")
            mode = 4;
        else if (type == "range")
            mode = 5;
        else
            REPORT_ERROR(1005,
                static_cast< std::string >("Threshold: mode not supported (" +
                type + ")"));

        long ret = 0;

        for (int i=0; i<size; i++)
            switch (mode)
            {
            case 1:
                if (ABS(data[i]) > a)
                    ret++;
                break;
            case 2:
                if (ABS(data[i]) < a)
                    ret++;
                break;
            case 3:
                if (data[i] > a)
                    ret++;
                break;
            case 4:
                if (data[i] < a)
                    ret++;
                break;
            case 5:
                if (data[i] >= a && data[i] <= b)
                    ret++;
            break;
        }

        return ret;
    }

    /** Substitute a value by another.
     * @ingroup Utilities
     *
     * Substitute an old value by a new one. The accuracy is used to say if
     * the value in the array is equal to the old value. Set it to 0 for
     * perfect accuracy.
     */
    void substitute(T oldv,
                    T newv,
                    double accuracy=XMIPP_EQUAL_ACCURACY)
    {
        for (int i=0; i<size; i++)
            if (ABS(data[i] - oldv) <= accuracy)
                data[i] = newv;
    }

    /** Binarize.
     * @ingroup Utilities
     *
     * This functions substitutes all values in a volume which are greater
     * than val+accuracy by 1 and the rest are set to 0. Use threshold to get a
     * very powerful binarization.
     */
    void binarize(double val=0, double accuracy=XMIPP_EQUAL_ACCURACY)
    {
        for (int i=0; i<size; i++)
            if (data[i] <= val + accuracy)
                data[i] = 0;
            else
                data[i] = 1;
    }

    /** ROUND n-dimensional.
     * @ingroup Utilities
     *
     * Applies a ROUND (look for the nearest integer) to each array element.
     */
    void ROUNDnD()
    {
        for (int i=0; i<size; i++)
            data[i] = ROUND(data[i]);
    }

    /** ROUND n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT ROUNDnD(const maT& a)
    {
        maT tmp(a);
        tmp.ROUNDnD();
        return tmp;
    }

    /** CEILING n-dimensional.
     * @ingroup Utilities
     *
     * Applies a CEILING (look for the nearest larger integer) to each
     * array element.
     */
    void CEILnD()
    {
        for (int i=0; i<size; i++)
            data[i] = CEIL(data[i]);
    }

    /** CEILING n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT CEILnD(const maT& a)
    {
        maT tmp(a);
        tmp.CEILnD();
        return tmp;
    }

    /** FLOOR n-dimensional.
     * @ingroup Utilities
     *
     * Applies a FLOOR (look for the nearest larger integer) to each
     * array element.
     */
    void FLOORnD()
    {
        for (int i=0; i<size; i++)
            data[i] = FLOOR(data[i]);
    }

    /** FLOOR n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT FLOORnD(const maT& a)
    {
        maT tmp(a);
        tmp.FLOORnD();
        return tmp;
    }

    /** ABS n-dimensional.
     * @ingroup Utilities
     *
     * Applies an ABS (absolute value) to each array element.
     */
    void ABSnD()
    {
        for (int i=0; i<size; i++)
            data[i] = ABS(data[i]);
    }

    /** ABS n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT ABSnD(const maT& a)
    {
        maT tmp(a);
        tmp.ABSnD();
        return tmp;
    }

    /** MAX n-dimensional.
     * @ingroup Utilities
     *
     * Each component of the result is the maximum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MAXnD(const maT& v1, const maT& v2, maT& result)
    {
        if (!v1.same_shape(v2))
            REPORT_ERROR(1007, "MAX: arrays of different shape");

        result.resize(v1);
        for (int i=0; i<result.size; i++)
            result.data[i] = MAX(v1.data[i], v2.data[i]);
    }

    /** MAX n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT MAXnD(const maT& v1, const maT& v2)
    {
        maT tmp;
        MAXnD(v1, v2, tmp);
        return tmp;
    }

    /** MIN n-dimensional.
     * @ingroup Utilities
     *
     * Each component of the result is the minimum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MINnD(const maT& v1, const maT& v2, maT& result)
    {
        if (!v1.same_shape(v2))
            REPORT_ERROR(1007, "MIN: arrays of different shape");

        result.resize(v1);
        for (int i=0; i<result.size; i++)
            result.data[i] = MIN(v1.data[i], v2.data[i]);
    }

    /** MIN n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT MINnD(const maT& v1, const maT& v2)
    {
        maT tmp;
        MINnD(v1, v2, tmp);
        return tmp;
    }

    /** Sqrt.
     * @ingroup Utilities
     *
     * Each component of the result is the square root of the original
     * component.
     */
    void SQRTnD()
    {
        for (int i=0; i<size; i++)
            data[i] = static_cast< T >(sqrt(static_cast< double >(data[i])));
    }

    /** Sqrt n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT SQRTnD(const maT& a)
    {
        maT tmp(a);
        tmp.SQRTnD();
        return tmp;
    }

    /** Sum of matrix values.
     * @ingroup Utilities
     *
     * This function returns the sum of all internal values.
     *
     * @code
     * double sum = m.sum();
     * @endcode
     */
    double sum() const
    {
        double sum=0;

        for (int i=0; i<size; i++)
            sum += data[i];

        return sum;
    }

    /** Sum of squared vector values.
     * @ingroup Utilities
     *
     * This function returns the sum of all internal values to the second
     * power.
     *
     * @code
     * double sum2 = m.sum2();
     * @endcode
     */
    double sum2() const
    {
        double sum=0;

        for (int i=0; i<size; i++)
            sum += data[i] * data[i];

        return sum;
    }

    /** Log10.
     * @ingroup Utilities
     *
     * Each component of the result is the log10 of the original components.
     */
    void self_log10()
    {
        for (int i=0; i<size; i++)
            data[i] = static_cast< T >(log10(static_cast< double >(data[i])));
    }

    /** Log10 n-dimensional.
     * @ingroup Utilities
     *
     * The same as before but the result is returned.
     */
    friend maT log10(const maT& a)
    {
        maT tmp(a);
        tmp.self_log10();
        return tmp;
    }

    /** Compute center of mass.
     * @ingroup Utilities
     *
     * If a mask is provided it must be of the same dimension of the object
     * and of type int (i.e., matrix2D<int> *). Only those logical indexes
     * within the object whose mask value is 1 (also in logical indexes)
     * are taken into account.
     */
    void center_of_mass(matrix1D< double >& center, void* mask=NULL);

