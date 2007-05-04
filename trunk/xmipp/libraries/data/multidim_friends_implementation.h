/** Core array by scalar.
 */
template<typename T>
void core_array_by_scalar(const maT& op1,
                          const T& op2,
                          maT& result,
                          char operation)
{
    for (int i=0; i<result.size; i++)
        switch (operation)
        {
        case '+':
            result.data[i] = op1.data[i] + op2;
            break;
        case '-':
            result.data[i] = op1.data[i] - op2;
            break;
        case '*':
            result.data[i] = op1.data[i] * op2;
            break;
        case '/':
            result.data[i] = op1.data[i] / op2;
            break;
        case '^':
            result.data[i] = static_cast< T >(pow(
            static_cast< double >(op1.data[i]), static_cast< double >(op2)));
            break;
        }
}

/** Scalar by array.
 */
template<typename T>
void core_scalar_by_array(const T& op1,
                          const maT& op2,
                          maT& result,
                          char operation)
{
    for (int i=0; i<result.size; i++)
        switch (operation)
        {
        case '+':
            result.data[i] = op1 + op2.data[i];
            break;
        case '-':
            result.data[i] = op1 - op2.data[i];
            break;
        case '*':
            result.data[i] = op1 * op2.data[i];
            break;
        case '/':
            result.data[i] = op1 / op2.data[i];
            break;
        case '^':
            result.data[i] = static_cast< T >(pow(
            static_cast< double >(op1), static_cast< double >(op2.data[i])));
            break;
        }
}

/** Array by array.
 */
template<typename T>
void core_array_by_array(const maT& op1,
                         const maT& op2,
                         maT& result,
                         char operation)
{
    for (int i=0; i<result.size; i++)
        switch (operation)
        {
        case '+':
            result.data[i] = op1.data[i] + op2.data[i];
            break;
        case '-':
            result.data[i] = op1.data[i] - op2.data[i];
            break;
        case '*':
            result.data[i] = op1.data[i] * op2.data[i];
            break;
        case '/':
            result.data[i] = op1.data[i] / op2.data[i];
            break;
        case '^':
            result.data[i] = static_cast< T >(pow(
            static_cast< double >(op1.data[i]),
            static_cast< double >(op2.data[i])));
            break;
        }
}

/** Equality for normal data types.
 */
template<typename T>
bool operator==(const maT& op1, const maT& op2)
{
    return op1.equal(op2);
}
