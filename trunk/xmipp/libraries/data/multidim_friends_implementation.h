/* Core array by scalar ---------------------------------------------------- */
template <class T>
void core_array_by_scalar(const maT &op1, const T &op2,
   maT &result, char operation) {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) + op2; break;
            case '-':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) - op2; break;
            case '*':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) * op2; break;
            case '/':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) / op2; break;
            case '^':
               MULTIDIM_ELEM(result,i)=
                  (T) pow((double)MULTIDIM_ELEM(op1,i),(double)op2); break;
         }
   }

/* Scalar by array --------------------------------------------------------- */
template <class T>
void core_scalar_by_array(const T &op1, const maT &op2,
   maT &result, char operation) {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=op1 + MULTIDIM_ELEM(op2,i); break;
            case '-':
               MULTIDIM_ELEM(result,i)=op1 - MULTIDIM_ELEM(op2,i); break;
            case '*':
               MULTIDIM_ELEM(result,i)=op1 * MULTIDIM_ELEM(op2,i); break;
            case '/':
               MULTIDIM_ELEM(result,i)=op1 / MULTIDIM_ELEM(op2,i); break;
            case '^':
               MULTIDIM_ELEM(result,i)=(T)
                  pow((double)op1,(double)MULTIDIM_ELEM(op2,i)); break;
         }
   }

/* Array by array ---------------------------------------------------------- */
template <class T>
void core_array_by_array(const maT &op1, const maT &op2,
   maT &result, char operation) {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(result)
         switch (operation) {
            case '+':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) +
                                       MULTIDIM_ELEM(op2,i); break;
            case '-':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) -
                                       MULTIDIM_ELEM(op2,i); break;
            case '*':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) *
                                       MULTIDIM_ELEM(op2,i); break;
            case '/':
               MULTIDIM_ELEM(result,i)=MULTIDIM_ELEM(op1,i) /
                                       MULTIDIM_ELEM(op2,i); break;
            case '^':
               MULTIDIM_ELEM(result,i)=
                  (T) pow((double)MULTIDIM_ELEM(op1,i),
                          (double)MULTIDIM_ELEM(op2,i)); break;
         }
   }

/* Equality for normal data types ------------------------------------------ */
template <class T>
   bool operator == (const maT &op1, const maT &op2) {return op1.equal(op2);}
