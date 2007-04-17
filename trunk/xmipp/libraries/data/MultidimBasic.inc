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
 * MERCHANTABILITY or FITNESS FOR Aouble_ PARTICULAR PURPOSE.  See the       
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

/* ------------------------------------------------------------------------- */
/* MULTIDIM BASIC                                                            */
/* ------------------------------------------------------------------------- */

/* Core array by scalar ---------------------------------------------------- */
template <>
void core_array_by_scalar< complex<double> >(const maTC &op1,
   const complex<double> &op2, maTC &result, char operation) {
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
                  pow(MULTIDIM_ELEM(op1,i),op2); break;
         }
   }

/* Scalar by array --------------------------------------------------------- */
template <>
void core_scalar_by_array< complex<double> >(const complex<double> &op1,
   const maTC &op2, maTC &result, char operation) {
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
               MULTIDIM_ELEM(result,i)=pow(op1,MULTIDIM_ELEM(op2,i)); break;
         }
   }

/* Array by array ---------------------------------------------------------- */
template <>
void core_array_by_array< complex<double> >(const maTC &op1, const maTC &op2,
   maTC &result, char operation) {
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
               MULTIDIM_ELEM(result,i)=pow(MULTIDIM_ELEM(op1,i),
                                       MULTIDIM_ELEM(op2,i)); break;
         }
   }
