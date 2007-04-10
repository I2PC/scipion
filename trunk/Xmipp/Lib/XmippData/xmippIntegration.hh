/***************************************************************************
 *
 * Authors:     Javier Rodríguez Falces (jrodriguez@cnb.uam.es)
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

#ifndef XMIPPINTEGRATION_H
#define XMIPPINTEGRATION_H
using namespace std;

// TODO remove this, never open namespaces in header files
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <stdio.h>


/** @defgroup Numerical_interation Numerical interation 
   This code performs numeric integrations as described in the 
   Numerical Recipes Book, in particular it  implements the 
   "trapezoidal" (Trapeze) and the Romberg integraation. Both are
   designed for smoothly variant functions. 
 */

/** doubleFunction.
   @ingroup Numerical_interation
   Auxiliary class for numerical integration
  
*/
class doubleFunction{ 						//abstract class
public:
	virtual double operator()()=0; 			//pure virtual function
	virtual ~doubleFunction(){}; 				//virtual destructor
};


/**Fast integration routine. 
  Interpolations are made with lines
   @ingroup Numerical_interation
  
  Example of use.
  
  1) Define function to Numerical_interation as class:
@code
class Func1: public doubleFunction{ 			//actual function to be Numerical_interationd
public: 						//This should be in testinteg
	double x;
	double cte1,cte2;
	virtual double operator()() 			//overloads pure virtual
	{return sqrt(1+cte1*cte1*sin(cte2*x)*sin(cte2*x));}
};
@endcode

 2) In the main code
 
@code
#include <XmippData/xmippIntegration.hh>
Func1 cosine;   //cosine surface
cosine.cte1=fabs(cteA*cteA*cteB*cteB*vx*vx);
cosine.cte2=fabs(cteB*vx);
Trapeze Trap(cosine,cosine.x,inte_low,inte_high);
integralt= Trap();
@endcode

*/
 class Trapeze: public doubleFunction{
	double s;
	doubleFunction &func; 				//the generic function to be Numerical_interationd
	double a,b; 						//integral limits
	double &x; 							//integration variable
	double EPS; 						//desired accuracy
	int JMAX; 							//2**(JMAX)=max number of func. evaluation
public:
	virtual double operator()(); 				//overloads pure virtual function
    /**
	@param min integration lower limit
	@param max integration upper limit
	@param precision maximum error allowed
	@param max_iter maximum number of iterations
    */
	double operator()(double min,double max, 	//redo integral with different
		double precision=1.0e-7,int max_iter=20){ //parameters
		a=min;b=max;EPS=precision;
		JMAX=max_iter; return (*this)();};

    /**
	@param f pointer to function to be integrated
	@param var  integration variable
	@param min integration lower limit
	@param max integration upper limit
	@param precision maximum error allowed
	@param max_iter maximum number of iterations
    */
	Trapeze(doubleFunction &f,double &Var,double min,double max,
			double precision=1.0e-7,int max_iter=20):func(f),x(Var)
		{a=min;b=max;EPS=precision;JMAX=max_iter;};
    /** workhorse that doublely does the integral
	*/
	double Trap(int n); 				
};

/** More accurate integration
   @ingroup Numerical_interation
  More accurate integration than Trapeze with smaller truncation error
  than Trapeze method (interpolation are made with polynomials)

  Example of use.
  
  1) Define function to Numerical_interation as class:
@code
class Func1: public doubleFunction{ 			//actual function to be Numerical_interationd
public: 						//This should be in testinteg
	double x;
	double cte1,cte2;
	virtual double operator()() 			//overloads pure virtual
	{return sqrt(1+cte1*cte1*sin(cte2*x)*sin(cte2*x));}
};
@endcode

 2) In the main code
 
@code
#include <XmippData/xmippIntegration.hh>
Func1 cosine;   //cosine surface
cosine.cte1=fabs(cteA*cteA*cteB*cteB*vx*vx);
cosine.cte2=fabs(cteB*vx);
Romberg Rom(cosine,cosine.x,inte_low,inte_high);
integralt= Rom();
@endcode


 */
 class Romberg : public doubleFunction{
	double s;
	doubleFunction &func; 			//the function to be Numerical_interationd
	double a,b; 					//integral limits
	double &x; 						//integration variable
	double EPS; 					//desired accuracy
public:
	virtual double operator()();
    /**
	@param min integration lower limit
	@param max integration upper limit
	@param precision maximum error allowed
	@param max_iter maximum number of iterations
    */
	double operator()(double min,double max,double precision=1.0e-7)
		{a=min;b=max;EPS=precision;return (*this)();};
    /**
	@param f pointer to function to be integrated
	@param var  integration variable
	@param min integration lower limit
	@param max integration upper limit
	@param precision maximum error allowed
	@param max_iter maximum number of iterations
    */
	Romberg(doubleFunction &f,double &Var,double min,double max,
			double precision=1.0e-7):func(f),x(Var)
	{a=min;b=max;EPS=precision;};
    /** workhorse that doublely does the integral
	*/
	double midpnt(int n); 			
};

void polint(double xa[],double ya[],int n,double x,double &y,double &dy);

double *my_vector(int nl,int nh);

void error(const char *s);
///@}

#endif

 
 
 
 
