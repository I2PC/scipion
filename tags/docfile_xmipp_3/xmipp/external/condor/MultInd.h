/*

CONDOR 1.06 - COnstrained, Non-linear, Direct, parallel Optimization 
              using trust Region method for high-computing load, 
              noisy functions
Copyright (C) 2004 Frank Vanden Berghen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation version 2
of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

If you want to include this tools in any commercial product, 
you can contact the author at fvandenb@iridia.ulb.ac.be

*/

//
//	class Multiindex
//
#ifndef _INCLUDE_MULTIND_H_
#define _INCLUDE_MULTIND_H_

#include "VectorInt.h"

class MultInd;

class MultIndCache {
  public:
    MultIndCache();
    ~MultIndCache();
    MultInd *get(unsigned _dim, unsigned _deg);
  private:
    MultInd *head;
};

#ifndef __INSIDE_MULTIND_CPP__
extern MultIndCache cacheMultInd;
#endif

class MultInd {
friend class MultIndCache;
public:
    unsigned dim, deg;

    unsigned *lastChanges();
    unsigned *indexesOfCoefInLexOrder();
    
    MultInd(unsigned d=0);
    ~MultInd();

    void resetCounter();
    MultInd& operator++();       // prefix
    MultInd& operator++( int ) { return this->operator++(); } // postfix
//    unsigned &operator[]( unsigned i) {return coeffDeg[i];};
    inline operator unsigned*() const { return coeffDeg; };
    MultInd& operator=( const MultInd &P );
    bool operator==( const MultInd& m );
    unsigned index() {return indexV;};
    unsigned len();
  
  // Print it
    void print();

private:
    MultInd( unsigned _dim, unsigned _deg );
    void fullInit();
    void standardInit();

    VectorInt lastChangesV, indexesOfCoefInLexOrderV;
    unsigned *coeffDeg, *coeffLex, indexV;
    
    static unsigned *buffer, maxDim;
    // to do the cache:
    MultInd *next;
};


#endif 
