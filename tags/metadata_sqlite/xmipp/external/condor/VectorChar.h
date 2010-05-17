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

#ifndef _INCLUDE_VectorCHAR_H
#define _INCLUDE_VectorCHAR_H

#include <stdlib.h>

class VectorChar
{
  protected:
    int np,extention;
  public: 
    char *p;
    const int &n;

    VectorChar(): np(0), extention(0), p(NULL), n(np) {};
    VectorChar(int _n);
    VectorChar(int _n, int _ext);
    VectorChar(int _n, char *d);
    VectorChar( const VectorChar& P );
    VectorChar(VectorChar *v);
    ~VectorChar();

    void extend();
    void setSize(int _n);
    void exactshape();
    void print();
    
    // only use the following method at your own risks!
    void prepareExtend(int new_extention);
    
//    int &operator [](int i) { return p[i]; };
    inline int sz() {return np;};
    char operator==( const VectorChar& );
	VectorChar& operator=( const VectorChar& P );
    operator char*() const { return p; };
//    operator unsigned*() const { return (unsigned*)p; };
//    int &operator[]( unsigned i) {return p[i];};

    void set(char c);

  private:
    void alloc();

};

#endif

