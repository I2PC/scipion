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

/* ************************************************************************* */
/* HISTOGRAMS                                                                */
/* ************************************************************************* */
// An histogram is defined as a vector of double to make easier later
// work with them. Although all of its fields are public only those
// in the Library documentation should be used from outside

#ifndef _HISTOGRAM_HH
   #define _HISTOGRAM_HH

#include "xmippMatrices1D.hh"
#include "xmippMatrices2D.hh"
#include "xmippMatrices3D.hh"
// Histograms 1D ===========================================================
/**@name Histograms*/
//@{
/** Histograms with 1 parameter.
   This class of histograms are the usual ones where we want to count the
   number of elements within a certain range of a variable, then we make
   the histogram of that variable. The range is divided into small
   subranges within which the values will be grouped. Any value outside
   the global range will not be counted in the histogram.
   To see exactly which is the division between subranges let's have a
   look on the following example where an histogram between 0 and 2
   is computed with 5 steps.
   \begin{verbatim}

   
   [......)                        [.......]
   [      )                        [       ]
   [      )[......)                [       ]
   [      )[      )        [......)[       ]
   [      )[      )[......)[      )[       ]
   [      )[      )[      )[      )[       ]
   [   0  )[   1  )[   2  )[   3  )[    4  ]
   [      )[      )[      )[      )[       ]
   |---------|---------|---------|---------|
  0.0       0.5       1.0       1.5       2.0
   
   \end{verbatim}
   The border points are 0.0, 0.4, 0.8, 1.2, 1.6 and 2.0. The brackets
   and parenthesis try to represent where the border point belongs to, and
   the numbers within the bars are the index of each bar whithin the
   histogram. The height of each bar is the number of times that a value
   within that subrange has been inserted. Be careful that this is not
   a probability density function (pdf), to be so it should be divided
   by the total number of values inserted.

   The following example shows how to work with the histograms. In it
   we will compute which is the central range within which the 95% of
   the values of a matrix are comprised. This example could be
   even simplified by using the function compute_hist but it has been
   kept like this to show the idea behind the histograms
   
   \begin{verbatim}
   // Variable definition
   histogram1D     hist;
   matrix2D        A(50,50);
   double           min_val, max_val;
   double           eff0, effF;
   
   // Matrix initialisation
   A.init_random(0,100);
   
   // Histogram initialisation
   min_val=A.min();
   max_val=A.max();
   hist.init(min_val,max_val,200);
   
   // Histogram computation
   for (int i=STARTINGY(A); i<=FINISHINGY(A); i++)
      for (int j=STARTINGX(A); j<=FINISHINGX(A); j++)
         hist.insert_value(MAT_ELEM(A,i,j));
   
   // Effective range computation
   eff0=hist.percentil( 2.5);
   effF=hist.percentil(97.5);
   cout << "The effective range goes from " << eff0
        << " to " << effF << endl;
   \end{verbatim}
   
*/
class histogram1D: public matrix1D<double> {
  public:
    // Structure ...........................................................
    double     hmin;        // minimum value of the histogram
    double     hmax;        // maximum value of the histogram 
    double     step_size;   // size of step
    int        no_samples;  // No. of points inside the histogram

    // Procedures ..........................................................
    /** Empty constructor.
        Creates an empty histogram. Before using it you must
        initialise it with init.
        \\ Ex: histogram1D hist;*/
    histogram1D() {clear();}
    
    /** Copy constructor.
        Makes an exact copy of the given histogram into another histogram.
        \\Ex: histogram1D hist2(hist1); */
    histogram1D(const histogram1D &H) {clear(); *this=H;}
    
    /** Empties an histogram.
        Remember to initialise it before using it again.
        \\Ex: hist.clear(); */
    void clear();

   /** Assignment. */
   histogram1D & operator = (const histogram1D &H);

    /** Initialisation of the histogram.
        This is the operation which allows the histogram to be used.
        This should be performed before inserting any value in it. The
        information given to this initialisation is the range within which
        the values will be counted, and the number of steps (discrete bars)
        in this range.If the value is outside this range it will not be
        taken into account although we have asked for its insertion in
        the histogram.
        \\Ex: hist.init(-3,3,100);
        \\---> 100 steps in the range -3...3 */
    void init(double min_val, double max_val, int n_steps);
    
    /** Insert a value within histogram.
        The right interval is chosen according to the initialisation of
        the histogram and the count of elements in that interval is
        incremented by 1. If the value lies outside the global range of
        the histogram nothing is done.
        \\ Ex: hist.insert_value(3); */
    void insert_value(double val);
    
    /** Returns the percentil value.
        This function returns the value within the range for which a
        given percent mass of the total number of elements are below it.
        For instance, if we have 120 values distributed between 0 and 45,
        and we ask for the percentil of 60%, the returned value is that
        within 0 and 45 for which 120*0.6=72 elements are smaller than it.
        \\ Ex: percentil60=hist.percentil(60);*/
    double percentil(double percent_mass);
    
    /** Mass below.
        Returns the number of points which are below a certain value */
    double mass_below(double value);

    /** Mass above.
        Returns the number of points which are above a certain value */
    double mass_above(double value) {return no_samples-mass_below(value);}

    /** Show an histogram.
        The first column is the value associated to each histogram
        measure. The second one is the histogram measure. */
    friend ostream& operator << (ostream &o, const histogram1D &hist);

    /** Write an histogram to disk. */
    void write(const FileName &fn) _THROW;
    
    /**@name Access functions
       This functions are not very common to use, and they allow access
       to the histogram elements or the histogram itself as a vector of
       doubles which we can work with. */
    //@{
    /** Value --> Index.
        Given a value it returns the code of the interval where it should
        be counted. If it is outside the global range of the histogram
        it returns -1
        \\ Ex: hist.val2index(12.3,interval_for_it); */
    void val2index(double v, int &i) const
       {if (v==hmax) i=stepNo()-1;
        else i=(int)FLOOR((v-hmin)/step_size);
        if (i<0 || i>=stepNo()) i=-1;}
   
    /** Index --> Value.
        Given the code of one interval, this function returns the value
        of its starting point (its left border point). If the intervals
        are defined as [a,b), this function returns a.
        \\Ex: hist.index2val(0,beginning_of_interval_0);*/
    void index2val(double i,double &v) const
       {v=hmin+i*step_size;}
    
    /** Minimum value where the histogram is defined.
        \\ Ex: cout << "Minimum value for histogram " << hist.min() << endl;*/
    double hist_min() const {return hmin;}

    /** Maximum value where the histogram is defined.
        \\ Ex: cout << "Maximum value for histogram " << hist.max() << endl;*/
    double hist_max() const {return hmax;}

    /** Step size for the histogram.
        \\ Ex: cout << "Step size of the histogram " << hist.step() << endl;*/
    double step() const {return step_size;}

    /** Number of steps in the histogram.
        \\ Ex: cout << "No. Steps in the histogram " << hist.stepNo() << endl;*/
    int stepNo() const {return XSIZE(*this);}

    /** Number of samples introduced in the histogram.
        \\ Ex: cout << "No. Samples in the histogram " << hist.sampleNo() << endl;*/
    double sampleNo() const {return no_samples;}
    //@}
};

// Functions associated to Histograms 1D ===================================
/**@name Functions related to histograms 1D
*/
//@{
/** Compute histogram of a vector within its minimum and maximum value.
    Given an array as input, this function returns its histogram within
    the minimum and maximum of the array, in this way all the values
    in the array are counted. The array can be of any numerical type
    (short int, int, double, ...) and dimension.
    The number of steps must always be given.
    \\ Ex: histogram1D hist; compute_hist(v,hist,100); */
template <class T>
   void compute_hist(T &array, histogram1D &hist, int no_steps)
   {double min, max; array.compute_double_minmax(min, max);
    compute_hist(array, hist, min, max, no_steps);}

/** Compute histogram of the array within two values.
    Given a array as input, this function returns its histogram within
    two values, the array values outside this range are not counted.
    This can be used to avoid the effect of outliers which causes a
    "compression" in the histogram.
    The array can be of any numerical type
    (short int, int, double, ...). The number of steps must always be
    given.
    \\ Ex: histogram1D hist; compute_hist(v,hist,0,1,100); */
template <class T>
   void compute_hist(const T &v, histogram1D &hist, 
      double min, double max, int no_steps)
   {hist.init(min,max,no_steps);
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
      hist.insert_value(MULTIDIM_ELEM(v,i));}

/** Compute the detectability error between two pdf's.
    The input histograms are expressed as probability density functions
    representing two different objects with the same parameter of
    variation, for instance, the protein and background both defined
    by a grey-level. The first histogram would be the histogram of
    the grey-levels associated to voxels which we know to belong to
    the protein area, while the second histogram is for the grey-level
    of voxels which we know to belong to the background. Then the
    intersection area between the both associated pdf's represents
    the probability of comitting an error when classifying a voxel, ie,
    the probability of saying that a voxel is protein when it really is
    background and viceversa. This function allows you to compute this
    probability error when the two histograms are provided. Be careful
    that this probability usually has got a very low value.
    \\Ex: detect_error=detectability_error(hist1,hist2); */
double detectability_error(const histogram1D &h1, const histogram1D &h2);
//@}

// Histograms 2D ===========================================================
/** Histograms with 2 parameters.
   The histogram with 2 parameters can be regarded as an approximation
   to the joint probability density function of two variables (just
   dividing the histogram by its total mass). Ie, the 2D histogram
   is the count of times that a certain combination of values in two
   variables has ocurred. For example, this 2D histograms can be used
   to plot the projection distribution over the topological sphere,
   in such situation only the first two Euler angles are interesting and
   we could plot how many projections are there with first angle equal
   to 45º and second equal to 0º, and so on covering the whole range for
   each angle.
   
   The 2D histogram is defined, as in the 1D case, by the respective ranges
   for the two variables and the number of intervals on each range.
   The distribution and limits of intervals are just the 2D extension
   of the graph shown in histograms 1D.
   
   @see histogram1D
*/
class histogram2D: public matrix2D<double> {
  public:
    // Structure ...........................................................
    double     imin;         // minimum value of the i axis
    double     imax;         // maximum value of the i axis 
    double     istep_size;   // size of step
    double     jmin;         // minimum value of the j axis
    double     jmax;         // maximum value of the j axis 
    double     jstep_size;   // size of step
    int       no_samples;   // No. of points inside the histogram

    // Procedures ..........................................................
    /** Empty constructor.
        Creates an empty histogram. Before using it you must
        initialise it with init.
        \\ Ex: histogram2D hist; */
    histogram2D() {clear();}
    
    /** Copy constructor.
        Makes an exact copy of the given histogram into another histogram.
        \\Ex: histogram2D hist2(hist1); */
    histogram2D(const histogram2D &H) {*this=H;}
    
    /** Empties an histogram.
        Remember to initialise it before using it again.
        \\Ex: hist.clear(); */
    void clear();

    /** Assignment. */
    histogram2D & operator = (const histogram2D &H);

    /** Initialisation of the histogram.
        This is the operation which allows the histogram to be used.
        This should be performed before inserting any value in it. The
        information given to this initialisation is the range for
        each variable within which
        the values will be counted, and the number of steps (discrete bars)
        in these ranges.If the value is outside the defined ranges it
        will not be
        taken into account although we have asked for its insertion in
        the histogram.
        \\Ex: hist.init(0,90,100,0,360,200);
        \\---> 100 steps in the range V=0...90 and 200 steps for U=0...360 */
    void init(double imin_val, double imax_val, int in_steps,
       double jmin_val, double jmax_val, int jn_steps);

    /** Insert a value within histogram.
        The right interval is chosen according to the initialisation of
        the histogram and the count of elements in that interval is
        incremented by 1. If the value lies outside the global range of
        the histogram nothing is done. Notice that two values are needed,
        these are the two values of the two variables in our example of
        the projection at with first Euler angle=45 and second=0, the
        insertion of this projection in the 2D histogram would be like
        in the following example.
        \\ Ex: hist.insert_value(45,0); */
    void insert_value(double v, double u);

    /** Show an histogram.
        The first column and second column are the (X,Y) coordinates of
        each histogram measure. The third one is the histogram measure. */
    friend ostream& operator << (ostream &o, const histogram2D &hist);

    /** Write an histogram to disk. */
    void write(const FileName &fn) _THROW;
    
    /**@name Access functions
       This functions are not very common to use, and they allow access
       to the histogram elements or the histogram itself as a matrix of
       doubles which we can work with. */
    //@{
    /** Value --> Index.
        Given two values for the two variables it returns the code of
        the interval where it should
        be counted. If it is outside the global range of the histogram
        it returns -1 for that variable (i or j). The following example
        tells you that the interval corresponding to (45,0) is that
        with code (i,j), i.e, you could access to hist()(i,j) to know
        how many projections are there in the same interval as this
        projection.
        \\ Ex: hist.val2index(45,0,i,j); */
    void val2index(double v, double u, int &i, int &j) const
       {if (v==imax) i=IstepNo()-1; else i=(int)FLOOR((v-imin)/istep_size);
        if (u==jmax) j=JstepNo()-1; j=(int)FLOOR((u-jmin)/jstep_size);
        if (i<0 || i>=IstepNo()) i=-1;
        if (j<0 || j>=JstepNo()) j=-1;}

    /** Index --> Value.
        Given the code of one interval, this function returns the value
        of its starting point (its left-top border point, ie, its
        lowest corner). If the intervals
        are defined as [v0,vF) and [u0,uF), this function returns the
        point (v0,u0)
        \\Ex: hist.index2val(5,1,v,u);*/
    void index2val(double i, double j, double &v, double &u) const
       {v=imin+i*istep_size;
        u=jmin+j*jstep_size;}
   
    /** Minimum i value where the histogram is defined.
        \\ Ex: cout << "Minimum value for histogram " << hist.Imin() << endl;*/
    double Ihist_min() const {return imin;}

    /** Maximum i value where the histogram is defined.
        \\ Ex: cout << "Maximum value for histogram " << hist.Imax() << endl;*/
    double Ihist_max() const {return imax;}

    /** Step size in i for the histogram.
        \\ Ex: cout << "Step size of the histogram " << hist.Istep() << endl;*/
    double Istep() const {return istep_size;}

    /** Number of steps in i in the histogram.
        \\ Ex: cout << "No. Steps in the histogram " << hist.IstepNo() << endl;*/
    int IstepNo() const {return YSIZE(*this);}

    /** Minimum j value where the histogram is defined.
        \\ Ex: cout << "Minimum value for histogram " << hist.Jmin() << endl;*/
    double Jhist_min() const {return jmin;}

    /** Maximum j value where the histogram is defined.
        \\ Ex: cout << "Maximum value for histogram " << hist.Jmax() << endl;*/
    double Jhist_max() const {return jmax;}

    /** Step size in j for the histogram.
        \\ Ex: cout << "Step size of the histogram " << hist.Jstep() << endl;*/
    double Jstep() const {return jstep_size;}

    /** Number of steps in j in the histogram.
        \\ Ex: cout << "No. Steps in the histogram " << hist.JstepNo() << endl;*/
    int JstepNo() const {return XSIZE(*this);}

    /** Number of samples introduced in the histogram.
        \\ Ex: cout << "No. Samples in the histogram " << hist.sampleNo() << endl;*/
    int sampleNo() const {return no_samples;}
   //@}
};

// Functions associated to Histograms 2D ===================================
/**@name Functions related to histograms 2D
    \\The vectors can be of any numerical type
    (short int, int, double, ...), but both of the same type.
    \\Vectors must be of the same shape, the first element of v1 and
    the first of v2 define the position were the first point will be
    inserted in the histogram, then the second of v1 and of v2, ...
    That is, the elements of v1 and v2 serve as coordinates within the
    histogram.
    \\The number of steps must always be given.
*/
//@{
/** Compute histogram of two arrays within their minimum and maximum values.
    Given two arrays as input, this function returns their joint histogram
    within
    the minimum and maximum of the arrays, in this way all the values
    in the arrays are counted. Both arrays must have the same shape */
template <class T>
   void compute_hist(const T &v1, const T &v2,
      histogram2D &hist, int no_steps1, int no_steps2)
      {double min1, max1; v1.compute_double_minmax(min1, max1);
       double min2, max2; v2.compute_double_minmax(min2, max2);
       compute_hist(v1,v2,hist,min1,max1,min2,max2,no_steps1,no_steps2);}

/** Compute histogram of two arrays within given values.
    Given two arrays as input, this function returns their joint histogram
    within the specified values, all the values lying outside are not counted */
template <class T>
   void compute_hist(const T &v1, const T &v2, histogram2D &hist, 
      double m1, double M1, double m2, double M2, int no_steps1, int no_steps2)
      _THROW
   {if (!v1.same_shape(v2))
      REPORT_ERROR(1,"compute_hist: v1 and v2 are of different shape");
    hist.init(m1,M1,no_steps1,m2,M2,no_steps2);
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v1)
       hist.insert_value(MULTIDIM_ELEM(v1,i),MULTIDIM_ELEM(v2,i));}
//@}
//@}
#endif
