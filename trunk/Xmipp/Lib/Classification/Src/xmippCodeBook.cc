/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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

//-----------------------------------------------------------------------------
// xmippCB.cc
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

#include "../xmippCodeBook.hh"

//-----------------------------------------------------------------------------

/**
 * This class implements a codebook.
 * A codebook is a set of examples (each of them being a vector). These
 * examples are usually labeled, ie, classified (the codebook is calibrated),
 * but it is not necessary (the codebook is NOT calibrated). A codebook is the
 * result obtained after one has trained some kind of algorithms. The way to
 * classify a data once the algorithm has been trained is to look for the
 * example in the code book that best matches the data to classify. Then, the
 * same label of the example in the codebook is associated to the data wanted
 * to be classified (if it is a calibrated codebook), or the example itself is
 * returned (if it is a NO calibrated codebook) indicating with this that the
 * data belongs to the same 'class' that the returned example.
 */


  /**
   * Constructor.
   * Constructs a codebook with initial code vectors at zero.
   * from an unsigned integer to instantiate the template
   * @param _n       Number of vectors
   * @param _size    Size of code vectors
   * @param _lower   Lower value for random elements
   * @param _upper   Upper value for random elements
   * @param _cal     Calibrated or not, that is, a CB with class labels or not
   */
  xmippCB::xmippCB (unsigned _n, unsigned _size, bool _cal = false )
   : xmippCDSet<xmippVector, xmippLabel>(),
     xmippCTSet<xmippVector, xmippLabel>(_cal)
    {
      // Fill vectors with zeros
      theItems.resize(_n);
      //if (calibrated())
	  theTargets.resize(_n, "");
      for ( unsigned i=0 ; i<_n ; i++) {
        xmippVector v;
	v.resize(_size, 0);
	theItems[i] = v;
      }
  };



  /**
   * Constructor.
   * Constructs a codebook with random initial code vectors.
   * from an unsigned integer to instantiate the template
   * @param _n       Number of vectors
   * @param _size    Size of code vectors
   * @param _lower   Lower value for random elements
   * @param _upper   Upper value for random elements
   * @param _cal     Calibrated or not, that is, a CB with class labels or not
   */
  xmippCB::xmippCB (unsigned _n, unsigned _size, xmippFeature _lower = 0, xmippFeature _upper = 1,
         bool _cal = false )
   : xmippCDSet<xmippVector, xmippLabel>(),
     xmippCTSet<xmippVector, xmippLabel>(_cal)
    {
      theItems.resize(_n);
      //if (calibrated())
	  theTargets.resize(_n, "");
      // Assign random vectors
      for ( unsigned i=0 ; i<_n ; i++) {
        vector<xmippFeature> v;
	v = randomVector(_size, _lower, _upper);
	theItems[i] = v;
      }
  };

  /**
   * Constructor.
   * Constructs a codebook with initial code vectors taken randomly from
   * the training file.
   * from an unsigned integer to instantiate the template
   * @param _n       Number of vectors
   * @param _ts	     Training set; will be used to get initial values
   */
  xmippCB::xmippCB (unsigned _n, const xmippCTVectors& _ts )
    : xmippCDSet<xmippVector, xmippLabel>(), xmippCTSet<xmippVector, xmippLabel>(_ts.calibrated())
  {
    // Take random samples
//    xmippUniform<unsigned> chosen( 0, _ts.size() -1 );
    randomize_random_generator();

    theItems.resize(_n);
    //if (_ts.calibrated())
	theTargets.resize(_n, "");
    // Assign random vectors
    for ( unsigned i=0 ; i<_n ; i++) {
        vector<xmippFeature> v;
	int index = (int) rnd_unif(0, _ts.size() -1);
	v = _ts.theItems[index];
	theItems[i] = v;
    }

  };

  /**
   * Constructs a code book given a stream
   * @param _is  The input stream
   * @exception  runtime_error  If there are problems with the stream
   */
  xmippCB::xmippCB(istream& _is) : xmippCDSet<xmippVector, xmippLabel>(), xmippCTSet<xmippVector, xmippLabel>(_is)
  {
     readSelf(_is);
  };

  /**
   * Returns the code vector that represents the input in the codebook
   * @param _in    Sample to classify
   */
  xmippVector& xmippCB::test(const xmippVector& _in) const
  {
    // eval the first one to init best & bestDist
    vector<xmippVector>::const_iterator i = itemsBegin();
    vector<xmippVector>::const_iterator best = i;
    double bestDist = (double) eDist(*i, _in);

    // eval the rest
    for (i++ ; i<itemsEnd() ; i++)
    {
      double dist = (double) eDist(*i, _in);
      if (dist < bestDist)
      {
        bestDist = dist;
        best = i;
      }
    }

    return (xmippVector&)*best;
  };

  /**
   * Returns the index to the code vector that represents the input in the codebook
   * @param _in  Sample to classify
   */
  unsigned xmippCB::testIndex(const xmippVector& _in) const
  {
    // eval the first one to init best & bestDist
    vector<xmippVector>::const_iterator i = itemsBegin();
    vector<xmippVector>::const_iterator best = i;
    double bestDist = (double) eDist(*i, _in);

    // eval the rest
    unsigned bestIndex = 0, index = 0;
    for (i++ ; i<itemsEnd() ; i++)
    {
      double dist = (double) eDist(*i, _in);
      if (dist < bestDist)
      {
        bestDist = dist;
        best = i;
        bestIndex = index;
      }
      index++;
    }

    return bestIndex;
  };

  /**
   * Returns the index of the codevector closest to an input.
   * This is the method used to classify inputs
   * @param _ts  Training set
   * @param _in  Index to the Sample to be classified
   */
  unsigned xmippCB::winner(const xmippCTVectors& _ts, unsigned _in) const
  {
    return testIndex(_ts.theItems[_in]);
  };


  /**
   * Fills the classifVectors with the list of the best input vectors associated to it.
   * @param _ts  Sample list to classify
   */
  void xmippCB::classify(const xmippCTVectors* _ts) 
  {
    	classifVectors.clear(); // clear previous classification.
    	classifVectors.resize(size());
    	aveDistances.clear(); // clear previous classification.
    	aveDistances.resize(size());
	for (unsigned j=0 ; j<_ts->size() ; j++)
	  classifVectors[testIndex(_ts->theItems[j])].push_back(j);

	for (unsigned i=0 ; i< size() ; i++) {
	  double aveDist = 0; 
	  for (unsigned j=0 ; j< classifVectors[i].size() ; j++)
	     aveDist += (double) eDist(theItems[i], _ts->theItems[classifVectors[i][j]]);
	  if (classifVectors[i].size() != 0)   
	     aveDist /= (double) classifVectors[i].size();   
	  aveDistances[i] = (double) aveDist;
	}
	
  };


  /**
   * Prints the histogram values of each Fuzzy codevector. 
   * @param _os  The the output stream
   */
  void xmippCB::printHistogram(ostream& _os) const
  {
     _os << "1 " << size() << endl;
     for (int j = 0; j < size(); j++)
       _os << j << " " << classifSizeAt(j) << endl;
  };
  

  /**
   * Prints the Average Quantization Error of each codevector. 
   * @param _os  The the output stream
   */
  void xmippCB::printQuantError(ostream& _os) const
  {
     _os << "1 " << size() << endl;
     for (int j = 0; j < size(); j++)
       _os << j << " " << aveDistances[j] << endl;
  };

  /**
   * Returns the list of input vectors associated to this code vector.
   * @_index  code vector index
   */
  const vector< unsigned>& xmippCB::classifAt(const unsigned& _index) const
  {
      if (_index < 0 || _index > classifVectors.size())
	  {
		ostrstream msg;
		msg << "index out of range";
		throw runtime_error(msg.str());
	  }
	  return classifVectors[_index];
  };

   /**
   * Returns the number of input vectors associated to this code vector.
   * @_index  code vector index
   */
  unsigned xmippCB::classifSizeAt(const unsigned& _index) const
  {
      if (_index < 0 || _index > classifVectors.size())
	  {
		ostrstream msg;
		msg << "index out of range";
		throw runtime_error(msg.str());
	  }
	  return classifVectors[_index].size();
  };


  /**
   * Returns the label associated to an input
   * @param _in  Sample to classify
   */
  xmippLabel xmippCB::apply(const xmippVector& _in) const
  {
    return theTargets[testIndex(_in)];
  };
  

  /**
   * Calibrates the code book
   * @param _ts   The calibrated training set
   * @param _def  Default target for non-calibrated vectors
   * @exception runtime_error  If the training set is not calibrated
   */
  void xmippCB::calibrate (xmippCTVectors& _ts, 
			  xmippLabel _def = "")
  {
    // set the default label
    for (vector<xmippVector>::const_iterator i = itemsBegin() ;
         i<itemsEnd() ; i++)
      theTargets[i - itemsBegin()] = _def;
    if (_ts.calibrated()) {
       for (unsigned j=0 ; j<_ts.size() ; j++)
          theTargets[testIndex(_ts.theItems[j])] = _ts.theTargets[j];
       calibrated(true );
    } else
       calibrated(false);
   };

   /**
   * Returns the index of the codevector closest to an input.
   * This is the method used to classify inputs
   * @param _in  Sample to classify.
   */
  unsigned xmippCB::output(const xmippVector& _in) const
  {
    return testIndex(_in);
  };


  /**
   * Standard output for a codebook
   * @param _os The output stream
   */
  void xmippCB::printSelf(ostream& _os) const 
  {
	xmippCTSet<xmippVector, xmippLabel>::printSelf(_os);
  };

  /**
   * Standard input for a codebook
   * @param _is The input stream
   */
  void xmippCB::readSelf (istream& _is)
  {

    #ifndef _NO_EXCEPTION
    try
    {
    #endif
      clear();
      string line;

      // Determines the number of rows and columns in the training set
      
      long dim, size;       
      _is >> dim; 
      _is >> line;
      if (!sscanf(line.c_str(),"%d",&size)) {
         int x, y;
	 _is >> x;
	 _is >> y;
	 size = x*y;
      } 
      getline(_is, line);
      theItems.resize(size);
      theTargets.resize(size);

      for (int i = 0; i < size; i++) {
    	 vector<xmippFeature> v;
         v.resize(dim);
         for (int j=0; j < dim; j++){
	    xmippFeature var;  
	    _is >> var;
	    v[j] = var;
	 }
         getline(_is, line);
	 theItems[i] = v;
	 theTargets[i] = line;
      }
    #ifndef _NO_EXCEPTION
    }
    catch (exception& e)
    {
      ostrstream msg;
      msg << e.what() << endl << "Error reading the code book";
      throw runtime_error(msg.str());
    }
    #endif
  };


  /** 
   * Reads the classif vectors from a stream.
   * @param _is  The input stream
   */
    void xmippCB::readClassifVectors(istream& _is) {
    	    int dim;
    	    _is >> dim;
    	    classifVectors.resize(dim);
    	    for (int i = 0; i < classifVectors.size(); i++)
    		    _is >> classifVectors[i];  
    }


  /** 
   * Writes the classif vectors to a stream
   * @param _os  The output stream
   */
    void xmippCB::writeClassifVectors(ostream& _os) const {
    	    _os << classifVectors.size() << endl;
    	    for (int i = 0; i < classifVectors.size(); i++)
    		    _os << classifVectors[i] << endl;  
    }


  /**
   * Saves the xmippCB class into a stream. 
   * this method can be used to save the status of the class.
   * @param _os The output stream
   */
  void xmippCB::saveObject(ostream& _os) const
  {
    writeClassifVectors(_os); 
    xmippCTSet<xmippVector, xmippLabel>::saveObject(_os);
  };


  /**
   * Loads the xmippCB class from a stream. 
   * this method can be used to load the status of the class.
   * @param _is The output stream
   */
  void xmippCB::loadObject(istream& _is)
  {
     clear();
     readClassifVectors(_is); 
     xmippCTSet<xmippVector, xmippLabel>::loadObject(_is);
  };


  /** 
   *	UnNormalize all features in the codebook 
   * 	@param _varStats The normalization information
   */
  
  void xmippCB::unNormalize(const vector<xmippCTVectors::statsStruct>&  _varStats) {
        if (_varStats.size() != theItems[0].size()) {
		ostrstream msg;
		msg << "Normalization information does not coincide with codebook structure";
		throw runtime_error(msg.str());
	}
	for(unsigned it = 0; it < size(); it++) {
		for (unsigned i = 0; i < theItems[0].size(); i++) {
			if (!isnan(theItems[it][i]))
				theItems[it][i] = theItems[it][i]*_varStats[i].sd + _varStats[i].mean;
		}
	}
  }


  /** 
   *	Normalize all features in the codebook 
   * 	@param _varStats The normalization information
   */
  
  void xmippCB::Normalize(const vector<xmippCTVectors::statsStruct>&  _varStats) {
        if (_varStats.size() != theItems[0].size()) {
		ostrstream msg;
		msg << "Normalization information does not coincide with codebook structure";
		throw runtime_error(msg.str());
	}
	for(unsigned it = 0; it < size(); it++) {
		for (unsigned i = 0; i < theItems[0].size(); i++) {
			if (!isnan(theItems[it][i])) {
			    if (_varStats[i].sd != 0)
				theItems[it][i] = (theItems[it][i] - _varStats[i].mean)/_varStats[i].sd;
			}
		}
	}
  }


//-----------------------------------------------------------------------------
