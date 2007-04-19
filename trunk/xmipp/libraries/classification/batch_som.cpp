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

/*-----------------------------------------------------------------------------
 xmippBatchSOM.cc
 Implements Kohonen Self-Organizing Feature Maps by using Batch SOM
-----------------------------------------------------------------------------*/

#include "batch_som.h"

 /**
  * Construct a BatchSOM from the code vectors in a stream
  * @param _is  The stream
  */
  xmippBatchSOM::xmippBatchSOM(istream& _is): xmippSOM(_is)
  {
    readSelf(_is);
  };



  /**
   * Trains the SOM
   * @param _som  The som to train
   * @param _ts   The training set
   */
  void xmippBatchSOM::train (xmippMap& _som, const xmippCTVectors& _ts) const
  {


    unsigned long t=0;

    int verbosity = listener->getVerbosity();
    if (verbosity)
  	listener->OnReportOperation((string) "Batch Training Kohonen SOM....\n");
    if (verbosity == 1 || verbosity == 3)
      	listener->OnInitOperation(somNSteps);

    SomIn aveVector(_som.theItems[0].size());	
    vector<unsigned> tmpVector;

    while (t < somNSteps)
    {
	    _som.classify(&_ts);
	    // Check for each SOM unit
	    for (unsigned it = 0; it < _som.size(); it++) {
	  	for (unsigned a = 0; a < aveVector.size(); a++)
	  	   aveVector[a] = 0.;
	  	long total = 0;
	  	   // Collects a list of the input vectors assigned to the neighborhood
	  	vector<unsigned> neig=_som.neighborhood(_som.indexToPos(it), ceil(somRadius(t, somNSteps)));
	  	for (vector<unsigned>::iterator itt = neig.begin();itt<neig.end();itt++)
	  	{
	  		tmpVector =  _som.classifAt(*itt);	
	  		for (unsigned j=0 ; j < tmpVector.size() ; j++) {
	  			SomIn v = _ts.theItems[tmpVector[j]];
	  			for (unsigned a = 0; a < v.size(); a++)
	  				aveVector[a] += v[a];
	  			total++;
	  		}

	  	 }
		 if (total != 0) {	
	  	 	for (unsigned a = 0; a < aveVector.size(); a++)
	  	     		aveVector[a] /= (xmippFeature) total;
	  	 	_som.theItems[it] = aveVector;
		 }
	    }

      	    if (verbosity == 1 || verbosity == 3)
      		listener->OnProgress(t);
      	    if (verbosity >= 2) {
	   	char s[100];
	   	sprintf(s, "Iteration %d of %d.\n", t+1, somNSteps);
  	   	listener->OnReportOperation((string) s);
      	    }	
	    t++;
	
    } // while t < somSteps

    if (verbosity == 1 || verbosity == 3)
        listener->OnProgress(somNSteps);

  };


