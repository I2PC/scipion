/***************************************************************************
 *
 * Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/funcs.h>
#include <data/args.h>
#include <data/image.h>
#include <data/threads.h>
#include <data/histogram.h>

typedef struct threadProcessPlaneArgs
{        
    unsigned int myID;
    unsigned int numThreads;
    double sigma_s;
    double sigma_r;
    Image<double> *inputVolume;
    Image<double> *resultVolume;
    bool fast;
};        

void * thread_process_plane( void * args )
{
    threadProcessPlaneArgs * thrParams = (threadProcessPlaneArgs *) args;
    
    unsigned int myID = thrParams->myID;
    unsigned int numThreads = thrParams->numThreads;
    double sigma_s = thrParams->sigma_s;
    double sigma_r = thrParams->sigma_r;
    Image<double> * inputVolume = thrParams->inputVolume;
    Image<double> * resultVolume = thrParams->resultVolume;
    bool fast = thrParams->fast;
    
    if( !fast )
    {
        sigma_s /= 3.0;
        sigma_s = CEIL( sigma_s );
        sigma_r /= 3.0;
    }
    else
    {
        sigma_s = CEIL( sigma_s );
    }
    
    double curr_I;
    double inv_2_sigma_s_2 = 1.0/( 2.0* sigma_s * sigma_s);
    double inv_2_sigma_r_2 = 1.0/( 2.0* sigma_r * sigma_r);
    double sigma_r_2 = sigma_r * sigma_r;
    int _3_sigma_s = 3 * sigma_s;
    double _3_sigma_r = 3.0 * sigma_r;
	int y_min = inputVolume->data.startingY();
    int y_max = inputVolume->data.finishingY();
    int x_min = inputVolume->data.startingX();
    int x_max = inputVolume->data.finishingX();
    int z_min = inputVolume->data.startingZ();
    int z_max = inputVolume->data.finishingZ();
	int curr_x, curr_y, curr_z;
	int prev_x, prev_y, prev_z, prev_I;
    double error;
    int plane=0;

    int myFirstY, myLastY;
    
    int numElems = y_max-y_min+1;
    numElems /= numThreads;
    
    myFirstY = myID * numElems + y_min;
    
    if( myID == numThreads -1 )
        myLastY = y_max;
    else
        myLastY = myFirstY + numElems - 1;
    
    if( myID == 0 )
    {
        std::cerr << "Progress (over a total of " << z_max - z_min + 1 << " slices): ";
    }
    
    if( fast )
    {
        for (int k=z_min; k<=z_max; k++) 
        {
            if( myID == 0 ) std::cerr << k + z_min << " ";
            for (int i=myFirstY; i<=myLastY; i++) 
            {
                for (int j=x_min; j<=x_max; j++)
	            {
                    // x == j
		            // y == i
		            // z == k
                    
                    int xc = j;
                    int yc = i;
                    int zc = k;
                    
                    int xcOld, ycOld, zcOld;
                    double YcOld;
                    double Yc = A3D_ELEM( inputVolume->data, k, i, j);
                    int iters =0;
                    double shift;
                    int num;
                    double mY;
                    
                    do
		            {
                        xcOld = xc;
                        ycOld = yc;
                        zcOld = zc;
                        
                        double mx=0;
                        double my=0;
                        double mz=0;
                        mY=0;
                        num=0;
                       
			            for( int z_j = -sigma_s ; z_j <=sigma_s ; z_j++ )
                        {
                            int z2 = zc + z_j;
                            if( z2 >= z_min && z2 <= z_max)
                            {
                                int z_j_2 = z_j * z_j;  // Speed-up
                            
        		                for( int y_j = -sigma_s ; y_j <= sigma_s ; y_j++ )
        		                {
                                    int y2 = yc + y_j;
		    		                if( y2 >= y_min && y2 <= y_max)
			                        {
                                        int y_j_2 = y_j * y_j; // Speed-up
                                        
        				                for( int x_j = -sigma_s ; x_j <= sigma_s ; x_j++ )
		    			                {
                                            int x2 = xc + x_j;
			        		                if( x2 >= x_min && x2 <= x_max)
					                        {
                                            
                                                if( x_j*x_j+y_j_2+z_j_2 <= sigma_s*sigma_s )
                                                {
                                                    double Y2 = A3D_ELEM( inputVolume->data, z2, y2, x2);
                                                    double dY = Yc - Y2;
                                                                                                        
                                                    if( dY*dY <= sigma_r_2 )
                                                    {
                                                        mx += x2;
                                                        my += y2;
                                                        mz += z2;
                                                        mY += Y2;
                                                        num++;
                                                    }
                                                }
	    					                }
    	    			                }
    	    		                }
			                    }
                            }
                        }

                        double num_ = 1.0/num;
                        Yc = mY*num_;
                        xc = (int) (mx*num_+0.5);
                        yc = (int) (my*num_+0.5);
                        zc = (int) (mz*num_+0.5);
                        int dx = xc-xcOld;
                        int dy = yc-ycOld;
                        int dz = zc-zcOld;
                        
                        double dY = Yc-YcOld;
                        
                        shift = dx*dx+dy*dy+dz*dz*dY*dY;
                        iters++;
		            }while(shift>3 && iters<100);

                    A3D_ELEM( resultVolume->data, k,i,j ) = Yc;
                }
	        }
        }    
    }
    else
    {    
        for (int k=z_min; k<=z_max; k++) 
        {
            if( myID == 0 ) std::cerr << k << " ";
            for (int i=myFirstY; i<=myLastY; i++) 
            {
                for (int j=x_min; j<=x_max; j++)
	            {
                    // x == j
		            // y == i
		            // z == k

		            curr_x = j;
		            curr_y = i;
                    curr_z = k;
		            curr_I = A3D_ELEM( inputVolume->data, k, i, j);

                    double sum_denom;
                    double I_sum;
		            do
		            {
			            // Check neighbourhood
			            double x_sum = 0, y_sum = 0, z_sum = 0;
		                sum_denom = 0;
                        I_sum = 0;

                        for( int z_j = curr_z - _3_sigma_s ; z_j <= curr_z + _3_sigma_s ; z_j++ )
                        {
                            if( z_j >= z_min && z_j <= z_max)
                            {
        		                for( int y_j = curr_y - _3_sigma_s ; y_j <= curr_y + _3_sigma_s ; y_j++ )
        		                {
		    		                if( y_j >= y_min && y_j <= y_max)
			                        {
        				                for( int x_j = curr_x - _3_sigma_s ; x_j <= curr_x + _3_sigma_s ; x_j++ )
		    			                {
			        		                if( x_j >= x_min && x_j <= x_max)
					                        {
                                                double I_j = A3D_ELEM( inputVolume->data, z_j, y_j, x_j);

                                                if( ABS( I_j - curr_I ) <= _3_sigma_r )
			    				                {							
                                 	                // Take this point into account
					        		                double eucl_dist = (curr_x - x_j)*(curr_x - x_j)+(curr_y - y_j)*(curr_y - y_j)+(curr_z - z_j)*(curr_z - z_j);
        							                double aux_0 = exp(-(eucl_dist) * inv_2_sigma_s_2);
		    						                double aux_1 = exp(-((curr_I-I_j)*(curr_I-I_j)) * inv_2_sigma_r_2 );

                                                    double aux = aux_0*aux_1;

			        				                x_sum += x_j*aux;
					        		                y_sum += y_j*aux;
                                                    z_sum += z_j*aux;
							                        I_sum += I_j*aux;

                                  	                sum_denom += aux;
		    	    			                }
	    					                }
    	    			                }
    	    		                }
			                    }
                            }
                        }

			            prev_x = ROUND(curr_x);
			            prev_y = ROUND(curr_y);
                        prev_z = ROUND(curr_z);
			            prev_I = curr_I;

			            curr_x = ROUND(x_sum / sum_denom);
			            curr_y = ROUND(y_sum / sum_denom);
			            curr_z = ROUND(z_sum / sum_denom);
			            curr_I = I_sum / sum_denom;

                        error = ABS(prev_x - curr_x) + ABS(prev_y - curr_y) + ABS(prev_z - curr_z);

		            }while(error>0);

                    A3D_ELEM( resultVolume->data, k,i,j ) = curr_I;
                }
	        }
        }
    }
    return NULL;
}

int main(int argc, char **argv)
{
	FileName inputFile, outputFile;
	double sigma_r,sigma_s;
    unsigned int numThreads; 
    bool fast;
    
	try
    {
        if ( checkParameter( argc, argv, "-i") )
            inputFile = getParameter( argc, argv, "-i");
        if ( checkParameter( argc, argv, "-o") )
            outputFile = getParameter( argc, argv, "-o");
        if  ( checkParameter( argc, argv, "-hr") )
            sigma_r = textToDouble( getParameter( argc, argv, "-hr") );
        if  ( checkParameter( argc, argv, "-hs") )
            sigma_s = textToDouble( getParameter( argc, argv, "-hs") );
        if  ( checkParameter( argc, argv, "-fast") )
            fast = true;
        else
            fast = false;
        if  ( checkParameter( argc, argv, "-thr") )
            numThreads = textToInteger( getParameter( argc, argv, "-thr") );
        else
            numThreads = 1;
    }
    catch (Xmipp_error)
    {
        std::cout << "Mean_shift: Filters an image based on the mean-shift algorithm" << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "-i                : Input file" << std::endl;
        std::cout << "-o                : Output file" << std::endl;
        std::cout << "-hr               : Sigma for the range domain" << std::endl;
        std::cout << "-hs               : Sigma for the spatial domain" << std::endl;

        exit(1);
    }

	Image<double> inputVolume;
    Image<double> resultVolume;
    
	inputVolume.read(inputFile);
	inputVolume.data.setXmippOrigin();

    resultVolume = inputVolume;
    resultVolume.data.setXmippOrigin();

    pthread_t * th_ids = new pthread_t[numThreads];
    threadProcessPlaneArgs * th_args = new threadProcessPlaneArgs[numThreads];
        
    for( int nt = 0; nt < numThreads ; nt++ )
    {
        th_args[nt].myID = nt;
        th_args[nt].numThreads = numThreads;
        th_args[nt].sigma_s = sigma_s;
        th_args[nt].sigma_r = sigma_r;
        th_args[nt].inputVolume = &inputVolume;
        th_args[nt].resultVolume = &resultVolume;
        th_args[nt].fast = fast;
    
        pthread_create( (th_ids+nt), NULL, thread_process_plane, (void *)(th_args+nt));
    }
    
    for( int nt = 0; nt < numThreads ; nt++ )
    {
        pthread_join( *(th_ids+nt), NULL);
    }
    
    resultVolume.write( outputFile );
    
    std::cout << "DONE!" << std::endl;
}
