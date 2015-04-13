/***************************************************************************
 *
* Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
* Updated by:  J.M. de la Rosa Trevin  (jmdelarosa@cnb.csic.es)
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

#include "mean_shift.h"

struct ThreadProcessPlaneArgs
{
    unsigned int myID;
    unsigned int numThreads;
    double sigma_s;
    double sigma_r;
    MultidimArray<double> *input;
    MultidimArray<double> *output;
    bool fast;
};

void * thread_process_plane( void * args )
{
    ThreadProcessPlaneArgs * thrParams = (ThreadProcessPlaneArgs *) args;

    unsigned int myID = thrParams->myID;
    unsigned int numThreads = thrParams->numThreads;
    double sigma_s = thrParams->sigma_s;
    double sigma_r = thrParams->sigma_r;
    MultidimArray<double> &input = *thrParams->input;
    MultidimArray<double> &output = *thrParams->output;
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
    int _3_sigma_s = (int)(3 * sigma_s);
    double _3_sigma_r = 3.0 * sigma_r;
    int y_min = input.startingY();
    int y_max = input.finishingY();
    int x_min = input.startingX();
    int x_max = input.finishingX();
    int z_min = input.startingZ();
    int z_max = input.finishingZ();
    int curr_x, curr_y, curr_z;
    int prev_x, prev_y, prev_z;
    double error;

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
            if( myID == 0 )
            {
                if( z_min < 0 )
                    std::cerr << k - z_min << " ";
                else
                    std::cerr << k + z_min << " ";
            }

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
                    double YcOld=0;
                    double Yc = A3D_ELEM( input, k, i, j);
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

                        for( int z_j = -(int)sigma_s ; z_j <=(int)sigma_s ; z_j++ )
                        {
                            int z2 = zc + z_j;
                            if( z2 >= z_min && z2 <= z_max)
                            {
                                int z_j_2 = z_j * z_j;  // Speed-up

                                for( int y_j = -(int)sigma_s ; y_j <= (int)sigma_s ; y_j++ )
                                {
                                    int y2 = yc + y_j;
                                    if( y2 >= y_min && y2 <= y_max)
                                    {
                                        int y_j_2 = y_j * y_j; // Speed-up

                                        for( int x_j = -(int)sigma_s ; x_j <= (int)sigma_s ; x_j++ )
                                        {
                                            int x2 = xc + x_j;
                                            if( x2 >= x_min && x2 <= x_max)
                                            {

                                                if( x_j*x_j+y_j_2+z_j_2 <= sigma_s*sigma_s )
                                                {
                                                    double Y2 = A3D_ELEM( input, z2, y2, x2);
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
                    }
                    while(shift>3 && iters<100);

                    A3D_ELEM( output, k,i,j ) = Yc;
                }
            }
        }
    }
    else
    {
        for (int k=z_min; k<=z_max; k++)
        {
            if( myID == 0 )
                std::cerr << k << " ";
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
                    curr_I = A3D_ELEM( input, k, i, j);

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
                                                double I_j = A3D_ELEM( input, z_j, y_j, x_j);
                                                double I_dist=fabs(I_j - curr_I);

                                                if( I_dist <= _3_sigma_r )
                                                {
                                                    // Take this point into account
                                                    double eucl_dist = (curr_x - x_j)*(curr_x - x_j)+(curr_y - y_j)*(curr_y - y_j)+(curr_z - z_j)*(curr_z - z_j);
                                                    double aux = exp(-(eucl_dist*inv_2_sigma_s_2) - (I_dist*I_dist*inv_2_sigma_r_2));
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

                        prev_x = (int)round(curr_x);
                        prev_y = (int)round(curr_y);
                        prev_z = (int)round(curr_z);

                        double isum_denom=1.0/sum_denom;
                        curr_x = (int)round(x_sum*isum_denom);
                        curr_y = (int)round(y_sum*isum_denom);
                        curr_z = (int)round(z_sum*isum_denom);
                        curr_I = I_sum*isum_denom;

                        error = fabs(prev_x - curr_x) + fabs(prev_y - curr_y) + fabs(prev_z - curr_z);
                    }
                    while(error>0);

                    A3D_ELEM( output, k,i,j ) = curr_I;
                }
            }
        }
    }
    return NULL;
}

/* Define params ------------------------------------------------------------------- */
void MeanShiftFilter::defineParams(XmippProgram *program)
{
    program->addParamsLine("== Mean shift ==");
    program->addParamsLine("[--mean_shift <hr> <hs> <iter=1>] : Filter based on the mean-shift algorithm");
    program->addParamsLine("                                  :+ hs: Sigma for the range domain");
    program->addParamsLine("                                  :+ hr: Sigma for the spatial domain");
    program->addParamsLine("                                  :+ iter: Number of iterations to be used");
    program->addParamsLine("      alias -t;");
    program->addParamsLine("[--thr <n=1>]                     : Number of processing threads");
    program->addParamsLine("[--fast]                          : Use faster processing (avoid gaussian calculations)");
    program->addParamsLine("[--save_iters]                    : Save result image/volume for each iteration");

}

void MeanShiftFilter::readParams(XmippProgram *program)
{
    sigma_r = program->getDoubleParam("--mean_shift", 0);//hr
    sigma_s = program->getDoubleParam("--mean_shift", 1);//hs
    iters = program->getIntParam("--mean_shift", 2);//iters
    fast = program->checkParam("--fast");
    numThreads = program->getIntParam("--thr");
    save_iters = program->checkParam("--save_iters");
    verbose = program->verbose;
}

void MeanShiftFilter::apply(MultidimArray<double> &input)
{
    MultidimArray<double> output(input);

    pthread_t * th_ids = new pthread_t[numThreads];
    ThreadProcessPlaneArgs * th_args = new ThreadProcessPlaneArgs[numThreads];

    for( int iter = 0 ; iter < iters ; ++iter )
    {
      if (verbose)
        std::cout << formatString("Running iteration %d/%d", iter+1, iters) << std::endl;

        for( int nt = 0; nt < numThreads ; nt++ )
        {
            th_args[nt].myID = nt;
            th_args[nt].numThreads = numThreads;
            th_args[nt].sigma_s = sigma_s;
            th_args[nt].sigma_r = sigma_r;
            th_args[nt].input = &input;
            th_args[nt].output = &output;
            th_args[nt].fast = fast;

            pthread_create( (th_ids+nt), NULL, thread_process_plane, (void *)(th_args+nt));
        }

        for( int nt = 0; nt < numThreads ; nt++ )
        {
            pthread_join( *(th_ids+nt), NULL);
        }

        if( save_iters )
        {
            FileName fn_aux = "xxx"; //fixme
            fn_aux = fn_aux.insertBeforeExtension( std::string("_iter_") + integerToString(iter) );

            if (verbose)
              std::cout << "Saving intermidiate file: " << fn_aux << std::endl;

            output.write( fn_aux );
        }

        //update input for next iteration
        input = output;
    }
}
