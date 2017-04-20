/***************************************************************************
 *
 * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
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

#include "../../../libraries/reconstruction_adapt_cuda/gpu_rotate_image.h"


//RUN_XMIPP_PROGRAM(ProgGpuRotateImage)

/*double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}*/



int main(int argc, char** argv) {

	struct timeval t_inigpu, t_fingpu;
	double secsgpu;
	gettimeofday(&t_inigpu, NULL);

	ProgGpuRotateImage program;
	program.read(argc, argv);
    int errorCode = program.tryRun();

    gettimeofday(&t_fingpu, NULL);
    secsgpu = timeval_diff(&t_fingpu, &t_inigpu);
    printf("Total Rotate Image GPU: %.8g milisegundos\n", secsgpu * 1000.0);

    //AJ to analyze the obtained computing times
    double timeCpuMemIn=0;
    double timeCpuMemOut=0;
    double timeGpuMemInSpline=0;
    double timeGpuFiltSpline=0;
    double timeGpuKernelSpline=0;
    double timeGpuMemOutSpline=0;
    double timeGpuMemInLinear=0;
    double timeGpuKernelLinear=0;
    double timeGpuMemOutLinear=0;
    if(mytimes[0].calcTime==true){
    	for(int i=0; i<mytimes[0].size; i++){
    		timeCpuMemIn+=mytimes[i].secs_cpu_mem_in;
    		timeCpuMemOut+=mytimes[i].secs_cpu_mem_out;
    		if(mytimes[0].spline==true){
    			timeGpuMemInSpline+=mytimes[i].secs_gpu_mem_inS;
    			timeGpuFiltSpline+=mytimes[i].secs_gpu_filtS;
    			timeGpuKernelSpline+=mytimes[i].secs_gpu_kernelS;
    			timeGpuMemOutSpline+=mytimes[i].secs_gpu_mem_outS;
    		}else{
    			timeGpuMemInLinear+=mytimes[i].secs_gpu_mem_in;
    			timeGpuKernelLinear+=mytimes[i].secs_gpu_kernel;
    			timeGpuMemOutLinear+=mytimes[i].secs_gpu_mem_out;
    		}
    	}

        printf("Total CPU time for memory in-transactions: %.8g miliseconds\n", timeCpuMemIn * 1000.0);
        printf("Total CPU time for memory out-transactions: %.8g miliseconds\n", timeCpuMemOut * 1000.0);
        if(mytimes[0].spline==true){
        	printf("Total GPU time for memory in-transactions (Spline): %.8g miliseconds\n", timeGpuMemInSpline * 1000.0);
        	printf("Total GPU time for filtering kernel (Spline): %.8g miliseconds\n", timeGpuFiltSpline * 1000.0);
        	printf("Total GPU time for interpolation kernel (Spline): %.8g miliseconds\n", timeGpuKernelSpline * 1000.0);
        	printf("Total GPU time for memory out-transactions (Spline): %.8g miliseconds\n", timeGpuMemOutSpline * 1000.0);
        }else{
        	printf("Total GPU time for memory in-transactions (Linear): %.8g miliseconds\n", timeGpuMemInLinear * 1000.0);
        	printf("Total GPU time for interpolation kernel (Linear): %.8g miliseconds\n", timeGpuKernelLinear * 1000.0);
        	printf("Total GPU time for memory out-transactions (Linear): %.8g miliseconds\n", timeGpuMemOutLinear * 1000.0);
        }
        printf("  \n");

    	timeCpuMemIn/=mytimes[0].size;
    	timeCpuMemOut/=mytimes[0].size;
    	if(mytimes[0].spline==true){
    		timeGpuMemInSpline/=mytimes[0].size;
    		timeGpuFiltSpline/=mytimes[0].size;
			timeGpuKernelSpline/=mytimes[0].size;
			timeGpuMemOutSpline/=mytimes[0].size;
    	}else{
    		timeGpuMemInLinear/=mytimes[0].size;
			timeGpuKernelLinear/=mytimes[0].size;
			timeGpuMemOutLinear/=mytimes[0].size;
    	}

    	printf("Average CPU time for memory in-transactions: %.8g miliseconds\n", timeCpuMemIn * 1000.0);
    	printf("Average CPU time for memory out-transactions: %.8g miliseconds\n", timeCpuMemOut * 1000.0);
    	if(mytimes[0].spline==true){
    		printf("Average GPU time for memory in-transactions (Spline): %.8g miliseconds\n", timeGpuMemInSpline * 1000.0);
    		printf("Average GPU time for filtering kernel (Spline): %.8g miliseconds\n", timeGpuFiltSpline * 1000.0);
    		printf("Average GPU time for interpolation kernel (Spline): %.8g miliseconds\n", timeGpuKernelSpline * 1000.0);
    		printf("Average GPU time for memory out-transactions (Spline): %.8g miliseconds\n", timeGpuMemOutSpline * 1000.0);
    	}else{
    		printf("Average GPU time for memory in-transactions (Linear): %.8g miliseconds\n", timeGpuMemInLinear * 1000.0);
    		printf("Average GPU time for interpolation kernel (Linear): %.8g miliseconds\n", timeGpuKernelLinear * 1000.0);
    		printf("Average GPU time for memory out-transactions (Linear): %.8g miliseconds\n", timeGpuMemOutLinear * 1000.0);
    	}

    }


    return errorCode;

}


