#include <vector>
#include <sstream>
#include <fstream>
#include <time.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/video/video.hpp"
#include "opencv2/gpu/gpu.hpp"

#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <data/normalize.h>


using namespace std;
using namespace cv;
#ifdef GPU
using namespace cv::gpu;
#endif
class ProgOpticalAligment: public XmippProgram
{

	public:	    
	    String fname;
	    String foname;
	    int winSize;
	    int gpuDevice; 


	void defineParams()
	{
	        addUsageLine ("Align moviews using optical flow");
	        addParamsLine("     -i <inMoviewFnName>          : input movie File Name");
	        addParamsLine("     -o <outAverageMoviewFnName>  : output aligned micrograhp File Name");
	        addParamsLine("     [--winSize <int=150>]     : window size for optical flow algorithm");
	#ifdef GPU
	        addParamsLine("     [--gpu <int=0>]         : GPU device to be used");
	#endif
	} 
        void readParams()	{
	    fname     = getParam("-i");
	    foname    = getParam("-o");
	    winSize   = getIntParam("--winSize");
	#ifdef GPU
	    gpuDevice = getIntParam("--gpu");
	#endif
	} 
        void run()
        {
		main2();
        }

	int saveMat( const string& filename, const Mat& M)
	{
	    if (M.empty()){
	       return 0;
	    }
	    ofstream out(filename.c_str(), ios::out|ios::binary);
	    if (!out)
	       return 0;

	    int cols = M.cols;
	    int rows = M.rows;
	    int chan = M.channels();
	    int eSiz = (M.dataend-M.datastart)/(cols*rows*chan);

	    // Write header
	    out.write((char*)&cols,sizeof(cols));
	    out.write((char*)&rows,sizeof(rows));
	    out.write((char*)&chan,sizeof(chan));
	    out.write((char*)&eSiz,sizeof(eSiz));

	    // Write data.
	    if (M.isContinuous()){
	       out.write((char *)M.data,cols*rows*chan*eSiz);
	    }
	    else{
	       return 0;
	    }
	    out.close();
	    return 1;
	}

	int readMat( const string& filename, Mat& M)
	{
	    ifstream in(filename.c_str(), ios::in|ios::binary);
	    if (!in){
	       //M = NULL_MATRIX;
	       return 0;
	    }
	    int cols;
	    int rows;
	    int chan;
	    int eSiz;

	    // Read header
	    in.read((char*)&cols,sizeof(cols));
	    in.read((char*)&rows,sizeof(rows));
	    in.read((char*)&chan,sizeof(chan));
	    in.read((char*)&eSiz,sizeof(eSiz));

	    // Determine type of the matrix 
	    int type = 0;
	    switch (eSiz){
	    case sizeof(char):
		 type = CV_8UC(chan);
		 break;
	    case sizeof(float):
		 type = CV_32FC(chan);
		 break;
	    case sizeof(double):
		 type = CV_64FC(chan);
		 break;
	    }

	    // Alocate Matrix.
	    M = Mat(rows,cols,type,Scalar(1));  

	    // Read data.
	    if (M.isContinuous()){
	       in.read((char *)M.data,cols*rows*chan*eSiz);
	    }
	    else{
	       return 0;
	    }
	    in.close();
	    return 1;
	}

	void xmipp2Opencv(const MultidimArray<double> &xmippArray, Mat &opencvMat)
	{
	   int h=YSIZE(xmippArray);
	   int w=XSIZE(xmippArray);
	   opencvMat=cv::Mat::zeros(h, w,CV_32FC1);
	   FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
	      opencvMat.at<float>(i,j)=DIRECT_A2D_ELEM(xmippArray,i,j);      
	}

	void opencv2Xmipp(const Mat &opencvMat, MultidimArray<double> &xmippArray)
	{
	   int h=opencvMat.rows;
	   int w=opencvMat.cols;
	   xmippArray.initZeros(h, w);
	   FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
	      DIRECT_A2D_ELEM(xmippArray,i,j)=opencvMat.at<float>(i,j);      
	}

	void convert2Uint8(Mat opencvDoubleMat, Mat &opencvUintMat)
	{
	   Point minLoc,maxLoc;
	   double min,max;
	   Mat doubleMat;
	   cv::minMaxLoc(opencvDoubleMat,&min,&max,&minLoc,&maxLoc,noArray());
	   opencvDoubleMat = (opencvDoubleMat-min)*(255.0/(max-min));
	   opencvDoubleMat.convertTo(opencvUintMat, CV_8U);
	   //opencvDoubleMat.convertTo(opencvUintMat, CV_8U, 255.0/(max - min), -min * 255.0/(max - min));   
	}

	void computeAvg(const FileName &movieFile, int begin, int end, MultidimArray<double> &avgImg)
	{

	   int N=(end-begin)+1;
	   ImageGeneric movieStack;
	   MultidimArray<double> imgNormal, sumImg;

	   movieStack.readMapped(movieFile,begin);
	   movieStack().getImage(imgNormal);   
	   sumImg=imgNormal;
	     
	   for (int i=begin;i<end;i++)
	   {
	      movieStack.readMapped(movieFile,i+1);
	      movieStack().getImage(imgNormal);
	      sumImg+=imgNormal;            
	   }
	   avgImg=sumImg/double(N);
	}

	int main2()
	{
	    // XMIPP structures are defined here
	    MultidimArray<double> preImg, avgCurr, avgStep, mappedImg;
	    ImageGeneric movieStack, movieStackNormalize;
	    Image<double> II;
	    FileName movieFile = fname, normalizeFile;
	    ArrayDim aDim;
	    
	    // String for generating the file names
	    std::string theNumberString;
	    std::ostringstream ostr;
	    string fn;
	    ofstream fx, fy;
	    
	    // For measuring times
	    clock_t tStart, tStart2;    

	    #ifdef GPU
	    // Matrix that we required in GPU part
	    GpuMat d_flowx, d_flowy, d_dest;
	    GpuMat d_avgcurr, d_preimg, d_mapx, d_mapy;
	    #endif

	    // Matrix required by Opencv
	    Mat flowx, flowy, mapx, mapy, flow,dest; 
	    Mat avgcurr, avgstep, preimg, preimg8, avgcurr8;
	    Mat planes[] = {flowx, flowy};    
	    
	    // Object for optical flow 
	    
	    int imagenum, cnt=2, div=0, stackcnt=0;
	    int h, w, idx, levelNum, levelCounter=1;
	    
	    float* row_ptrx, * row_ptry;
	    size_t elem_step;     

	    normalizeFile = "normalized-"+fname;
	    movieStack.read(movieFile,HEADER);
	    movieStack.getDimensions(aDim);
	    imagenum=aDim.ndim;
	    h=aDim.ydim;
	    w=aDim.xdim;
	    levelNum=sqrt(double(imagenum));    

	 #ifdef GPU
	    FarnebackOpticalFlow d_calc;
	    gpu::setDevice(0); 

	  
	    // Initialize the parameters for optical flow structure
	    d_calc.numLevels=6;
	    d_calc.pyrScale=0.5;
	    d_calc.fastPyramids=false;
	    d_calc.winSize=winSize;
	    d_calc.numIters=1;
	    d_calc.polyN=5;
	    d_calc.polySigma=1.1;
	    d_calc.flags=0;
	 #endif  
	    // Initialize variables with zero
	    mapx=cv::Mat::zeros(h, w,CV_32FC1);
	    mapy=cv::Mat::zeros(h, w,CV_32FC1);
	    elem_step=mapx.step / sizeof(float);
		       
	    tStart2=clock();
	    
	    // Compute the average of the whole stack
	    computeAvg(movieFile, 1, imagenum, avgCurr);
	    II()=avgCurr;
	    //TODO: Is this debuging info?
            //Can we remove it
            II.write("avgimage.mrc");           
	  
	    while (div!=1)
	    {
	       div = int(imagenum/cnt);            
	       avgStep.initZeros(h, w);      
	       cout<<"Level "<<levelCounter<<"/"<<levelNum<<" of the pyramid is under processing"<<std::endl;
	       tStart = clock();       
	       idx = 0;
	       
	       if (div==1)
		  cnt=imagenum;
	
	       for (int i=0;i<cnt;i++)
	       {
	       
		  if (div==1)
		  {
		     movieStack.readMapped(movieFile,i+1);
		     movieStack().getImage(preImg);	  
		  }
		  else
		  {
		     if (i==cnt-1)
		        computeAvg(movieFile, (cnt-1)*div+1, imagenum, preImg);	     
		     else
		        computeAvg(movieFile, i*div+1, (i+1)*div, preImg);
		  }          
		  
		  xmipp2Opencv(avgCurr, avgcurr);
		  xmipp2Opencv(preImg, preimg);
		  convert2Uint8(avgcurr,avgcurr8);
		  convert2Uint8(preimg,preimg8);
	#ifdef GPU	  
		     d_avgcurr.upload(avgcurr8);
		     d_preimg.upload(preimg8);	  
	
		     d_calc(d_avgcurr, d_preimg, d_flowx, d_flowy);
		     d_flowx.download(flowx);
		     d_flowy.download(flowy);
		     	  
		     d_avgcurr.release();
		     d_preimg.release();
		     d_flowx.release();
		     d_flowy.release();
	#else
		     calcOpticalFlowFarneback(avgcurr8, preimg8, flow, 0.5, 6, winSize, 1, 5, 1.1, 0);
		     split(flow, planes);
		     flowx = planes[0]; flowy = planes[1];	
	#endif	  
		  flowx.convertTo(mapx, CV_32FC1);
		  flowy.convertTo(mapy, CV_32FC1);	  		  	  

		  for( int row = 0; row < mapx.rows; row++ )
		     for( int col = 0; col < mapx.cols; col++ )
		     {
			mapx.at<float>(row,col) += (float)col;
			mapy.at<float>(row,col) += (float)row;
		
		     }
	#ifdef GPU

		  {
		     d_mapx.upload(mapx);	     
		     d_mapy.upload(mapy);
		     d_preimg.upload(preimg);
		     gpu::remap(d_preimg,d_dest,d_mapx,d_mapy,INTER_CUBIC);
		     d_dest.download(dest);
		  
		     d_dest.release();
		     d_preimg.release();
		     d_mapx.release();
		     d_mapy.release();
		  }
	#else	  	  
		     remap(preimg, dest, mapx, mapy, INTER_CUBIC);	  
	#endif
		  	  
		  opencv2Xmipp(dest, mappedImg);	  
		  avgStep+=mappedImg;  
	       }
	       avgCurr =  avgStep/cnt;    
	       cout<<"Processing level "<<levelCounter<<"/"<<levelNum<<" has been finished"<<std::endl;
	       printf("Processing time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	       cnt = cnt * 2;  
	       levelCounter++;
	    }
	    II()=avgCurr;
	    II.write(foname);
	    printf("Total Processing time: %.2fs\n", (double)(clock() - tStart2)/CLOCKS_PER_SEC);    
	    return 0;
	}
};

int main(int argc, char *argv[])
{
    ProgOpticalAligment prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
