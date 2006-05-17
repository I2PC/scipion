//IMPORTANT: Coord. Savings are in System:
//                Origin: Image Center
//                Axes Negatively Oriented

#include "exp_shift_computation.hh"




#define BORDER_POINT 10000

//////////////////////////// AUXILIAR FUNCTIONS
inline float min( float f1, float f2 )
{
	return f1<f2 ? f1 : f2;
}




// Comparison function for qsort of Lat. Pts 
// in increasing x values.
// Needed for Computation of DELAUNAY Triang.
/*int XYZCompare(const void *v1,const void *v2)
{
      XYZ *p1,*p2;
      p1 = (XYZ *) v1;
      p2 = (XYZ *) v2;

      return (p1->x - p2->x);
        
}*/

int LatPtCompare(const void *v1,const void *v2)
{
      LatPoint *p1,*p2;
      p1 = (LatPoint *) v1;
      p2 = (LatPoint *) v2;

      return (int)(p1->x - p2->x);
        
}

void Usage(char *argv[]);

//////////////////////////////////////////////////////////
////////////////////////// MAIN 

int main(int argc,char **argv)
{
  using namespace std;
  /////////////////////// Command Line Parameters
  // I/O Files
  string filein;
  string fileout;
  // THRESHOLD OF CROSS-CORRELATION PEAK HEIGHT CALCULATED AS;
  //     max of croscorrelation * FACTOR (READ FROM UNIT 5)
   double cc_peak_factor;
   
  // Get command line parameters ------------------------------------------
   try {
      filein  = get_param(argc, argv, "-i");
      fileout = get_param(argc, argv, "-o");
      cc_peak_factor =  AtoF(get_param(argc,argv,"-cc_peak_factor","0.0"));
      
  } catch (Xmipp_error XE) {Usage(argv);}
  

	//Peaks Coord.
	 /** Experimental Lattice and Peaks Coordinates */
         CCLattice_IO ExpLat;
         /**  Experimental Displacements */
         vector <LatPoint> INCR_coord;
  
	cout<<"Reading"<<endl;
	//Read .cor file
	ExpLat.read(filein.c_str());
    
	cout<<"PeaksCorrespondance"<<endl;
	//Uniform sistematic sampling
        PeaksCorresp(ExpLat,INCR_coord,cc_peak_factor);
 
	cout<<"BendingInterp"<<endl;
	//Displacement Interpolation
        BendingInterp(INCR_coord);

	cout<<"Savings"<<endl;
    //Save Lat. displacements
    SaveLatBending(fileout.c_str(),ExpLat,INCR_coord);

	return 0;

}

void Usage (char *argv[]) {
    cout << "Purpose:\n"
         << "    Computes 2D experimental distorsion at unit cell center\n"
         << "    The input is a .cor (MRC) file\n"
         << "Usage:" << argv[0] <<" -i filein -o fileout " 
	 << " -cc_peak_factor cc_peak_factor " << endl << endl
         << "\t-i               :  Input Correlation file" << endl
         << "\t-o               :  Output Distortion file" << endl
         << "\t-cc_peak_factor  : crosscorrelation thershold (0-1)" << endl
	;    
    exit(1);

}
/////////////////////////////////////////////////////
//////////////////////////// PEAKS CORRESPONDANCE

////Peaks Correspondance. 
void PeaksCorresp(CCLattice_IO & ExpLat,vector <LatPoint> &INCR_coord, double cc_peak_factor)
{
   double auxX,auxY;
   LatPoint valCoord;
    //Threshold over CCvalue
   double Th=cc_peak_factor* ExpLat.cc_max;
   
   
   for (int kk = 0; kk < ExpLat.MRC_Xcoord.size(); kk++)
   {
   	
       //Compute ideal position in regular grid	
       valCoord.i=ExpLat.MRC_Xindex[kk];
       valCoord.j=ExpLat.MRC_Yindex[kk];
       valCoord.x=ExpLat.O[0] + ExpLat.MRC_Xindex[kk]*ExpLat.a(0) + ExpLat.MRC_Yindex[kk]*ExpLat.b(0);
       valCoord.y=ExpLat.O[1] + ExpLat.MRC_Xindex[kk]*ExpLat.a(1) + ExpLat.MRC_Yindex[kk]*ExpLat.b(1);
       //Deviation of experimental position
       valCoord.Incrx=ExpLat.MRC_Xcoord[kk] - valCoord.x;
       valCoord.Incry=ExpLat.MRC_Ycoord[kk] - valCoord.y;
       valCoord.Interp=(ExpLat.MRC_CCcoheficiente[kk]<Th);
       INCR_coord.push_back(valCoord);
      }

   #undef DEBUG
 //  #define DEBUG
   #ifdef DEBUG
   for (int kk = 0; kk < INCR_coord.size(); kk++)
	 cout <<  INCR_coord[kk].x << " " << INCR_coord[kk].y << " " 
	     <<  INCR_coord[kk].Incrx << " " << INCR_coord[kk].Incry << endl;   

   #endif
   #undef DEBUG
}


//////////////////////////////////////////////////////////////
///////////////////////////// INTERPOLATION
void BendingInterp(vector <LatPoint> &INCR_coord)
{

	int k,n;
	int N=INCR_coord.size();
	int Tind;
	//Del. Triang.
	vector <ITRIANGLE> LatTri;
        //Interpolation Parameters
        float w[3],del;
	float x[3],y[3],Incrx[3],Incry[3];
	float xi,yi;
        LatPoint TargetPt;
	
        //Lat. Triang.
	LatTriang(INCR_coord, LatTri);

        for(k=0;k<N;k++){
          if(INCR_coord[k].Interp){
             TargetPt.x=INCR_coord[k].x;
             TargetPt.y=INCR_coord[k].y;
             TargetPt.Incrx=0;
             TargetPt.Incry=0;
             LinInterp(INCR_coord,TargetPt,LatTri);
             INCR_coord[k]=TargetPt;
            }
        }
}

// Displacement Interpolation from Triangulation of Irregular Grid
void LinInterp(vector <LatPoint> & INCR_coord,LatPoint &TargetPt, vector <ITRIANGLE> &  LatTri)
{
	//Interpolation Parameters
        int Tind;  
	int n;
        float w[3],del;
	float x[3],y[3],Incrx[3],Incry[3];
	float xi,yi;
	  
	       
	           
        //#undef TIMES
        //#define TIMES
        #ifdef TIMES
	struct tms before,after;
	times(&before);
	#endif
	
	        
                //find nearest Triangle
		Tind=FindNearest(INCR_coord,TargetPt, LatTri);
		#ifdef TIMES
                times(&after);
		cout << "Triang. search time " << after.tms_utime-before.tms_utime << endl;
		#endif
		
		if(Tind>0)
		{
		//Lattice Position of Target
		 xi=TargetPt.x;
                 yi=TargetPt.y;
		//Lattice Position of Interpolant nearest Pts
	         x[0]=INCR_coord[LatTri[Tind].p1].x;
                 y[0]=INCR_coord[LatTri[Tind].p1].y;
		 Incrx[0]=INCR_coord[LatTri[Tind].p1].Incrx;
		 Incry[0]=INCR_coord[LatTri[Tind].p1].Incry;
	         x[1]=INCR_coord[LatTri[Tind].p2].x;
                 y[1]=INCR_coord[LatTri[Tind].p2].y;
		 Incrx[1]=INCR_coord[LatTri[Tind].p2].Incrx;
		 Incry[1]=INCR_coord[LatTri[Tind].p2].Incry;
	         x[2]=INCR_coord[LatTri[Tind].p3].x;
                 y[2]=INCR_coord[LatTri[Tind].p3].y;
		 Incrx[2]=INCR_coord[LatTri[Tind].p3].Incrx;
		 Incry[2]=INCR_coord[LatTri[Tind].p3].Incry;
	
	         
		//Linear Interpolation (from MatLab griddata)
		//Barycentric Coord. 
		//IMPORTANT: Need that (x[k],y[k]) define a triangle (del!=0)
	       del=(x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
	        w[2]=((x[0]-xi)*(y[1]-yi)-(x[1]-xi)*(y[0]-yi))/del;
		  w[1]=((x[2]-xi)*(y[0]-yi)-(x[0]-xi)*(y[2]-yi))/del;
		  w[0]=((x[1]-xi)*(y[2]-yi)-(x[2]-xi)*(y[1]-yi))/del;
	
		//Displacement Interpolation

		 for(n=0;n<3;n++)
		 {
		  TargetPt.Incrx=TargetPt.Incrx+Incrx[n]*w[n];
		  TargetPt.Incry=TargetPt.Incry+Incry[n]*w[n];
		 }
		}
		else
		{
		 TargetPt.Incrx=BORDER_POINT;
		 TargetPt.Incry=BORDER_POINT;
		}
}

//////////////////////////////////////////////// TRIANGULATION HANDLINNG ////////////////////////////////////

/////////////Delaunay Triangulation of Lat. Pts
void LatTriang(vector <LatPoint> &INCR_coord, vector <ITRIANGLE> & LatTri)
{
   int i;
   int N=INCR_coord.size();
   int ntri = 0;
   int incr=1;
   ITRIANGLE *v=NULL;
   XYZ *p = NULL;
   LatPoint * pLat=NULL;
   p=(XYZ *) malloc((N+3)*sizeof(XYZ));
   v = (ITRIANGLE *) malloc(3*N*sizeof(ITRIANGLE));
   pLat=(LatPoint *)malloc(N*sizeof(LatPoint));

  

   //Sort vertixes in x increasing order
   for(i=0;i<N;i++)
   {
	   pLat[i]=INCR_coord[i];
   }
   qsort(pLat,N,sizeof(LatPoint),LatPtCompare);
   //Copy sorted values to INCR_coord vector
   for(i=0;i<N;i++)
   {
	   INCR_coord[i]=pLat[i];
   }
  

   //Set Lat. Pts Coord.
   for(i=0;i<N;i++)
   {
        p[i].x=INCR_coord[i].x;
	p[i].y=INCR_coord[i].y;
	p[i].z=0;
   }
  
   Triangulate(N,p,v,&ntri);
  

   //Copy Triangulation to vector format
   for(i=0;i<ntri;i++)
   {
	   LatTri.push_back(v[i]);
   }
   

   //free memory
   free(pLat);
   free(v);
   free(p);

}

/////////////////// Computation of Nearest Triangl. to TargetPt
// Returns index of triangle in vector LatTri
int FindNearest(vector <LatPoint> &INCR_coord, LatPoint &TargetPt, vector <ITRIANGLE> &  LatTri)
{

	int k,Vind,t;
	XYZ Pi,P1,P2,P3;
	float x1,x2,x3,y1,y2,y3;
	int N=LatTri.size();
    
    //Target point position
    Pi.x=TargetPt.x;
    Pi.y=TargetPt.y;


	//Nearest Triangle Computation
	for(k=0;k<N;k++)
	{
        
    	//Triangle Vertex position
	  Vind=LatTri[k].p1;
          P1.x=INCR_coord[Vind].x;
	  P1.y=INCR_coord[Vind].y;
	  Vind=LatTri[k].p2;
          P2.x=INCR_coord[Vind].x;
	  P2.y=INCR_coord[Vind].y;
	  Vind=LatTri[k].p3;
          P3.x=INCR_coord[Vind].x;
	  P3.y=INCR_coord[Vind].y;

        //Inside Condition
		x1 = P1.x-Pi.x; y1 = P1.y-Pi.y; 
        x2 = P2.x-Pi.x; y2 = P2.y-Pi.y; 
        x3 = P3.x-Pi.x; y3 = P3.y-Pi.y; 
        t =  (x1*y2 > x2*y1) + (x2*y3 > x3*y2) + (x3*y1 > x1*y3);
       
       // printf("Inside Condition of Triangle %d is %d:\n",k,t);

        if((t==3) || (t==0)) return k;
	}

   return -1;
}

	
//////////////////////////////////////////////////////////////
//////////////////////////////  I/O FUNCTIONS

// Read *.cor file. Currently using CCLattice_IO class ReadMRCCord function

/* void ReadMRCCord(const char * filename, LatParam &L,vector <float> &MRC_Xcoord,vector <float> &MRC_Ycoord)
{
    
	FILE * fid;
        float valXCoord,valYCoord,val;
        char s[100];

	fid=fopen(filename,"r");

	if(fid == NULL)
		exit(1);

        //Skip first 5 lines
        for(int k=0;k<5;k++)
             fgets(s, sizeof s, fid);

        
	//Lattice Param.
	fscanf(fid,"%d\t%d",&L.dim[0],&L.dim[1]);
	fscanf(fid,"%f\t%f",&valXCoord,&valYCoord);
	L.O.push_back(valXCoord);
    L.O.push_back(valYCoord);
    fscanf(fid,"%f\t%f",&valXCoord,&valYCoord);
	L.a.push_back(valXCoord);
    L.a.push_back(valYCoord);
	fscanf(fid,"%f\t%f",&valXCoord,&valYCoord);
    L.b.push_back(valXCoord);
    L.b.push_back(valYCoord);
    fscanf(fid,"%d\t%d",&L.Na[0],&L.Na[1]);
	fscanf(fid,"%d\t%d",&L.Nb[0],&L.Nb[1]);

        //Skip Misterious number
        fscanf(fid,"%f\n",&valXCoord);


   	//Peaks coordinates
	while(!feof(fid))
	{
      fscanf(fid,"%f\t%f\t%f",&valXCoord,&valYCoord,&val);
	  if(valXCoord*valYCoord*val)
	  {   MRC_Xcoord.push_back(valXCoord);	
		  MRC_Ycoord.push_back(valYCoord);	
	  }
	}


    fclose(fid);

}*/

// Save Lat. displacements in format
//  i j x y Incrx Incry
// with OUR coord system set to:
//                Origin: Image Center
//                Axes Negatively Oriented
void  SaveLatBending(const char * fileout,CCLattice_IO & ExpLat,vector <LatPoint> & INCR_coord)
{

	FILE * fid;
        int N=INCR_coord.size();
	int k;

	fid=fopen(fileout,"w");

	if(fid == NULL)
			return;

	//Lattice Info
        fprintf(fid,"; Dimensions: dimx dimy Ox Oy \n");
	fprintf(fid,";%d %d %f %f \n", ExpLat.dim[0],ExpLat.dim[1],ExpLat.O[0],ExpLat.O[1]);
	fprintf(fid,"; Lat. Vectors: ax ay bx by \n");
	fprintf(fid,";%f  %f  %f  %f \n",ExpLat.a(0),ExpLat.a(1),ExpLat.b(0),ExpLat.b(1));

	fprintf(fid,"; NumFila NumElem  h   k  X   Y IncrementoX  IncrementoY\n");
	//Correspondance Coord.
    for(k=0;k<N;k++)
	{   
        fprintf(fid,"%7d%3d",k+1,6);
	// MRC coord: 
	// Origin: Left bottom corner of image
	// Axes Positively Oriented
	// OUR coord:
	// Origin: Image Center
	// Axes Negatively Oriented

        fprintf(fid,"%4d%4d%9.2f%9.2f%9.3f%9.3f\n",INCR_coord[k].i,INCR_coord[k].j,INCR_coord[k].x,INCR_coord[k].y,INCR_coord[k].Incrx,INCR_coord[k].Incry);
	}
     
	
    fclose(fid);
}
