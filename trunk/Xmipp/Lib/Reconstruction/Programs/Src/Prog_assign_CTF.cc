/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.uam.es)
 *              Carlos Oscar Sánchez Sorzano
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

#include "../Prog_assign_CTF.hh"

#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippMicrograph.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippHeader.hh>

// Number of digits for file names
# define NUMBER_OF_DIGITS    5  

/****************************************************************************

   NAME:        RegionNumber
   
   DESCRIPTION:	This function gives the region number assigned to a particle's
   		coordinates. 
		It's assummed that Regions start in 1 and advance horizontally.
		-------------------
		| 1 | 2 | 3 | 4 | |
		-------------------
		| 5 | 6 | 7 | 8 | |				
		-------------------
		|   |   |   |   | |
		-------------------
		
		Coordinates that are outside the Regions are assigned the nearest
		one. 
			
   INPUT:       X, Y - Coordinates of the particle
   		N_horizontal, N_vertical - Size of Regions horizontally and 
				 							vertically
		Xdiv, Ydiv - Number of Regions in horizontal and vertical
                directions
				
   OUTPUT:      The region number. 
   
****************************************************************************/			 
int RegionNumber(int X,int Y,int N_horizontal,int N_vertical, int Xdiv,int Ydiv)
{
  int Vertical  =((Y/N_vertical  +1) >  Ydiv) ? (Y/N_vertical)   : (Y/N_vertical+1);
  int Horizontal=((X/N_horizontal+1) >  Xdiv) ? (X/N_horizontal) : (X/N_horizontal+1);
  return (Vertical-1)*Xdiv+Horizontal;
}


/****************************************************************************

   NAME:        CenterOfRegion
   
   DESCRIPTION:	Given a region number, and a direction, this function returns
   		the pixel of its center.
		Also needed parameters are the size of the region in the 
		considered direction, and the number of horizontal divisions.				
			
   INPUT:       RegionNumber - an integer
  		direction - a char that could be 'x' or 'y' 
		size - size of region in the selected direction
		Xdiv - number of horizontal divisions
								
   OUTPUT:      The pixel where the center is, in the selected direction 
   
****************************************************************************/			 
int CenterOfRegion(int RegionNumber,char direction,int size,int Xdiv)
{
	if(direction=='x') 
	{
		// RegionNumber % Xdiv is done to calculate the column
		return size/2+((RegionNumber-1) % Xdiv)*size;
	}
	else if(direction=='y')
	{
		// RegionNumber / Xdiv is done to calculate the row
		return size/2+((RegionNumber-1) / Xdiv)*size;
	}
}

/****************************************************************************

   NAME:        DetermineRegionsToInterpolate
   
   DESCRIPTION:	This function gives the Regions of the micrograph to interpolate,
   				given a particle's position. It takes into account three cases:
  			-If the particle's center is between four region's centers. 
			-If it is only between two region's centers.
   			-If the particle in near the corners of the micrograph, only the region
   			to which it belongs is needed.
			
   INPUT:        X, Y - Coordinates of the particle
   				 N_horizontal, N_vertical - Size of Regions horizontally and 
				 							vertically
				Xdiv, Ydiv - Number of Regions in horizontal and vertical directions
				
   OUTPUT:      The region numbers are stored in Regions argument. A zero
   				value means no region.
   
****************************************************************************/			 
void  DetermineRegionsToInterpolate(int X,int Y,
	int N_horizontal,int N_vertical, 
	int Xdiv,int Ydiv,matrix1D<int> &Regions)
{
	Regions.resize(4);
	Regions.init_zeros();
	
	// First determine the region to which the particle belongs
	int N=RegionNumber(X,Y,N_horizontal,N_vertical,Xdiv,Ydiv);
	// The compute the coordinates of the center of that region
	int Xc=CenterOfRegion(N,'x',N_horizontal,Xdiv);
	int Yc=CenterOfRegion(N,'y',N_vertical,Xdiv);
	// Add Regions according to the position of the particle's center
	// respect to the region's center.
	if(X>=Xc) 
	{		
		// If the particle's center is greater than the last region's center
		if(X>(N_horizontal/2+(Xdiv-1)*N_horizontal))
		{
			// Add only the last region, wich is actually the particle's region
			Regions(0)=N;
		}
		else
		{
			// Add the actual region and the following one
			Regions(0)=N;
			Regions(1)=N+1;
		}
	}
	// IF The particle's center is lower than the region's center
	else 
	{
		// And lower than the first region's center
		if(X<(N_horizontal/2))
		{
			// Add only the first region, wich is actually the particle's region
			Regions(0)=N;
		}
		else
		{
			// Add the actual region and the preceding one
			Regions(0)=N-1;
			Regions(1)=N;
		}
	}
	// If the particle's center is greater than the region's one 
	if(Y>=Yc) 
	{
		// And greater than the last region's center
		if(Y>(N_vertical/2+(Ydiv-1)*N_vertical))
		{
			// Do nothing, because this region has been added before 
		}
		else
		{
			// Add the Regions of the next row, provided the corresponding
			// ones in the actual row where added. That is:
			if(Regions(0)!=0) Regions(2)=Regions(0)+Xdiv;
			if(Regions(1)!=0) Regions(3)=Regions(1)+Xdiv;
		}
	}
	// IF The particle's center is lower than the region's center
	else 
	{
		// And lower than the first region's center
		if(Y<(N_vertical/2))
		{
			// Do nothing, because this region has been added before 
		}
		else
		{
			// Add the Regions of the preceding row, provided the corresponding
			// ones in the actual row where added. That is:
			if(Regions(0)!=0) Regions(2)=Regions(0)-Xdiv;
			if(Regions(1)!=0) Regions(3)=Regions(1)-Xdiv;
		}
	}
}

/****************************************************************************

   NAME:        InterpolateRegionsCTF
   
   DESCRIPTION:	Does bilinear or linear interpolation between the CTF files
   				for the Regions considered.
  			If the particle's center is between four region's centers, then a
   			bilinear interpolation is performed. If it is only between two region's
   			centers, a linear interpolation is donde between the two images.
   			If the particle in near the corners of the micrograph, only the region
   			to which it belongs is needed.
			
   INPUT:        X, Y coordinates of the particle's center
   				 N_horizontal, N_vertical - Size of Regions horizontally and 
				 							vertically
				Xdiv, Ydiv - Number of Regions in horizontal and vertical directions
 				 Regions to perform the interpolation - A zero value in this
				 matrix means no region
				 fn_CTFfile - root name for files containing the CTF of the
				 		 	selected Regions
				 Result - The resulting interpolated matrix
								
   OUTPUT:      The interpolated matrix is stored in Result 
   
****************************************************************************/			 
void InterpolateRegionsCTF(int X,int Y, int N_horizontal,int N_vertical, 
   int Xdiv,int Ydiv, matrix1D<int> &Regions,
   FileName fn_CTFfile, matrix2D<complex <double> > &Result)
{
	Result.resize(N_vertical,N_horizontal);
	Result.init_zeros();
	FourierImageXmipp CTFfile;
	FileName fn_1,fn_2,fn_3,fn_4;
    double w1,w2,w3,w4; // weights for interpolation

	// Count the number of zeros in Regions
	int Zeros=0;
	for(int i=0;i<Regions.get_dim();i++)
		if(Regions(i)==0) Zeros++;
	switch(Zeros)
	{
		case 3:  // We don't need to interpolate
				// Read the corresponding CTF model and return it
				fn_CTFfile.compose(fn_CTFfile.get_root().c_str(),Regions(0),"fft");
				CTFfile.read(fn_CTFfile);
				Result=CTFfile(); 
				break;
				
		case 2:	// There are two Regions. A linear interpolation is needed
				// If Regions(1) equals zero, then the Regions are vertical				
				if(Regions(1)==0)
				{
					int y1=CenterOfRegion(Regions(0),'y',N_vertical,Xdiv);
					int y2=CenterOfRegion(Regions(2),'y',N_vertical,Xdiv);
					w2=(double)(Y-y1)/(double)(y2-y1);
					w1=1-w2;

					fn_1.compose(fn_CTFfile.get_root().c_str(),Regions(0),"fft");
					fn_2.compose(fn_CTFfile.get_root().c_str(),Regions(2),"fft");
				}
				// Otherwise they are horizontal
				else
				{
					int x1=CenterOfRegion(Regions(0),'x',N_horizontal,Xdiv);
					int x2=CenterOfRegion(Regions(1),'x',N_horizontal,Xdiv);
					w2=(double)(X-x1)/(double)(x2-x1);
					w1=1-w2;
					fn_1.compose(fn_CTFfile.get_root().c_str(),Regions(0),"fft");
					fn_2.compose(fn_CTFfile.get_root().c_str(),Regions(1),"fft");
				}

				// Interpolate the two files.
				CTFfile.read(fn_1);
				CTFfile.read(fn_2);
				FOR_ALL_ELEMENTS_IN_MATRIX2D(CTFfile())
					Result(i,j)+=w2*CTFfile(i,j);
				break;
		case 0: // A bilinear interpolation is needed
				/* Calcuate centers of Regions. It's assumed that the Regions
				 are
				 ---------
				 | 0 | 1 |
				 ---------
				 | 2 | 3 |
				 ---------
				 so:
				  the X-center of Region(0) is the same of Region(2)
				  the Y-center of Region(0) is the same of Region(1)
				*/
				int y1=CenterOfRegion(Regions(0),'y',N_vertical,Xdiv);
				int y2=CenterOfRegion(Regions(2),'y',N_vertical,Xdiv);
				int x1=CenterOfRegion(Regions(0),'x',N_horizontal,Xdiv);
				int x2=CenterOfRegion(Regions(1),'x',N_horizontal,Xdiv);
				// Calculate weights
				double alpha=(double)(Y-y1)/(double)(y2-y1);
				double beta =(double)(X-x1)/(double)(x2-x1);				
				w1=(1-alpha)*(1-beta);
				w2=beta*(1-alpha);
				w3=alpha*(1-beta);
				w4=alpha*beta;
				fn_1.compose(fn_CTFfile.get_root().c_str(),Regions(0),"fft");
				fn_2.compose(fn_CTFfile.get_root().c_str(),Regions(1),"fft");
				fn_3.compose(fn_CTFfile.get_root().c_str(),Regions(2),"fft");
				fn_4.compose(fn_CTFfile.get_root().c_str(),Regions(3),"fft");
				// Interpolate
				CTFfile.read(fn_1); Result+=w1*CTFfile();
				CTFfile.read(fn_2); Result+=w2*CTFfile();
				CTFfile.read(fn_3); Result+=w3*CTFfile();
				CTFfile.read(fn_4); Result+=w4*CTFfile();
				
				break;
	}
}

template <class T>
void FFT_idx2digfreq(matrix2D<T> &M, matrix1D<int> &idx, matrix1D<double> &freq)
{
	matrix1D<int> size(2);
	size(1)=M.ColNo(); size(0)=M.RowNo();
    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
	{
    	if(idx(i)<size(i)/2)
			freq(i)=-size(i)+idx(i);
		else
			freq(i)=idx(i);
		freq(i)/=size(i);
	}
}

/* Read parameters ========================================================= */
void Prog_assign_CTF_prm::read(const FileName &fn_prm, bool do_not_read_files) _THROW {
   // Read parameters for adjust CTF from input file
   adjust_CTF_prm.read(fn_prm);
   
   // Read specific parameters for this program from input file
   FILE *fh_param;
   if ((fh_param = fopen(fn_prm.c_str(), "r")) == NULL)
      REPORT_ERROR(1,(string)"assign_CTF: There is a problem "
            "opening the file "+fn_prm); 
                                      
   reversed          =check_param(fh_param,"reverse endian");
   downsampling      =AtoI(get_param(fh_param,"downsampling",0,"1"));
   N_horizontal      =AtoI(get_param(fh_param,"N_horizontal",0));
   N_vertical        =AtoI(get_param(fh_param,"N_vertical",0));
   particle_horizontal=AtoI(get_param(fh_param,"particle_horizontal",0));
   particle_vertical =AtoI(get_param(fh_param,"particle_vertical",0));
   if (!do_not_read_files) {
      image_fn          =get_param(fh_param,"image",0);
      selfile_fn        =get_param(fh_param,"selfile",0);
      picked_fn         =get_param(fh_param,"picked",0);
      ARMAfile          =get_param(fh_param,"ARMAfile",0);
      CTFfile           =get_param(fh_param,"CTFfile",0);
   }

   only_interpolate     =check_param(fh_param,"only_interpolate");
   compute_at_particle  =check_param(fh_param,"compute_at_particle");
   ARMA_averaging       =check_param(fh_param,"ARMA_averaging");
   Periodogram_averaging=check_param(fh_param,"Periodogram_averaging");
   fclose(fh_param);

   // Read ARMA parameters from input file
   if (!Periodogram_averaging)
      ARMA_prm.read(fn_prm);
}

/* Write parameters ========================================================= */
void Prog_assign_CTF_prm::write(const FileName &fn_prm) _THROW {
   ofstream fh_param;
   fh_param.open(fn_prm.c_str(),ios::out);
   if (!fh_param)
      REPORT_ERROR(1,(string)"assign_CTF: There is a problem "
            "opening the file "+fn_prm+" for write");
   fh_param << "# Assign CTF parameters\n";
   fh_param << "image="                << image_fn             << endl
            << "N_horizontal="         << N_horizontal         << endl
            << "N_vertical="           << N_vertical           << endl
            << "particle_horizontal="  << particle_horizontal  << endl
            << "particle_vertical="    << particle_vertical    << endl
            << "selfile="              << selfile_fn           << endl
            << "picked="               << picked_fn            << endl
            << "ARMAfile="             << ARMAfile             << endl
            << "CTFfile="              << CTFfile              << endl
   ;
   if (only_interpolate) fh_param << "only_interpolate=yes\n";
   if (compute_at_particle) fh_param << "compute_at_particle=yes\n";
   if (ARMA_averaging) fh_param << "ARMA_averaging=yes\n";
   if (Periodogram_averaging) fh_param << "Periodogram_averaging=yes\n";
  
   fh_param << endl;
   fh_param.close();

   adjust_CTF_prm.write(fn_prm,FALSE);
   ARMA_prm.write(fn_prm,FALSE);
}

//#define DEBUG
/* Main ==================================================================== */
void Prog_assign_CTF_prm::process() {
   ifstream PosFile(picked_fn.c_str());  // File with picked coordinates
   SelFile  SF(selfile_fn);  // Selfile
   FileName fn_out=image_fn.remove_all_extensions();
   ofstream OutputFile_ctf((fn_out+"_ctf.sel").c_str());

   /*****************************************************************************
   	   Read input micrograph
   /*****************************************************************************/
   Micrograph M_in;
   M_in.open_micrograph(image_fn,reversed);
   int bits=M_in.depth();
   int Ydim, Xdim; // Micrograph dimensions
   M_in.size(Xdim, Ydim);

   if (!only_interpolate) {
	// Counter of divisions;
	int div_Number=0;
	int N=0;
	/*****************************************************************************
	 Divide the micrograph into smaller ones, as many as the user desires
	 Perform an spectral estimation of every division by means of an ARMA model
	*****************************************************************************/
        // Compute the number of divisions
	if (compute_at_particle) {
	   //check if sel file is empty
	   if(SF.LineNo()==0){
	      cerr << "Prog_assign_CTF_prm: sel file " << SF.name() 
	           << "is empty " << "I can not go any further sorry" 
		   << endl;
	      exit(1);
	    }
	      	    
           string line;
           PosFile.clear();       
	   PosFile.seekg(0, ios::beg); // rewing file
           while(getline(PosFile,line))
               if (line[0]!='#') N++;

           // check that the number of entries in the pos file is the right one
	   if(SF.LineNo()!=N){
	      cerr << "Prog_assign_CTF_prm: number of entries in "
	           << "pos file: "<< picked_fn.c_str() 
		   << "(equal to: " << N << ") "
		   << " and sel file " 
		   << SF.name() << "(equal to: " 
		   << SF.LineNo() << ") "
		   << "is different"
	           << "I can not go any further sorry" 
		   << endl;
	      exit(1);
	   }

           //find out the particle dimensions
           if (particle_horizontal<=0||particle_vertical<=0)
	      SF.ImgSize(particle_vertical,particle_horizontal);
	} else {
	   if (ARMA_averaging || Periodogram_averaging)
	      N=(FLOOR(Ydim/(N_vertical/2))-1)*(FLOOR(Xdim/(N_horizontal/2))-1);
	   else
	      N=FLOOR(Ydim/N_vertical)*FLOOR(Xdim/N_horizontal);
	}

	#ifdef DEBUG
	   cout << "N=" << N << endl;
	#endif
	div_Number=N;
	int ARMA_Number=(ARMA_averaging || Periodogram_averaging)? 1:N;
        matrix1D<double> X_zeros(ARMA_Number), Y_zeros(ARMA_Number);

        // Open the output file for the CTF zeros
        FileName fn_OutputZeros=fn_out.insert_before_extension("_zeros");
        ofstream OutputZeros;
        OutputZeros.open(fn_OutputZeros.c_str(), ios::out);
        if (!OutputZeros)
           REPORT_ERROR(1551,(string)"Assign_CTF: Cannot open "+fn_OutputZeros+" for output");
	OutputZeros << "# Piece_name first_zero_in_X first_zero_in_Y fitting_error\n";

        N=1;
        // Save the original sampling rate
	// The adjustment is done with one sampling rate related with the downsampling
	// while the CTF generation should be done without downsampling
	double original_sampling_rate=adjust_CTF_prm.Tm;
	adjust_CTF_prm.Tm*=downsampling;
	int i,j;i=j=0;
	float fi,fj;
	// just in case you wonder what is the clear for:
	/* When reading the file for the first time, you probably continue until
        an end of file condition occurs.  Before doing anything else with the
        ifstream, you will have to clear the error condition:
        */

        PosFile.clear();       
	PosFile.seekg(0, ios::beg); // Start of file
	SF.go_beginning();
      	FourierImageXmipp Filter_ARMA_avg;
	while (N<=div_Number) {
           if (compute_at_particle) {
              string line;
              getline(PosFile,line);
              while (line[0]=='#')
        	  getline(PosFile,line);	    
              sscanf(line.c_str(),"%f %f",&fj,&fi);
	      i = (int) fi;
	      j = (int) fj;
	      #ifdef DEBUG
        	 cout << "line" << line << endl;
		 cout << "read from file (j,i)= (" << j << "," << i << ")" <<endl;
		 cout << "Particle file name: " << SF.get_current_file() << endl;
	      #endif

	      //j,i are the window center, we need the top-left corner
	      j -= (int) (N_horizontal/2);
	      i -= (int) (N_vertical/2);
 	      if (i<0) i=0;
 	      if (j<0) j=0;
 	      if (i>Ydim-N_horizontal) i=Ydim-N_horizontal-1;
 	      if (j>Xdim-N_vertical)   j=Xdim-N_vertical-1;
	   } else {
	      int Xstep=N_horizontal, Ystep=N_vertical;
	      if (ARMA_averaging || Periodogram_averaging)
	         {Xstep/=2; Ystep/=2;}
	      if(N!=1) j+=Xstep;
	      if (j>Xdim) {
		 j=0;
		 i+=Ystep;
	      }
	      if (i>Ydim)
	         REPORT_ERROR(1,"Prog_assign_CTF_prm: window is outside the micrograph");
	   }

           #define ijDEBUG
	   #ifdef ijDEBUG
	      cout << "Block " << N <<" (j,i)= (" << j << "," << i << ")" <<endl;
	   #endif
      	   // test if i,j is inside the micrograph
           // this should never happend
      	   if (i+N_vertical>Ydim || j+N_horizontal>Xdim) continue;

      	   // The downsampling is done by taking 1 out of 2,3, ... pixels
	   // in the extracted region
           // Extract the division image where to perform the ARMA model
    	   ImageXmipp ImgDivision(FLOOR(N_vertical/downsampling),
      	             	      	  FLOOR(N_horizontal/downsampling));

           // image to store the filter     
           cerr << "Reading micrograph region number " << N
	        << " out of " << div_Number << endl;
    	   for (int k=0; k<YSIZE(ImgDivision()); k++)
               for (int l=0; l<XSIZE(ImgDivision()); l++)
		   ImgDivision(k,l)=
      	              M_in(j+downsampling*l,i+downsampling*k);
           ImgDivision().statistics_adjust(0,1);

           FourierImageXmipp Filter_Out; Filter_Out().resize(ImgDivision());
	   if (!Periodogram_averaging) {
              // Perform an spectral estimation of the division by means of an ARMA model
	      cerr << "Performing ARMA spectral estimation\n";
      	      matrix2D<double> ARParameters,MAParameters;
      	      double dSigma=CausalARMA(ImgDivision(),ARMA_prm.N_AR,ARMA_prm.M_AR,
      		  ARMA_prm.N_MA,ARMA_prm.M_MA,ARParameters,MAParameters);

    	      // Get the AR filter coeficients in the fourier plane
    	      ARMAFilter(ImgDivision(),Filter_Out(),
		 ARParameters,MAParameters,dSigma);
	   } else {
              // Compute the periodogram
	      FourierTransform(ImgDivision(),Filter_Out());
	   }

      	   // Store the determined filter
	   if (!ARMA_averaging && !Periodogram_averaging) {
    	      ARMAfile.compose(ARMAfile.get_root().c_str(),N,"fft");
      	      ARMA_prm.fn_filter=ARMAfile; 
      	      Filter_Out.write(ARMA_prm.fn_filter);

	      // Estimate the CTF parameters of division 
	      adjust_CTF_prm.fn_ctf=ARMA_prm.fn_filter;
	      adjust_CTF_prm.fn_outroot=CTFfile.get_root()+ItoA(N,5);
	      adjust_CTF_prm.adjust(20)=adjust_CTF_prm.adjust(13)=
		 adjust_CTF_prm.adjust(0)=0;

  	      // The name of the parameters file that generates ROUT_Adjust_CTF
              if (compute_at_particle) {
		 FileName fn_current_particle=SF.get_current_file();
	         fn_current_particle = fn_current_particle.get_baseName();
		 adjust_CTF_prm.fn_out_CTF_parameters  =
		    fn_current_particle.add_prefix("ctf-")+".param";
      	      }	else 
		 adjust_CTF_prm.fn_out_CTF_parameters  =
      	             CTFfile.get_root()+ItoA(N,5)+".param";

	      // Compute the background and CTF
      	      cerr << "Estimating CTF parameters\n";
      	      adjust_CTF_prm.produce_side_info();
	      double fitting_error=ROUT_Adjust_CTF(adjust_CTF_prm);

   	      // Read parameters of the CTF from parameters file
	      XmippCTF pure_ctf; // An Xmipp Model of ctf
              pure_ctf.enable_CTFnoise=FALSE;
    	      pure_ctf.read(adjust_CTF_prm.fn_out_CTF_parameters);
	      pure_ctf.Tm=original_sampling_rate; // Change the sampling rate to the original one
	      pure_ctf.Produce_Side_Info();

      	      // Generate the pure CTF file determined from the parameters file	  
      	      cerr << "Generating CTF\n";
              FourierImageXmipp Img_pure_CTF;
              if (compute_at_particle)
 		 Img_pure_CTF().init_zeros(particle_vertical,particle_horizontal);
      	      else
      	         Img_pure_CTF().init_zeros(N_vertical,N_horizontal);
      	      pure_ctf.Generate_CTF(N_vertical, N_horizontal, Img_pure_CTF());

      	      // Save the image
              if (compute_at_particle) {
		 CTFfile=SF.get_current_file();
		 CTFfile = CTFfile.get_baseName();
		 CTFfile = CTFfile.add_prefix("ctf-")+".fft";
                 SF.next();		       
      	      } else
                 CTFfile.compose(CTFfile.get_root().c_str(),N,"fft");
      	      Img_pure_CTF.write(CTFfile);
	      
              // Compute CTF zeros as a checking point
              matrix1D<double> u(2), freq(2);
              VECTOR_R2(u,1,0); pure_ctf.zero(1,u,freq); X_zeros(N-1)=XX(freq);
              VECTOR_R2(u,0,1); pure_ctf.zero(1,u,freq); Y_zeros(N-1)=YY(freq);
              OutputZeros << CTFfile << " " << X_zeros(N-1) << " "
			  << Y_zeros(N-1) << " " << fitting_error << endl;
	   } else {
	      // Perform ARMA averaging
              FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(Filter_ARMA_avg())
	         MULTIDIM_ELEM(Filter_ARMA_avg(),i)=
		    abs(MULTIDIM_ELEM(Filter_ARMA_avg(),i));
	      Filter_Out()*=Filter_Out();
	      if (N==1) Filter_ARMA_avg()=Filter_Out();
	      else      Filter_ARMA_avg()+=Filter_Out();
	      if (N==1 && ARMA_averaging)
	         system("xmipp_show -img PPPARMAavg.fft -poll &");
	   }

      	   // Increment the division counter
           N++;
        } //while end
	M_in.close_micrograph();

	/*****************************************************************************
	  If ARMA_averaging, adjust the model to the average
	*****************************************************************************/
	FileName fn_ARMA_avg, fn_ARMA_avg_model;
	if (ARMA_averaging || Periodogram_averaging) {
	   FileName fn_root=image_fn.remove_directories();
	   fn_root=fn_root.without_extension();
	   Filter_ARMA_avg()/=div_Number;
	   if (Periodogram_averaging)
	      Filter_ARMA_avg()*=MULTIDIM_SIZE(Filter_ARMA_avg());
           FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(Filter_ARMA_avg())
	      MULTIDIM_ELEM(Filter_ARMA_avg(),i)=
		 sqrt(abs(MULTIDIM_ELEM(Filter_ARMA_avg(),i)));

	   // Write the Average file
	   fn_ARMA_avg=fn_root.add_prefix("ctf-")+".fft";
	   Filter_ARMA_avg.write(fn_ARMA_avg);

      	   // Write the average file amplitude
	   ImageXmipp save;
	   FFT_magnitude(Filter_ARMA_avg(),save());
	   CenterFFT(save(),true);
	   double min_positive=save().compute_max();
	   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(save()) {
	      double val=MULTIDIM_ELEM(save(),i);
	      if (val!=0) {
	         MULTIDIM_ELEM(save(),i)=10*log10(val*val);
	         if (val<min_positive) min_positive=val;
	      }
	   }
	   min_positive=10*log10(min_positive*min_positive)-1;
	   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(save())
	      if (MULTIDIM_ELEM(save(),i)==0) MULTIDIM_ELEM(save(),i)=min_positive;
	   save.write(fn_root.add_prefix("ctflog-")+".xmp");
	
	   // Estimate the CTF parameters
	   adjust_CTF_prm.fn_ctf=fn_ARMA_avg;
	   adjust_CTF_prm.fn_outroot=fn_root;
	   adjust_CTF_prm.adjust(20)=adjust_CTF_prm.adjust(13)=
	      adjust_CTF_prm.adjust(0)=0;
      	   adjust_CTF_prm.fn_out_CTF_parameters=
      	          adjust_CTF_prm.fn_outroot+".param";

	   // Compute the background and CTF
      	   cerr << "Estimating CTF parameters\n";
      	   adjust_CTF_prm.produce_side_info();
	   double fitting_error=ROUT_Adjust_CTF(adjust_CTF_prm);

   	   // Read parameters of the CTF from parameters file
	   XmippCTF pure_ctf; // An Xmipp Model of ctf
           pure_ctf.enable_CTFnoise=FALSE;
    	   pure_ctf.read(adjust_CTF_prm.fn_out_CTF_parameters);
	   pure_ctf.Tm=original_sampling_rate; // Change the sampling rate to the original one
	   pure_ctf.Produce_Side_Info();

      	   // Generate the pure CTF file determined from the parameters file	  
      	   cerr << "Generating CTF\n";
           FourierImageXmipp Img_pure_CTF;
           Img_pure_CTF().init_zeros(N_vertical,N_horizontal);
      	   pure_ctf.Generate_CTF(N_vertical, N_horizontal, Img_pure_CTF());

      	   // Save the image
	   fn_ARMA_avg_model=fn_root.add_prefix("ctfmodel-")+".fft";
      	   Img_pure_CTF.write(fn_ARMA_avg_model);

           // Generate the CTF at the particle scale
           Img_pure_CTF().init_zeros(particle_vertical,particle_horizontal);
      	   pure_ctf.Generate_CTF(particle_vertical,particle_horizontal,
	      Img_pure_CTF());
	   fn_ARMA_avg_model=fn_root.add_prefix("ctfmodelparticle-")+".fft";
      	   Img_pure_CTF.write(fn_ARMA_avg_model);

           // Compute CTF zeros as a checking point
           matrix1D<double> u(2), freq(2);
           VECTOR_R2(u,1,0); pure_ctf.zero(1,u,freq); X_zeros(0)=XX(freq);
           VECTOR_R2(u,0,1); pure_ctf.zero(1,u,freq); Y_zeros(0)=YY(freq);
           OutputZeros << CTFfile << " " << X_zeros(0) << " "
		       << Y_zeros(0) << " " << fitting_error << endl;
	}
        OutputZeros.close();

	/*****************************************************************************
	  Assign a CTF to every image in the sel file
	*****************************************************************************/
	if (!compute_at_particle) {
	   // Number of interger divisions in the micrograph
	   int Ydiv=Ydim/N_vertical;
	   int Xdiv=Xdim/N_horizontal; 

	   FourierImageXmipp  Filter_Out(N_vertical,N_horizontal);

	   // Process the Selfile
	   cerr << "Interpolating CTF for every particle ... " <<  endl;
	   SF.go_beginning();
	   while (!SF.eof()) {  
              string line;
              getline(PosFile,line);
              FileName file_with_ctf_for_current_particle;
              if (line[0]!='#') {
		 if (SF.Is_ACTIVE()) {
		    if (!ARMA_averaging && !Periodogram_averaging) {
                       // Read coordinates of the particle
                       int X,Y;
                       sscanf(line.c_str(),"%d %d",&X,&Y);

		       matrix1D<int> Regions(4);
		       DetermineRegionsToInterpolate(X,Y,N_horizontal,N_vertical,Xdiv,Ydiv,Regions);     
		       // Interpolate the CTF files for the necessary Regions to obtain
		       // a new CTF file specific for the particle's position  
		       CTFfile=CTFfile.get_root();
		       InterpolateRegionsCTF(X,Y,N_horizontal,N_vertical,Xdiv,Ydiv,
      	        	  Regions,CTFfile,Filter_Out());

      	               // Scale the filter to the size of the particle's images
      	               Filter_Out().self_scale_to_size(
			  particle_vertical,particle_horizontal);

      	               // Write the filter file and assign it to the particle's
		       // image in the output file.  
    		       FileName fn_current_particle=SF.get_current_file();
      	               file_with_ctf_for_current_particle=
      	                  fn_current_particle.add_prefix("ctf-")+".fft";
      	               Filter_Out.write(file_with_ctf_for_current_particle);
		    } else
      	               file_with_ctf_for_current_particle=fn_ARMA_avg_model;

      	            // Write in output file the name of the micrograph and
		    // the name of the associated CTF file 
                    OutputFile_ctf << file_with_ctf_for_current_particle
		                   << " 1\n";
		 }
              }
              SF.next();
	   }
        }
        OutputFile_ctf.close();
        PosFile.close();
   }
}
#undef DEBUG
