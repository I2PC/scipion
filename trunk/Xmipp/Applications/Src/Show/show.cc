/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
 *              Carlos Oscar Sanchez Sorzano (coss@cnb.uam.es)
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

#include <XmippData/xmippArgs.hh>
#include <XmippGraphics/show2D.hh>
#include <XmippGraphics/showSel.hh>
#include <XmippGraphics/showVol.hh>
#include <XmippGraphics/showSpectra.hh>
#include <XmippGraphics/showSOM.hh>
#include <XmippGraphics/showSpectraSOM.hh>
#include <qapplication.h>

void Usage();

int main( int argc, char **argv ) {
    int numCols, numRows, mode, ifirst;
    FileName fn_dat;
    bool poll;
    try {
       if (check_param(argc,argv,"-img"))          {mode=0; ifirst=position_param(argc,argv,"-img");}
       else if (check_param(argc,argv,"-sel"))     {mode=1; ifirst=position_param(argc,argv,"-sel");}
       else if (check_param(argc,argv,"-vol"))     {mode=2; ifirst=position_param(argc,argv,"-vol");}
       else if (check_param(argc,argv,"-spect"))   {mode=3; ifirst=position_param(argc,argv,"-spect");}
       else if (check_param(argc,argv,"-som"))     {mode=4; ifirst=position_param(argc,argv,"-som");}
       else if (check_param(argc,argv,"-spectsom")){
          mode=5; ifirst=position_param(argc,argv,"-spectsom");
	  fn_dat=get_param(argc,argv,"-din");
       } else
          REPORT_ERROR(1,"No mode (img/sel/vol) supplied");
       numCols = AtoI(get_param(argc, argv, "-w", "10"));
       numRows = AtoI(get_param(argc, argv, "-h", "10"));
       poll=check_param(argc,argv,"-poll");
   } catch (Xmipp_error) {Usage(); exit(1);}

   try {
   QApplication::setFont( QFont("Helvetica", 12) );
   QApplication a(argc,argv);

   int shown=0;
   for ( int i=ifirst+1; i<argc; i++ ) {
       if (!exists(argv[i])) {
          if (argv[i][0]=='-') break; // There is nothing else to show
          FileName fn;
          switch (mode) {
	     case 0: continue;
	     case 1: continue;
	     case 2:
	        fn=argv[i];
		if (fn[fn.length()-1]=='x' || fn[fn.length()-1]=='y') {
		   fn=fn.substr(0,fn.length()-1);
		   if (!exists(fn.c_str())) continue;
		} else continue;
                break;
	     case 3: continue;
	     case 4: break;
	     case 5: break;
	  }
       } 
       if (mode==0) {
          ImageViewer *showimg = new ImageViewer(argv[i], poll);
          showimg->loadImage( argv[i] );
          showimg->show();
          shown++;
       } else if (mode==1) {
          ShowSel *showsel = new ShowSel;
          showsel->showonlyactive = !check_param(argc,argv,"-showall");
	  showsel->initWithFile(numRows, numCols, argv[i]);
	  showsel->show();
          shown++;
       } else if (mode==2) {
          ShowVol *showvol = new ShowVol;
	  if (poll) showvol->setPoll();
	  showvol->initWithFile(numRows, numCols, argv[i]);
	  showvol->show();
          shown++;
       } else if (mode==3) {
          ShowSpectra *showspectra=new ShowSpectra;
	  showspectra->initWithFile(numRows, numCols, argv[i]);
	  showspectra->show();
          shown++;
       } else if (mode==4) {
          ShowSOM *showsom=new ShowSOM;
	  showsom->initWithFile(argv[i]);
	  showsom->show();
          shown++;
       } else if (mode==5) {
          ShowSpectraSOM *showspectrasom=new ShowSpectraSOM;
	  showspectrasom->initWithFile(argv[i],fn_dat);
	  showspectrasom->show();
          shown++;
       }
   }
   
   if (!shown) return 0;

   QObject::connect(qApp, SIGNAL(lastWindowClosed()), qApp, SLOT(quit()));
   return a.exec();
   } catch (Xmipp_error XE) {cerr << XE;}
}

void Usage() {
    cout << "Usage: showsel [options]\n"
         << "    -img <images> |       : Input images\n"
         << "    -sel <selfiles> |     : Input selfiles\n"
         << "    -vol <XmippVolumes>   : Add x or y to the filename\n"
	 << "                            to see slices in that direction\n"
         << "    -spect <datafile>     : Spectra .dat file\n"
	 << "    -som <SOM rootname>   : SOM images\n"
	 << "    -spectsom <SOM root>  : SOM spectra\n"
	 << "       -din <Original.dat>: Original data\n"
         << "   [-w]                   : width (default: 10)\n"
         << "   [-h]                   : height (default: 10)\n"
         << "   [-showall]          : only for sel mode, show all images even\n"
         << "                         if sel = -1\n"
         << "   [-poll]                : check file change, NOT for sel files\n"
    ;
}
