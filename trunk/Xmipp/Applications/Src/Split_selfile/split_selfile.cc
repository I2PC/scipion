/***************************************************************************
 *
 * Authors:     Sjors Scheres (scheres@cnb.uam.es)
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

#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippArgs.hh>

void Usage();
 
int main( int argc, char **argv ) {
  FileName fn_in,fn_out,fn_root;
  SelFile  SFin,SFout,SFtmp,SFtmp2;
  SelLine  line;
  bool     dont_randomize;
  int N;

  try {
    fn_in=get_param(argc,argv,"-i");
    N=AtoI(get_param(argc,argv,"-n","2"));
    fn_root=get_param(argc,argv,"-o","");
    dont_randomize=check_param(argc,argv,"-dont_randomize");
    if (fn_root=="") fn_root=fn_in.without_extension();
    SFin.read(fn_in);
  }  catch (Xmipp_error) {Usage(); exit(1);}

  try {

    if (!dont_randomize) SFtmp=SFin.randomize();
    else                 SFtmp=SFin;
    int Nsub=ROUND((double)SFtmp.ImgNo()/N);
    for (int i=0; i<N; i++) {
      SFout.clear();
      SFout.reserve(Nsub);
      SFtmp.go_beginning();
      SFtmp.jump_lines(Nsub*i);
      if (i==N-1) Nsub=SFtmp.ImgNo()-i*Nsub;
      for (int nn=0; nn<Nsub; nn++) {
	SFout.insert(SFtmp.current());
	SFtmp.NextImg();
      }
      SFtmp2=SFout.sort_by_filenames();
      SFout=SFtmp2;
      string num="_"+ItoA(i+1);
      fn_out=fn_root+num;
      fn_out+=".sel";
      SFout.write(fn_out);
    }


  } catch (Xmipp_error) {cerr <<"ERROR, exiting..."<<endl; exit(1);}

}
                                                                                                          
void Usage() {
    cout << "Usage: split_selfile [options]\n"
         << "    -i <selfile>            : Input selfile\n"
         << "  [ -n <int=2> ]            : Number of output selfiles\n"
         << "  [ -o <rootname=selfile> ] : Rootname for output selfiles\n"
         << "                              output will be: rootname_<n>.sel\n"
         << "  [ -dont_randomize ]       : Do not generate random groups\n"
    ;
}
