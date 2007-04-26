class docfile:
   """ handle doc files files with format
 ; Headerinfo columns: rot (1) , tilt (2), psi (3), Xoff (4), Yoff (5), Weight (6), Mirror (7)
 ; proj_match_class00001.xmp
    1 7     0.00000   40.00000    0.00000    0.00000    0.00000   76.00000    0.00000
 ; proj_match_class00002.xmp
    2 7     7.65957   40.00000    0.00000    0.00000    0.00000   14.00000    0.00000
 first line title
 even lines filenames
 odd lines data   
       Authors: Roberto Marabini,
       March 2007
   """
   def __init__(self,docfilename):
       self.docfilename =docfilename
       f = open(self.docfilename, 'r')
       self.lineLst = []
       #skip first line
       self.firstline=' ; Headerinfo'
       self.firstline=f.readline()
       
       for line in f:
           line = line.strip()
           if not line.startswith(';'):
               mylist = (line.split())
               mylist.append(filename)
               self.lineLst.append(mylist)
           else:
               #store file name
               filename=(line.split())[1]
          
       f.close()
       
 
   def print_docfile(self):
       for i in range(len(self.lineLst)):
           print self.lineLst[i]
 
   def maximum_of_column(self,column):
       import string
       maximum=-9999999.
       for i in range(len(self.lineLst)):
           tmp = float(self.lineLst[i][column]) 
           if  tmp > maximum:
              maximum = tmp
       return maximum
       
   def minimum_of_column(self,column):
       import string
       maximum=+9999999.
       for i in range(len(self.lineLst)):
           tmp = float(self.lineLst[i][column]) 
           if  tmp < maximum:
              maximum = tmp
       return maximum
       
   def compare_7(self,a, b):
        return cmp(float(a[7]), float(b[7]))
        
   def sort(self,col):
        if (col==7):
           self.lineLst.sort(self.compare_7)
        else:
           print "sort not implemented"
           exit(1)   

   def write(self, docfilename):
       fh = open(docfilename, 'w')
       newline=self.firstline
       size = len(self.lineLst[0])
       
       for i in range(len(self.lineLst)):
          newline = newline + ' ; ' + self.lineLst[i][size-1] + '\n'
          newline = newline + ' ' + str(i) + ' ' + str(size) 
          for j in range(0,size-2):
              newline=newline + ' ' +self.lineLst[i][j]
          newline=newline + '\n'
       fh.writelines(newline)
       fh.close()
       
   def write_several(self,docfilename,bin,col,min,max):
   
       for h in range(0,bin):
            
           fh = open(docfilename+str(h)+'.doc', 'w')
           newline=self.firstline
           size = len(self.lineLst[0])-1
           _bin=(max-min)/bin
           i=0
           for ii in range(len(self.lineLst)):
              if float(self.lineLst[ii][col]) < min+h*_bin or\
                 float(self.lineLst[ii][col]) > min+(h+1)*_bin:
                     continue
              i=i+1
              newline = newline + ' ; ' + self.lineLst[ii][size] + '\n'
              newline = newline + ' ' + str(i) + ' ' + str(size-2) 
              for j in range(2,size):
                  newline=newline + ' ' +self.lineLst[ii][j]
              newline=newline + '\n'
           fh.writelines(newline)
           fh.close()
           import os
           if i==0: os.unlink(docfilename+str(h)+'.doc')
    
   def check_angle_range(self):
           for ii in range(len(self.lineLst)):
              rot=float(self.lineLst[ii][2])
              tilt=float(self.lineLst[ii][3])
              if tilt > 90.: 
                 tilt = -(tilt-180)
                 rot  += 180. 
              if tilt < -90.: 
                 tilt = -(tilt+180)
                 rot  -= 180. 
              self.lineLst[ii][2] = str(rot)  
              self.lineLst[ii][3] = str(tilt)
                 
      
