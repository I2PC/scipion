class selfile:
   """ handle self files
       Authors: Roberto Marabini,
                Sjors Scheres
       March 2007
   """

   def __init__(self):
       self.sellines=[]

   # Reads selfile from disc
   def read(self,selfilename):
       self.selfilename=selfilename
       fh=open(self.selfilename,'r')
       lines=fh.readlines()
       self.sellines=[]
       for line in lines:
           args=line.split()
           if (len(args)>1):
               self.sellines.append([args[0],args[1]])
       fh.close()

   # Fills selfile content
   def set(self,lines):
       self.sellines=lines

   # Returns the size [xdim,ydim] of the images inside the selfile
   def imgSize(self):
      import spider_header
      header=spider_header.spiderheader((self.find_first_active_image())[0])
      return header.nx,header.ny

   # counts number of active entries in selfile
   def length(self):
       i=0
       state=-1
       for name,state in self.sellines: 
           if state=='1': i = i + 1
       return i
   
   # counts number of entries in selfile
   def lenght_even_no_actives(self):
       i=0
       state=-1
       for name,state in self.sellines: 
           if state=='1' or state=='-1': i = i + 1
       return i
      
   # Insert an image "name" with state "state" in the selfile
   def insert(self,name,state):
       self.sellines.append([name,state])

   # Appends selfile content
   def append(self,lines):
       self.sellines+=lines

   # Prints all lines to screen
   def printlines(self):
       for args in self.sellines:
           print args

   # Writes selfile to disc
   def write(self,selfilename):
       lines=[]
       for name,state in self.sellines:
           newline=name+' '+state+'\n'
           lines.append(newline)
       fh=open(selfilename,'w')
       fh.writelines(lines)
       fh.close()

   # Finds first active image and returns name and state
   # If no image is active the returned name='' and state=-1
   def find_first_active_image(self):
       dname=''
       dstate='-1'
       for name,state in self.sellines:
           if state=='1':
               return name,state
       return dname,dstate

   # Adds 1 directory at the beginning of the path
   def add_1directory_begin(self,directory):
       newlines=[]
       for name,state in self.sellines:
           name=directory+'/'+name
           newlines.append([name,state])
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Adds 1 directory at the end of the path
   def add_1directory_end(self,directory):
       newlines=[]
       for name,state in self.sellines:
           parts=name.split('/')
           name=parts[-1]
           parts.pop()
           parts.append(directory)
           parts.append(name)
           name= '/'.join(parts)
           newlines.append([name,state])
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Removes N directories from beginning of path
   def remove_Ndirectories_begin(self,N):
       newlines=[]
       for name,state in self.sellines:
           parts=name.split('/')
           for i in range(N):
               parts.pop(0)
           name= '/'.join(parts)
           newlines.append([name,state])
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Removes N directories from end of path
   def remove_Ndirectories_end(self,N):
       newlines=[]
       for name,state in self.sellines:
           parts=name.split('/')
           img=parts[-1]
           for i in range(N+1):
               parts.pop()
           parts.append(img)
           name= '/'.join(parts)
           newlines.append([name,state])
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Replaces "strin" with "strout" in name
   def replace_string(self,strin,strout):
       newlines=[]
       for name,state in self.sellines:
           name=name.replace(strin,strout)
           newlines.append([name,state])
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Copies all particles in the selfile to "directory"
   # "directory" is created if it does not exist
   # A new selfile is made in the current directory
   def copy_sel(self,directory):
       import os
       import shutil
       newlines=[]
       if not os.path.exists(directory):
           os.makedirs(directory)
       for name,state in self.sellines:
           shutil.copy(name,directory)
           parts=name.split('/')
           name=directory+'/'+parts[-1]
           newlines.append([name,state])
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Given two pattern selfiles with equal file order (pattern1,pattern2)
   # and a selfile that is a subset of pattern1,
   # return a selfile that is the corresponding subset from pattern2
   def make_corresponding_subset(self,pattern1,pattern2):
       newlines=[]
       for name,state in self.sellines:
           for i in range(len(pattern1.sellines)):
               name1,state1=pattern1.sellines[i]
               if (name1==name):
                   name2,state2=pattern2.sellines[i]
                   newlines.append([name2,state])
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Given two  selfiles s1 and s2 create s1 U s2
   # the order will be first entry of s1 then first entry of s2 etc..
   #same size for s1 and s2 is assumed
   def intercalate_union(self,selfile2):
       newlines=[]
       i=0
       for name1,state1 in self.sellines:
               name2,state2=selfile2.sellines[i]
               newlines.append([name1,state1])
               newlines.append([name2,state2])
               i = i+1
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Given two  selfiles s1 and s2 create s1 U s2
   # the order will be first entry of s1 then first entry of s2 etc..
   #same size for s1 and s2 is assumed
   def intercalate_union_3(self,selfile2,selfile3):
       newlines=[]
       i=0
       for name1,state1 in self.sellines:
               name2,state2=selfile2.sellines[i]
               name3,state3=selfile3.sellines[i]
               newlines.append([name1,state1])
               newlines.append([name2,state2])
               newlines.append([name3,state3])
               i = i+1
       newsel=selfile()
       newsel.set(newlines)
       return newsel

if __name__ == '__main__':
   mysel=selfile()
   mysel.read('test.sel')
   pat1=selfile()
   pat2=selfile()
   pat1.read('unt.sel')
   pat2.read('til.sel')
   newsel=mysel.make_corresponding_subset(pat1,pat2)
   newsel.printlines()
#   newsel=mysel.copy_sel('images')
#   newsel.printlines()
#   newsel=newsel.replace_string('xmp','app')
#   newsel.printlines()
#   newsel.write('new.sel')
   name,state=mysel.find_first_active_image()
   print name, state
