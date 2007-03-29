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
       fh=open(selfilename,'r')
       self.sellines=fh.readlines()
       fh.close()

   # Fills selfile content
   def set(self,lines):
       self.sellines=lines

   # Appends selfile content
   def append(self,lines):
       self.sellines+=lines

   # Prints all lines to screen
   def printlines(self):
       for line in self.sellines:
           print line[:-1]

   # Writes selfile to disc
   def write(self,selfilename):
       fh=open(selfilename,'w')
       fh.writelines(self.sellines)
       fh.close()

   # Finds first active image and returns name and state
   # If no image is active the returned name='' and state=-1
   def find_first_active_image(self):
       dname=''
       dstate='-1'
       for line in self.sellines:
           name,state=line[:-1].split(" ")
           if state=='1':
               return name,state
       return dname,dstate

   # Adds 1 directory at the beginning of the path
   def add_1directory_begin(self,directory):
       newlines=[]
       for line in self.sellines:
           line=directory+'/'+line
           newlines.append(line)
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Adds 1 directory at the end of the path
   def add_1directory_end(self,directory):
       newlines=[]
       for line in self.sellines:
           parts=line.split('/')
           name=parts[-1]
           parts.pop()
           parts.append(directory)
           parts.append(name)
           line= '/'.join(parts)
           newlines.append(line)
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Removes N directories from beginning of path
   def remove_Ndirectories_begin(self,N):
       newlines=[]
       for line in self.sellines:
           parts=line.split('/')
           for i in range(N):
               parts.pop(0)
           line= '/'.join(parts)
           newlines.append(line)
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Removes N directories from end of path
   def remove_Ndirectories_end(self,N):
       newlines=[]
       for line in self.sellines:
           parts=line.split('/')
           name=parts[-1]
           for i in range(N+1):
               parts.pop()
           parts.append(name)
           line= '/'.join(parts)
           newlines.append(line)
       newsel=selfile()
       newsel.set(newlines)
       return newsel

   # Replaces "strin" with "strout"
   def replace_string(self,strin,strout):
       newlines=[]
       for line in self.sellines:
           line=line.replace(strin,strout)
           newlines.append(line)
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
       for line in self.sellines:
           name,state=line[:-1].split(" ")
           shutil.copy(name,directory)
           parts=name.split('/')
           line=directory+'/'+parts[-1]+' 1\n'
           newlines.append(line)
       newsel=selfile()
       newsel.set(newlines)
       return newsel

if __name__ == '__main__':
   mysel=selfile()
   mysel.read('test.sel')
   newsel=mysel.copy_sel('images')
   newsel.write('new.sel')
   name,state=mysel.find_first_active_image()
   print name, state
