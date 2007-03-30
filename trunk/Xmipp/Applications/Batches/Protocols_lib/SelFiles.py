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
       lines=fh.readlines()
       self.sellines=[]
       for line in lines:
           self.sellines.append(line[:-1].split(" "))
       fh.close()

   # Fills selfile content
   def set(self,lines):
       self.sellines=lines

   # Appends selfile content
   def append(self,lines):
       self.sellines+=lines

   # Prints all lines to screen
   def printlines(self):
       for name,state in self.sellines:
           print name,state

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

if __name__ == '__main__':
   mysel=selfile()
   mysel.read('test.sel')
#   newsel=mysel.copy_sel('images')
#   newsel.printlines()
#   newsel=newsel.replace_string('xmp','app')
#   newsel.printlines()
#   newsel.write('new.sel')
   name,state=mysel.find_first_active_image()
   print name, state
