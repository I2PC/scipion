class ctfdat:
   """ handle ctfdat files
       Authors: Carlos Oscar Sánchez Sorzano
       May 2007
   """

   def __init__(self):
       self.lines=[]

   # Append a ctf
   def append(self,fnProjection,fnCTF):
      self.lines+=[fnProjection+' '+fnCTF+'\n']

   # Change directories
   def changeDirectory(self,projectionDir,CTFDir):
       newctfdat=ctfdat()
       for line in self.lines:
           # Read a line
           aux=line.split()
	   fnProjection=aux[0]
	   fnCTF=aux[1]
           splitFnProjection=fnProjection.split('/')
           splitFnCTF=fnCTF.split('/')

	   # Create corresponding line in the scaled ctfdat file
	   newctfdat.append(projectionDir+"/"+splitFnProjection[-1],
	      CTFDir+"/"+splitFnCTF[-1])
       return newctfdat

   # Produce a list with the ctfs
   def produceCTFList(self):
       CTFList=[]
       for line in self.lines:
           args=line.split()
           if not len(args)==2:
	       continue
	   fnCTF=args[1]
	   if not fnCTF in CTFList:
	      CTFList+=[fnCTF]
       return CTFList

   # Reads ctfdat from disk
   def read(self,filename):
       self.filename=filename
       fh=open(filename,'r')
       self.lines=fh.readlines()
       fh.close()

   # Write ctfdat to disk
   def write(self,filename):
       self.filename=filename
       fh=open(filename,'w')
       fh.writelines(self.lines)
       fh.close()
