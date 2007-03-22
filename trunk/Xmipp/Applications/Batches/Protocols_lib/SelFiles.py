class selfile:
   """ handle self files
       Author:Roberto Marabini, March 2007
   """

   def __init__(self,selfilename):
      import csv
      self.reader=csv.reader(file(selfilename), delimiter=" ")
   
   def find_first_active_image(self):
      name=''
      stte='-1'
      for row in self.reader:
         name,state = row
         if state=='1':
            return name,state
      return name,state
            
if __name__ == '__main__':
   mysel=selfile('../untSelect.sel')
   name,state=mysel.find_first_active_image()
   print name, state
