class spiderheader:
   """ handle spider headers
       Author:Roberto Marabini, March 2007
   """
   def __init__(self,filename):
      import sys
      import os
      import struct

      format_little = '<12f'
      format_big    = '>12f'
      fsize = struct.calcsize(format_little)
      raw_data_fh = open(filename,'rb')
      raw_data = raw_data_fh.read(fsize)
      data =struct.unpack(format_little,raw_data)
      self.mode = data[4]
      return_value = self.check_endianess()
      if(return_value==0):
         data =struct.unpack(format_big,raw_data)

      self.nz = data[0]
      self.ny = data[1]
      self.nx = data[11]

   # mode should be an small integer
   def check_endianess(self):
      if(self.mode==1 or self.mode==2 or self.mode==3): return 1
      return 0
   #set header position "position" with value "value" for
   #file   filename 

def set_header_position(filename,position=260,value=1):
         import struct 
         myfile = open(filename,'rb+')
         myfile.seek(position)
         myfile.write(struct.pack('f',value))
         myfile.close()

if __name__ == '__main__':
   myheader=spiderheader('../imagenes/imau08078.raw')
   print myheader.nx, myheader.ny, myheader.nz
