#!/usr/bin/env python
#---------------------------------------------------------------------------
# getComponentFromVector
#---------------------------------------------------------------------------
def getComponentFromVector(_vector,_iteration):
   listValues=getListFromVector(_vector)
   if _iteration<0: _iteration=0
   if _iteration<len(listValues): return listValues[_iteration]
   else:                          return listValues[len(listValues)-1]

#---------------------------------------------------------------------------
# getListFromVector
#---------------------------------------------------------------------------
def getListFromVector(_vector):
   intervalos=string.split(_vector)
   if len(intervalos)==0:
      raise RuntimeError,"Empty vector"
   listValues=[]
   for i in range(len(intervalos)):
      intervalo=intervalos[i]
      listaIntervalo=string.split(intervalo,'x')
      if len(listaIntervalo)==1:
	 listValues+=listaIntervalo
      elif len(listaIntervalo)==2:
         listValues+=[ listaIntervalo[1] ] * string.atoi(listaIntervalo[0])
      else:
         raise RuntimeError,"Unknown syntax: "+intervalos
   return listValues
