#!/bin/csh

set dirmodes = modes
set out = diag_arpack.eigenfacs
set headerfile = header.txt
set end = $1

echo " VECTOR    1       VALUE  2.2318E-08" >  $headerfile
echo " -----------------------------------" >> $headerfile
rm -f $out

set i = 0
beginloop:

if ($i < $end) then
@ i = $i + 1

cat $headerfile $dirmodes/vec.$i >> $out
goto beginloop
endif

rm -f $headerfile
