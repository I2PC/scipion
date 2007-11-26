#/bin/sh
svn log | grep "|" | awk '{split($5,a,"-"); persona[$3]=persona[$3]+1; printf "%s %s-%s\n",$3,a[1],a[2]}' > inter.txt

echo '#\!/usr/bin/perl\
\
# Line frequency count\
%date = ();\
%author = ();\
while (<>) {\
    @words=split;\
    DOLLARauthor{DOLLARwords[0]}++;\
    DOLLARdate{DOLLARwords[1]}++;\
}\
\
printf "Activity summary per months\\n";\
foreach DOLLARdatei ( sort (keys %date) ) {\
    printf "%5d %s\\n", DOLLARdate{DOLLARdatei}, DOLLARdatei;\
}\
\
printf "\\n";\
printf "Activity summary per people\\n";\
foreach DOLLARauthori ( sort (keys %author) ) {\
    printf "%5d %s\\n", DOLLARauthor{DOLLARauthori}, DOLLARauthori;\
}\
' | sed 's/DOLLAR/$/' | sed 's/DOLLAR/$/' | sed 's/DOLLAR/$/' > inter.pl

chmod 755 inter.pl
./inter.pl < inter.txt
rm inter.pl inter.txt
