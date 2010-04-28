#!/usr/bin/perl -w

# move pdb coordinate along normal mode with amplitude dq 

open(A,"$ARGV[0]");
open(B,"$ARGV[1]");
$dq=$ARGV[2];

$i=0;

# Read PDB file
while(<A>){
       chomp($_);
       $atom=substr($_,0,6);

       if ($atom =~ m/REMARK/) {
         print "$_\n" ;
       }
       else { 
         if ($atom =~ m/ATOM  /) {
            $i++;
            $name[$i]= substr($_,12,4);
            $altloc[$i]=substr($_,16,1);
            $resname[$i]= substr($_,17,3);
            $chainid[$i]=substr($_,21,1);
            $resseq[$i]= substr($_,22,4);
            $icode[$i]=substr($_,26,1);
            $x[$i] = substr($_,30,8);
            $y[$i] = substr($_,38,8);
            $z[$i] = substr($_,46,8);
            $occup[$i] = substr($_,54,6);
            $bfact[$i] = substr($_,60,6);
            $serial[$i]=$i;
            #$mass[$i]=$bfact[$i];
         }
       }
}
close A;

$natom=$i; #total number of atoms

$i=0;

# Read the mode file
while($lineb = <B>){
     $lineb =~ s/,/ /g;
     ($a1,$b1,$c1) = split(" ", $lineb);
      $i++;
      $vx[$i]=$a1;
      $vy[$i]=$b1;
      $vz[$i]=$c1;
}

close B;

# Create output deformed PDB
for ($i=1; $i <= $natom; $i++){
    $x1[$i]=$x[$i]+$dq*$vx[$i];
    $y1[$i]=$y[$i]+$dq*$vy[$i];
    $z1[$i]=$z[$i]+$dq*$vz[$i];

    printf ("ATOM  %5d %4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",$serial[$i],$name[$i],$altloc[$i],$resname[$i],$chainid[$i],$resseq[$i],$icode[$i],$x1[$i],$y1[$i],$z1[$i],$occup[$i],$bfact[$i],$name[$i]);
}
printf ("END");
