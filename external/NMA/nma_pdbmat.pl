#!/usr/bin/perl -w

open(A,"$ARGV[0]");

use Getopt::Long;

$debug = 1;

$pdb    = 'in.pdb';
$cutoff = 20.0;
$c_bond = 0;
$ignore_pair = undef;

# Read the input dat file
while($lineb = <A>){
     $lineb =~ s/,/ /g;
     ($a1,$b1) = split(" ", $lineb);
     $pdb=$a1;
     $cutoff=$b1;
     last;
}

if(defined $ignore_pair){
  ($a,$b) = split(' ',$ignore_pair);
  ($lb,$hb)=split('-',$a);
  foreach $n ($lb..$hb){
    $A{$n}=undef;
  }
  
  ($lb,$hb)=split('-',$b);
  foreach $n ($lb..$hb){
    $B{$n}=undef;
  }
}

print "group a ",join(' ',sort {$a<=>$b} keys %A),"\n" if $debug > 0;
print "group b ",join(' ',sort {$a<=>$b} keys %B),"\n" if $debug > 0;

open(OUT,"> matrice.sdijf") or die;
open(PDB,"$pdb") or die "failed to open $pdb, ";

print "cutoff      : $cutoff\n";
print "pdb         : $pdb\n";
print "ignore pair : $ignore_pair\n" if defined $ignore_pair;

$n=0;
while(<PDB>){
  next unless /^ATOM/;
  $n++;
  $x[$n] = substr($_,30,8);
  $y[$n] = substr($_,38,8);
  $z[$n] = substr($_,46,8);
}

$N = $n;

if($c_bond == 1){
  system("printf \"bond\nquit\n\" | rdparm 4ake.parm > blist 2>/dev/null");
  open(BOND,"blist") or die;
  while(<BOND>){
    $bond[$1][$2] = 1 if /\((\d+),(\d+)\)/;
  }
  
  print "# constant of bond     $c_bond\n";
  print "# constant of non bond 1.0\n";
}else{
  print "# sprint constant : 1.0\n";
}

$count=0;
$count_bond=0;
$count_nbond=0;
for $n (1..$N){
  for $m ($n+1..$N){

    if( defined $bond[$n][$m] ){
      $k = $c_bond;
      $count_bond++;
    }else{
      $k = 1;
      $count_nbond++;
    }
    $dx = $x[$n] - $x[$m];
    next if( abs($dx) > $cutoff);
    $dy = $y[$n] - $y[$m];
    next if( abs($dy) > $cutoff );
    $dz = $z[$n] - $z[$m];
    next if( abs($dz) > $cutoff );
    $rr = $dx**2 + $dy**2 + $dz**2;
    next if( $rr > $cutoff**2);

    if( defined $ignore_pair  &&
	(exists $A{$n} && exists $B{$m}) ||
	(exists $A{$m} && exists $B{$n}) ){
      print "cut $n <-> $m\n" if $debug >0;
      next;
    }

    $r = sqrt($rr);
    if ($r==0)
    {
       $drdx = 0.;
       $drdy = 0.;
       $drdz = 0.;
    }
    else
    {
       $drdx = $dx / $r;
       $drdy = $dy / $r;
       $drdz = $dz / $r;
    }
    
    $H[3*$n-2][3*$n-2] += $k*$drdx*$drdx;
    $H[3*$n-2][3*$n-1] += $k*$drdx*$drdy;
    $H[3*$n-2][3*$n  ] += $k*$drdx*$drdz;
    $H[3*$n-1][3*$n-1] += $k*$drdy*$drdy;
    $H[3*$n-1][3*$n  ] += $k*$drdy*$drdz;
    $H[3*$n  ][3*$n  ] += $k*$drdz*$drdz;
    
    $H[3*$m-2][3*$m-2] += $k*$drdx*$drdx;
    $H[3*$m-2][3*$m-1] += $k*$drdx*$drdy;
    $H[3*$m-2][3*$m  ] += $k*$drdx*$drdz;
    $H[3*$m-1][3*$m-1] += $k*$drdy*$drdy;
    $H[3*$m-1][3*$m  ] += $k*$drdy*$drdz;
    $H[3*$m  ][3*$m  ] += $k*$drdz*$drdz;
    
    $H[3*$n-2][3*$m-2] += -$k*$drdx*$drdx;
    $H[3*$n-2][3*$m-1] += -$k*$drdx*$drdy;
    $H[3*$n-2][3*$m  ] += -$k*$drdx*$drdz;
    $H[3*$n-1][3*$m-2] += -$k*$drdy*$drdx;
    $H[3*$n-1][3*$m-1] += -$k*$drdy*$drdy;
    $H[3*$n-1][3*$m  ] += -$k*$drdy*$drdz;
    $H[3*$n  ][3*$m-2] += -$k*$drdz*$drdx;
    $H[3*$n  ][3*$m-1] += -$k*$drdz*$drdy;
    $H[3*$n  ][3*$m  ] += -$k*$drdz*$drdz;
    
    $count++;
  }
}

print "# number of interactions $count, $count_bond(bond), $count_nbond(nonbond)\n";

$count=0;
for $n ( 0 .. ( 3*$N+2 ) ){
  for $m ($n..( 3*$N+2 )){
    if( defined $H[$n][$m] && $H[$n][$m] != 0 ){
      $count++;
    }
  }
}
print "# number of non-zero elements $count\n";

$count=0;
for $n ( 0 .. ( 3*$N+2 ) ){
  if( defined $H[$n][$n] && $H[$n][$n] != 0 ){
    $count+=$H[$n][$n];
    }
}
print "# matrix trace $count\n";

$count=0;
for $n ( 0 .. ( 3*$N+2 ) ){
  for $m ($n..( 3*$N+2 )){
    if( defined $H[$n][$m] && $H[$n][$m] != 0 ){
      printf(OUT "%6d %5d %15.12f\n",$n,$m,$H[$n][$m]);
    }
  }
}

sub distance {
  return sqrt( ($_[0]-$_[3])**2 + ($_[1]-$_[4])**2 + ($_[2]-$_[5])**2 );
}
