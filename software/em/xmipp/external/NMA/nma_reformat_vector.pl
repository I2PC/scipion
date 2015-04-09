#!/usr/bin/perl -w

@ARGV==2 or print "usage: reformat_vector.pl natom org_vector\n";

$natom = $ARGV[0];
$orgvec = $ARGV[1];

open(IN,$orgvec) or die;

$clm=0;
$row=0;
$vec=1;
open(OUT,">vec.$vec") or die;
while(){
    if( @F == 0 ){
	$_=<IN>;
	defined $_ or exit;
	chomp;
	@F = split;
    }

    if($row==$natom){
	$vec++;
	open(OUT,">vec.$vec") or die;
	$row=0;
    }
    
    $a = shift(@F);
    
    print OUT $a;
    if($clm==0 or $clm==1){
	print OUT " ";
	$clm++;
    }else{
	print OUT "\n";
	$clm=0;
	$row++;
    }
}
