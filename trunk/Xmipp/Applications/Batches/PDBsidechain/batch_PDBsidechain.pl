# This script reads a PDB file and creates a .descr file that
# can be used with create_phantom, project, etc. Each atom is replaced
# by a blob. Hidrogen atoms are not counted

# Get input parameters -----------------------------------------------------
$argc = @ARGV;
if ($argc < 2) {
   printf "Usage: PDB2descr \n";
   printf "          <in_file.pdb>          : Input PDB file\n";
   printf "          <out_file.pdb>         : Output PDB file\n";
   exit(1);
}

# Remove main chain --------------------------------------------------------
open(INFILE, "<$ARGV[0]");
open(OUTFILE,">$ARGV[1]");

while ( $tmp = <INFILE> ) {
  @separated  = split (/ +/,$tmp);
  if ($separated[0] =~ ATOM) {
     if (!($separated[2] eq "N") && !($separated[2] eq "CA") &&
         !($separated[2] eq "C") && !($separated[2] eq "O")) {
        printf OUTFILE "%s",$tmp;
     }
  } else {
     printf OUTFILE "%s",$tmp;
  }
}

close(INFILE);
close(OUTFILE);
