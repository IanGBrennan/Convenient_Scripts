#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# SMRT-raxMRP.pl v. 1.0 by Edward L. Braun, Department of Biology,
# University of Florida, Gainesville FL 32611 USA
#
# This script implements the SMRT-ML method described in:
#   DeGiorgio M., Degnan JH. 2010 Fast and consistent estimation 
#   of species trees using supermatrix rooted triples. Mol Biol
#   Evol. 27, 552-569. (doi:10.1093/molbev/msp250)
#
# It was modified from the script used in:
#   Sun K, Meiklejohn KA, Faircloth BC, Glenn TC, Braun EL, Kimball
#   RT. 2014. The evolution of peafowl and other taxa with ocelli
#   (eyespots): A phylogenomic approach. Proc. R. Soc. B 281, 20140823.
#   (doi:10.1098/rspb.2014.0823)
#
# if you use this script please cite the SMRT-ML paper (above) and
# the publication describing the first use of this script:
#   Meiklejohn KA, Faircloth BC, Glenn TC, Kimball RT, Braun EL.
#   submitted. Analysis of a rapid evolutionary radiation using 
#   ultraconserved elements (UCEs): Evidence for a bias in some 
#   multi-species coalescent methods.
# -- upon final publication full bibliographic material will be 
#    available from Braun's website and the "Early Bird" website:
#    - http://people.biology.ufl.edu/ebraun 
#    - http://biology.ufl.edu/earlybird
############################################################################

############################################################################
# Store locations of programs called by this program
#
#   Note that this program has two dependencies:
#      - RAxML -- store path to a raxml executable in $raxexec
#        it may be necessary to change the RAxML system call 
#        on line 202 depending upon the RAxML version used
#      - clann -- a suprtree program by Creevey and McInerney
#        that implements several methods
#   It should be possible to substitute other programs both for
#   the phylogenetic estimation and the supertree building. Users
#   interested in doing should change the names of the executables
#   and alter the system calls (lines 202 and 217) appropriately.
#    
############################################################################

my($progname) = $0;
my($raxexec) = "raxmlHPC"; # location of the raxml executable
my($clannexec) = "clann"; # location of the clann executable

############################################################################
# Initialize variables
############################################################################

my($iter);
my($jter);
my($kter);
my($lter);
my($mter);

if ( @ARGV != 3) {
	print "Usage:\n  \$ $progname <authority> <infile> <outfile>\n";
	print "  authority = file listing taxa (outgroup listed first)\n";
	print "  infile    = (relaxed) phylip format infile\n";
	print "  outfile   = base output file name\n";
	print "      raxml rooted triples written to <outfile>.tre\n";
	print "      MRP matrix written to <outfile>.MRP.nex\n";
	print "exiting...\n";
	exit;
}

my($authority)=$ARGV[0];
my($infile)=$ARGV[1];
my($outfile)=$ARGV[2];

############################################################################
# Read the input files
############################################################################
if (-e "$outfile.tre") {
	unlink("$outfile.tre"); # remove outfile.tre if it exists
}
if (-e "$outfile.MRP.nex") {
	unlink("$outfile.MRP.nex"); # remove outfile.MRP.nex if it exists
}

print "Reading the authority file... ";
open (my $AUTHF, $authority) or die "Could not open file $authority for input.\n";
my @taxlist = <$AUTHF>; # Read the authority file
close($AUTHF) or die "Could not close file $authority.\n";
my($taxnum) = $#taxlist + 1;
print "$taxnum taxa read:\n";

my($ogname) = $taxlist[0];
chomp($ogname);
print "  $ogname (outgroup)\n";
for ($iter=1; $iter<$taxnum; $iter++) {
	print "  $taxlist[$iter]";
}
print "\n";

print "Reading data file...";
open (my $INF, $infile) or die "Could not open file $infile for input.\n";
my @seqdata = <$INF>; # Read the infile file
close($INF) or die "Could not close file $infile.\n";
my($seqnum) = $#seqdata + 1;

my($ntax);
my($nchar);
($ntax,$nchar) = split(/\s+/, $seqdata[0]);
print "$ntax taxa and $nchar sites\n";

my @seqline;
my($ognum);
# identify the outgroup
for ($iter=0; $iter<$ntax; $iter++) {
	@seqline = split(/\s+/, $seqdata[$iter+1]);
#	print "$seqline[0]\n";
	if ( $seqline[0] eq "$ogname" ) {
		$ognum = $iter+1;
	}
}
# print "$ognum\n";

############################################################################
# Conduct analyses of all triples
############################################################################
my($itax)=$taxnum-2;
my($focaltaxon1);
my($focaltaxon2);
my($focaltaxon3);
# print "$itax\n";
my($ingroupnum) = $taxnum-1;
my($tripletnumber) = $ingroupnum * ($ingroupnum-1) * ($ingroupnum-2) / 6;
print "\nAnalyze $ingroupnum ingroup species ($tripletnumber triplets)\n";
print "\nAnalyzing rooted triplets:\n";

my($triplecount) = 0;
for ($iter=1; $iter<$itax; $iter++) {
	$focaltaxon1 = $taxlist[$iter]; # first focal taxon
	chomp($focaltaxon1);
	
	for ($jter=$iter; $jter<$itax; $jter++) {
		$focaltaxon2 = $taxlist[$jter+1]; # second focal taxon
		chomp($focaltaxon2);
		
		for ($kter=$jter; $kter<$itax; $kter++) {
			$focaltaxon3 = $taxlist[$kter+2]; # third focal taxon
			chomp($focaltaxon3);
			$triplecount++;
			print " triplet $triplecount:  $focaltaxon1 -- $focaltaxon2 -- $focaltaxon3";
			
			open (my $SF, ">scratch.infile") or die "Could not open file scratch.infile for output.\n";

			print $SF " 4 $nchar\n";
			# first print the outgroup
			print $SF "$seqdata[$ognum]";
			# then find the first ingroup taxon
			for ($lter=1; $lter<=$ntax; $lter++) {
				@seqline = split(/\s+/, $seqdata[$lter]);
				if ( $seqline[0] eq "$focaltaxon1" ) {
					print $SF "$seqdata[$lter]";
					$lter=$ntax;
				}
			}
			# then find the second ingroup taxon
			for ($lter=1; $lter<=$ntax; $lter++) {
				@seqline = split(/\s+/, $seqdata[$lter]);
				if ( $seqline[0] eq "$focaltaxon2" ) {
					print $SF "$seqdata[$lter]";
					$lter=$ntax;
				}
			}
			# then find the third ingroup taxon
			for ($lter=1; $lter<=$ntax; $lter++) {
				@seqline = split(/\s+/, $seqdata[$lter]);
				if ( $seqline[0] eq "$focaltaxon3" ) {
					print $SF "$seqdata[$lter]";
					$lter=$ntax;
				}
			}
			
			close($SF) or die "Could not close file scratch.infile.\n";
			
			print "  ...conducting RAxML analysis...";
			# RAxML call - check that it is appropriate for your raxml installation
			system("$raxexec -s scratch.infile -n scratchfile -p 1 --no-seq-check -m GTRGAMMA > screen");
			print "done\n";
			system("cat RAxML_bestTree.scratchfile >> $outfile.tre");
			system("rm RAxML*.scratchfile");
			
		} # $kter loop
	} # $jter loop
} # $iter loop

print "\nRun complete. Triplet treefiles written to $outfile.tre\n";

# execute the clann speciestree program
print "\nCoding MRP matrix using clann supertree program...";

# write a control file for clann to analyze the raxml data
open (my $RCLNF, ">clann.ctl") or die "Could not open file clann.ctl for output.\n";
print $RCLNF "execute '$outfile.tre'\n";
print $RCLNF "set criterion mrp\n";
print $RCLNF "hs\n";
print $RCLNF "quit\n";
close($RCLNF) or die "Could not close file clann.ctl\n";
		
# call clann and set up the paup file
system("$clannexec < clann.ctl > screen 2>&1");
system("sed -n '/#NEXUS/,/begin paup;/p' coding.nex > SMRT-ML-temp.0");
system("sed '\$d' SMRT-ML-temp.0 > $outfile.MRP.nex");

# append a paup block for the MRP analysis
open (my $RPBF, ">>$outfile.MRP.nex") or die "Could not open file $outfile.MRP.nex for output.\n";
print $RPBF "begin paup;\n";
print $RPBF "\tSet increase=auto;\n";
print $RPBF "\toutgroup $ogname;\n";
print $RPBF "\tPSet collapse=minBrlen;\n";
print $RPBF "\tHSearch addSeq=random nreps=10;\n";
print $RPBF "\tConTree / treeFile=$outfile.MRP.tre replace;\n";
print $RPBF "end;\n\n";
close($RPBF) or die "Could not close file $outfile.MRP.nex\n";

# clean up the directory
unlink("clann.ctl");
unlink("coding.nex");
unlink("scratch.infile");
unlink("screen");

print "done\n";
print "Execute $outfile.MRP.nex in PAUP*\n";