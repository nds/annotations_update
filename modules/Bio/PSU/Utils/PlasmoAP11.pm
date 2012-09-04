package Bio::PSU::Utils::PlasmoAP11;

# PlasmoAP version 1.1  - with 4 different possible outcomes!
# 11 September 2002
# Bernardo Foth

use strict;
use Data::Dumper;
use Carp;

# Bioperl for sequence FASTA file parsing
# use Bio::SeqIO;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name => value pairs (odd number supplied)";
  }

  my $self = {
	      seqobj   => undef,
	      seqId    => undef,
	      sequence => undef,
	      @_,
	     };

  bless ($self, $class);
  return $self;
}

####################### Definition of Variables #######################

# my $temp;
# my $sequence;
my $seqNumber = 0 ;

my $play;
my $playSub;
my $acid;
my $basic;
my $enrich;
my $enrichStart;
my $enrichEnd;
my $shift;

my $yes = "yes";
my $no = "no";

########################## Classic (old!!) values ###############################

my $lengthA = 12;
my $cutoffA = 2; ##### must be <= 2 acidic AA in first 12 AAs

my $lengthB = 17;
my $cutoffB = 1; ##### ratio acidic/basic must be <=1.0 in first 17 AAs (and there must be >= 1 basic AA)

my $lengthC = 17;
my $minimumKN = 6;
my $cutoffC = 40; ##### stretch with >=6 K+N per 17 AAs must occur at AA number <=40

my $cutoffD = 0.9; ##### ratio acidic/basic must be <=0.9 in KN-rich stretch from step C

####################################################################

my $valueA;
my $valueB;
my $valueC;
my $valueD;

my $verdictA;
my $verdictB;
my $verdictC;
my $verdictFinal;

my $finalScore;

####################### Preparations... #######################
# print "Your sequences are being processed...  (stop script by pressing apple+dot!)\n\n";

# open (INFILE, "infile.txt") or die ("Cannot open \"infile.txt\" for reading: $!");
# my $infile = $ARGV[0];
# my $sfasta = Bio::SeqIO->new (
#			      -file => $infile,
#			     );
# open (OUTFILE, ">outfile.txt") or die ("Cannot create or open \"outfile.txt\" for writing: $!"); 

####################### Analysis Start!! #######################
sub analysis {
        my $self = shift;
	my $seqobj;
	my $sequence;
	my $seqId;

	if (defined $self->{seqobj} && ($self->{seqobj}->isa ("Bio::PrimarySeq") || $self->{seqobj}->isa ("Bio::Seq"))) {
	  $seqobj   = $self->{seqobj};
	  $seqId    = $seqobj->display_id;
	  $sequence = $seqobj->seq;
	}
	elsif (defined $self->{seqId} && defined $self->{sequence}) {
	  $seqId    = $self->{seqId};
	  $sequence = $self->{sequence};
	}
	else {
	  print STDERR "ERROR - you need to give a Bioperl sequence or a sequence id and a protein sequence!!\n";
	  return "";
	}

	unless (length($sequence) == 0) {
		$sequence = uc($sequence);
		++$seqNumber;

		$finalScore = 0;

		############################################################
		########################  Parameters for CHECK1 ############
		############################################################

		$lengthB = 22;
		$cutoffB = 0.7;

		$lengthC = 40;
		$minimumKN = 9;
		$cutoffC = 40;

		$cutoffD = 0.9; 
		
		############### Analysis B CHECK1 #######################
		$verdictB = $no;		
		$play = substr($sequence, 0, $lengthB);
		$play =~ tr/DE/@/;
		$acid = $play =~ s/@/1/g;

		$play =~ tr/KRH/@/;
		$basic = $play =~ s/@/2/g;

		if ($basic > 0) {
			$valueB = $acid/$basic;	
			if ($valueB <= $cutoffB) {
				$verdictB = $yes;
			}
			else { $verdictB = $no; }
		}
		else { 
			$valueB = "infinite";
			$verdictB = $no;
		}

		if ($verdictB eq $yes) {
		
			############### Analysis C CHECK1 #######################   
			$verdictC = $no; 
			$play = substr($sequence, 0, ($cutoffC+$lengthC+10)); #add another 10 AAs just to be safe
			$play =~ tr/KN/@/;
			$shift = 0;
			$enrichStart = -1;
			$enrichEnd = -1;
	
			while ($shift < $cutoffC && $enrichStart == -1) {
				$playSub = substr($play, $shift, $lengthC);
				$enrich = $playSub =~ s/@/@/g;
				if ($enrich >= $minimumKN) {
					$enrichStart = $shift;
	
					$play = substr($sequence, $enrichStart, $lengthC);
		 
					$play =~ tr/DE/@/;
					$acid = $play =~ s/@/1/g;
			
					$play =~ tr/KRH/@/;
					$basic = $play =~ s/@/2/g;
		
					if ($basic > 0) {
						$valueD = $acid/$basic;	
						if ($valueD <= $cutoffD) {
							$verdictC = $yes;
						}
						else { $verdictC = $no;}
					}
					else { 
						$verdictC = $no;
					}
				}
				++$shift;
			}
		}

		############### Final Analysis CHECK1 #######################
		if (($verdictB eq $yes) && ($verdictC eq $yes) ) {
			$finalScore = $finalScore +2;
		}

		##############################################################
		###################  Parameters for CHECK2 ###################
	        ##############################################################

		$lengthA = 15;
		$cutoffA = 2;
		
		$lengthC = 40;
		$minimumKN = 9;
		$cutoffC = 40;

		$cutoffD = 0.6;

		############### Analysis A CHECK2 #######################
		$verdictA = $no;
		$play = substr($sequence, 0, $lengthA);
		$play =~ tr/DE/@/;
		$valueA = $play =~ s/@/@/g;
		if ($valueA < 1) {$valueA = 0;}
		if ($valueA <= $cutoffA) {
			$verdictA = $yes;
		}
		else { $verdictA = $no;	}

		if ($verdictA eq $yes) {

			############### Analysis C CHECK2 #######################   
			$verdictC = $no; 
			$play = substr($sequence, 0, ($cutoffC+$lengthC+10)); #add another 10 AAs just to be safe
			$play =~ tr/KN/@/;
			$shift = 0;
			$enrichStart = -1;
			$enrichEnd = -1;
	
			while ($shift < $cutoffC && $enrichStart == -1) {
				$playSub = substr($play, $shift, $lengthC);
				$enrich = $playSub =~ s/@/@/g;
				if ($enrich >= $minimumKN) {
					$enrichStart = $shift;
	
					$play = substr($sequence, $enrichStart, $lengthC);
		 
					$play =~ tr/DE/@/;
					$acid = $play =~ s/@/1/g;
			
					$play =~ tr/KRH/@/;
					$basic = $play =~ s/@/2/g;
		
					if ($basic > 0) {
						$valueD = $acid/$basic;	
						if ($valueD <= $cutoffD) {
							$verdictC = $yes;
						}
						else { $verdictC = $no;}
					}
					else { 
						$verdictC = $no;
					}
				}
				++$shift;
			}
		}

		############### Final Analysis CHECK2 #######################
		if (($verdictA eq $yes) && ($verdictC eq $yes) ) {
			$finalScore = $finalScore +2;
		}
		
		
		
		##############################################################		
		########## CHECK3: add one extra point if first charged#######
                ########## AA is basic (ie not acidic)!                #######
		##############################################################		

		$play = $sequence;
		$play =~ tr/DEKRH/11222/;
		$play =~ s/[A-z]//g;
		if (substr($play, 0, 1) == 2) { $finalScore++; }

		##############################################################		
		
	}



	######### Format output ##############		

	my $rate = $finalScore;
	if ($finalScore == 5) {$rate = "++"; }
	else { $rate =~ tr/01234/\-\-\-0\+/; }

	# print (OUTFILE "$seqId\t");
	# print (OUTFILE "$finalScore\t");
	# print (OUTFILE "$rate\n");

	my $result = "$seqId\t$finalScore\t$rate";

	return $result;

} # next sequence


####################### Cleaning up
# close INFILE or die ("Cannot close file : $!"); 
# close OUTFILE or die ("Cannot close file : $!"); 

# print "$seqNumber sequences processed...\n";
# print "\nDONE\n\n";

1;
