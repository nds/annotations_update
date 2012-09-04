package Bio::PSU::Utils::ORF_Finder;

use strict;
use Data::Dumper;
use Carp;

use Bio::PSU::Seq;

my $verbose = 0;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name => value pairs (odd number supplied)";
  }

  my $self = {
	      verbose  => 0,
	      @_,
	     };
  bless ($self, $class);
  $verbose = $self->{verbose};

  return $self;
}

sub getMaxORF {
  my $self = shift;
  my ($mRNA_seqobj, $is_partial) = @_;

  # in aminoacids
  my $min_ORF_length = 40;

  # Get the set of potential ORFs as translated sequences (string)

  my @ORFs = getORFs ($mRNA_seqobj, $min_ORF_length);

  my $ORF_sequence;
  if (not $is_partial) {
    
    if ($verbose) {
      print STDERR "complete mRNA mode\n";
    }

    # Working on full length mRNAs, so we expect to have the first Methionine
    
    my @ORFs_firstMet = ();
    foreach my $ORF_sequence (@ORFs) {
      $ORF_sequence = trimToFirstMet ($ORF_sequence, $min_ORF_length);
      push (@ORFs_firstMet, $ORF_sequence);
    }
    
    if (@ORFs_firstMet > 0) {
      $ORF_sequence = max (@ORFs_firstMet);
    }
  }
  else {

    # partial mRNA sequences, don't predict the first Methionine

    if ($verbose) {
      print STDERR "partial mRNA mode\n";
    }

    $ORF_sequence = max (@ORFs);

  }

  if (defined $ORF_sequence) {

    if ($verbose) {
      print STDERR "**********************\ngot this translation, \n$ORF_sequence\n**********************\n\n";
    }
    
    
    my $pept_seqobj = Bio::PSU::Seq->new (
					  -id  => $mRNA_seqobj->id,
					  -str => $ORF_sequence,
					 );
    
    return $pept_seqobj;
  }
  else {
    print STDERR "didn't find any ORF!\n";
    return undef;
  }
}

#####################
#
# Private methods
#
#####################

sub getORFs {
  my ($mRNA_seqobj, $min_ORF_length) = @_;
  my @ORFs = ();

  my $mRNA_reversed_seqobj = $mRNA_seqobj->revcom;
 
  # processing the 6-frame translation

  foreach my $mRNA_seqobj_temp ($mRNA_seqobj, $mRNA_reversed_seqobj) { 
    my $i = 1;
    while ($i < 4) {
      my $pept_seqobj = $mRNA_seqobj_temp->translate (1,$i,0);

      if ($verbose) {
	print STDERR "translated sequence:\n" . $pept_seqobj->str . "\n\n";
      }

      my @sub_ORFs = getORFsfromSequence ($pept_seqobj->str, $min_ORF_length);

      push (@ORFs, @sub_ORFs);

      $i++;
    }
  }

  return @ORFs;
}


sub getORFsfromSequence {
  my ($pept_sequence, $min_ORF_length) = @_;
  my @sub_ORFs = ();
  my $pept_sequence_original = $pept_sequence;
  my $sequence_length = length ($pept_sequence_original);

  # ie all the sub strings between two stop codon (*), above the cutoff length

  my $index_start = 0;

  if (not ($pept_sequence =~ /\*/)) {
    # no stop codon, return the whole sequence
    if (length ($pept_sequence) >= $min_ORF_length) {

      if ($verbose) {
	print STDERR "no stop codon in the translated sequence, returning the whole sequence\n";
      }

      return ($pept_sequence);
    }
  }

  # foreach stop codons
  while ($pept_sequence =~ m/\*/g) {
    my $index_stop = pos ($pept_sequence) - 1;
    my $length = $index_stop - $index_start;
    my $subsequence = substr ($pept_sequence_original, $index_start, $length);
    
    if ($verbose) {
      print STDERR "\tsubsequence:\n\t$subsequence\n";
    }
    
    if (length ($subsequence) >= $min_ORF_length) {
      push (@sub_ORFs, $subsequence);
    }
    
    $index_start = $index_stop + 1;

  }

  # Get the last subsequence, ie between the last stop codon and the end of the sequence

  my $subsequence = substr ($pept_sequence_original, $index_start, length ($pept_sequence_original));
  
  if ($verbose) {
    print STDERR "\tlast subsequence:\n\t$subsequence\n\n";
  }

  if (length ($subsequence) >= $min_ORF_length) {
    push (@sub_ORFs, $subsequence);
  }

  return @sub_ORFs;
}


sub max {
  my (@ORFs) = @_;
  my $max_ORF = "";

  foreach my $ORF (@ORFs) {
    if (length ($ORF) > length ($max_ORF)) {
      $max_ORF = $ORF;
    }
  }

  return $max_ORF;
}


sub trimToFirstMet {
  my ($ORF_sequence, $min_ORF_length) = @_;

  if ($verbose) {
    print STDERR "sequence length before trimming, " . length ($ORF_sequence) . "\n";
  }

  my $trimmed_ORF_sequence = $ORF_sequence;
  $trimmed_ORF_sequence =~ m/M/g;
  my $met_index = pos $trimmed_ORF_sequence;

  if (defined $met_index) {
    if ($verbose) {
      print STDERR "Methionine index, $met_index\n";
    }
    
    substr ($trimmed_ORF_sequence, 0, $met_index-1, "");
    
    if ($verbose) {
      print STDERR "sequence length after trimming, " . length ($trimmed_ORF_sequence) . "\n";
    }
    
    if (length ($trimmed_ORF_sequence) > $min_ORF_length) {
      return $trimmed_ORF_sequence;
    }
    else {
      if ($verbose) {
	print STDERR "after trimming the predicted translated sequence to the first Methionine, the sequence gets below the cutoff size!\n";
      }
      return $ORF_sequence;
    }
  }
  else {
    if ($verbose) {
      print STDERR "predicted translation doesn't contain any Methionine!\n";
    }
    return $ORF_sequence;
  }
}


1;
