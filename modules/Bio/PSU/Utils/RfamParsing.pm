package Bio::PSU::Utils::RfamParsing;

use strict;
use Data::Dumper;
use Carp;
use IO::File;

# PSU
#use Bio::PSU::Feature;
#use Bio::PSU::Range;

use Bio::Tools::HMMER::Domain;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name => value pairs (odd number supplied)";
  }

  my $self = {
      file      => undef,
      fh        => undef,
      debug     => 0,
      fileindex => 0,
      fragments => undef,
      @_,
  };

  bless ($self, $class);

  if (not defined $self->{fh} and defined $self->{file}) {
      my $fh = IO::File->new ($self->{file});
      $self->{fh} = $fh;
      print STDERR "$fh\n";
  }
  
   return $self;
}



sub next_domain{
    my $self  = shift;
    my $fh    = $self->{fh};
    my $debug = $self->{debug};
    my $fileindex = $self->{fileindex};
    my $fragments = $self->{fragments};
    
    if ($debug) {
	print STDERR "\nnext gene parsing...\n";
	print STDERR "starting at fileindex, $fileindex\n";
    }
    
    # initialisation stuff
    
    my $filepos = tell $fh;
    if ($debug) {
	print STDERR "filepos before repositioning: $filepos\n";
    }
    
    seek ($fh, $fileindex, 0) or die "can't seek on file, " . $self->{file} . "!\n";
    $filepos = tell $fh;
    if ($debug) {
	print STDERR "filepos after repositioning: $filepos\n";
    }
    
    if (eof ($fh)) {
	
	if ($debug) {
	    print STDERR "eof reached!\n";
	}
	
	return (undef, undef);
    }
    
    # parsing ...
    my $domain;

    if (my $line = $fh->getline){
	chomp $line;
	if ($debug){
	    print STDERR "$line\n";
	}
        $line =~ /(\S+)\s*Rfam\s+similarity\s+(\d+)\s*(\d+)\s*(\S+)\s*(\+|-)\s+\.\s+\S+;model_end=(\d+);model_start=(\d+);rfam-acc=(\S+);rfam-id=(\S+)/;
        my $id = $1;
        my $sstart = $2;
        my $send = $3;
	my $score = $4;
        my $acc = $8;
        my $dstart = $7;
        my $dend = $6;
        my $name = $9;

	my $strand = "+1";
	$strand = "-1" if ($5 eq '-');
	#my ($id, $sstart, $send, $acc, $dstart, $dend, $score, $name) = split ("\t", $line);


        if ($fragments){ # rfam-scan.pl now deals with this so don't use -f option
            my $frags;
            ($id, $frags) = split (":", $id);
            my ($fstart, $fend) = split ("-", $frags);
            $sstart += ($fstart - 1);
            $send += ($fstart - 1);
        }

	$domain = Bio::Tools::HMMER::Domain->new;

	#my $strand = 1;
	#if ($sstart > $send){
	#    $strand = -1;
	#    my @tmp = ($sstart, $send);
	#    $sstart = $tmp[1];
	#    $send = $tmp[0];
	#}	

	$domain->hmmacc($acc);
	$domain->hmmname($name);
	$domain->strand($strand);
	$domain->start($sstart);
	$domain->end($send);
	$domain->bits($score);
	$domain->hstart($dstart);
	$domain->hend($dend);
	$domain->seq_id($id);
    }
    
    $fileindex++;
    $self->{fileindex} = tell ($fh);
    return $domain;
}


1;
