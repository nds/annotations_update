package Bio::PSU::Utils::Dgpi;

use strict;
use Carp;

use Data::Dumper;

use IO::File;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name=>value pairs (odd number supplied)";
  }

  my $self     = {
		  fh => undef,
                  file => undef,
		  debug => 0,
		  @_,
		 };
  bless ($self, $class);

  if (not defined $self->{fh} and defined $self->{file}) {
    my $fh = IO::File->new ($self->{file});
    $self->{fh} = $fh;
  }

  return $self;
}

sub next_prediction {
  my $self = shift;
  my $fh = $self->{fh};
  my $debug =  $self->{debug};

  # seek ($fh, 0, 0);

  my %cleavages;
  my %dgpi;

  while (<$fh>) {
    my $line = $_;

    if ($line =~ /near/) {
      $line =~ /.+near\s(\d+)\s.+/;
      my $position = $1;
      
      if ($debug) {
	print STDERR "position: $position\n";
      }
      
      $dgpi{position} = $position;

      next;
    }
    
    if ($line =~ /\(score=/) {
      $line =~ /.+at\s(\d+)\s\(score=([^\)]+)\).+/;
      my $site = $1;
      my $potential_cleavage = $2;
      
      $cleavages{$site} = $potential_cleavage;
      
      if ($debug) {
	print STDERR "potential_cleavage: $potential_cleavage at $site\n";
      }
      
      next;
    }

    if ($line =~ /best/) {
      $line =~ /.+best\scleavage\ssite\sis\s(\d+)/;
      my $best_position = $1;

      my $end = $best_position+1;
      my $coord  = "$best_position,$end";

      # get the score
      my $score = $cleavages{$best_position};

      $dgpi{coords}   = $coord;
      $dgpi{cleavage} = $score;

      if ($debug) {
	print STDERR "coord: $coord\n";
	print STDERR "score: $score\n";
      }

      next;
    }
    
    if ($line =~ /is\sGPI-anchored/) {
      
      my $is_GPIanchored = 1;
      
      if ($debug) {
	print STDERR "it has a GPI-anchored signal\n";
      }

      # check next line

      $line = <$fh>;

      if ($line =~ /no\scleavage\ssite\sdetected/i) {
        if ($debug) {
  	  print STDERR "no cleavage site!\n";
        }
        $dgpi{has_cleavage_site} = 0;
      }
      elsif ($line =~ /there\sis\sa\spotential\scleavage\ssite/i) {
 	if ($debug) {
  	  print STDERR "predicted cleavage site\n";
 	}
	$dgpi{has_cleavage_site} = 1;
      }

      return \%dgpi;
    }
  }

  if ($debug) {
    print STDERR "no GPI-anchored site predicted\n";
  }
  
  return undef;

}

sub close_fh {
  my $self = shift;
  my $fh = $self->{fh};
  $fh->close;
}

1;
