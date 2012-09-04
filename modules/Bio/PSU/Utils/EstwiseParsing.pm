package Bio::PSU::Utils::EstwiseParsing;

use strict;
use Data::Dumper;
use Carp;

eval "require Bio::PSU::Utils::EstwiseHit";
# use Bio::PSU::Utils::EstwiseHit;

my $debug = 0;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name => value pairs (odd number supplied)";
  }

  my $self = {
	      estwiseOutput  => undef,
	      pfamMapping => undef,
	      iproMapping => undef,
	      goMapping   => undef,
	      @_,
	     };

  ###########################################
  # cache the files content into arrays !!
  ##
  
  # Pfam
  
  open MAPPING, "< $self->{pfamMapping}" or die "can't open file, $!\n";
  my @pfamMappingSet = ();
  while (<MAPPING>) {
    my $line = $_;
    chomp ($line);
    push (@pfamMappingSet, $line);
  }
  close MAPPING;
  $self->{pfamMappingSet} = [@pfamMappingSet];

  # interpro

  open MAPPING, "< $self->{iproMapping}" or die "can't open file, $!\n";
  my @iproMappingSet = ();
  while (<MAPPING>) {
    my $line = $_;
    chomp ($line);
    push (@iproMappingSet, $line);
  }
  close MAPPING;
  $self->{iproMappingSet} = [@iproMappingSet];

  # go

  open MAPPING, "< $self->{goMapping}" or die "can't open file, $!\n";
  my @goMappingSet = ();
  while (<MAPPING>) {
    my $line = $_;
    chomp ($line);
    push (@goMappingSet, $line);
  }
  close MAPPING;
  $self->{goMappingSet} = [@goMappingSet];

  #print STDERR "Pfam set Dump:\n" . Dumper (@pfamMappingSet) . "\n###########################################\n";
  #print STDERR "ipro set Dump:\n" . Dumper (@iproMappingSet) . "\n###########################################\n";
  #print STDERR "go set Dump:\n" . Dumper (@goMappingSet) . "\n###########################################\n";

  bless ($self, $class);
  return $self;
}

#########
# set the estwise output
# it can be a file or a directory containing a set of estwise output files
##

sub setEstwiseOutput {
  my $self = shift;
  my ($estwiseOutput) = @_;
  $self->{estwiseOutput} = $estwiseOutput;
}

sub getEstwiseOutput {
  my $self = shift;
  return $self->{estwiseOutput};
}

sub mapPfam_name_accession {
  my $self = shift;
  my ($pfam_name) = @_;
  my @pfamMappingSet = @{$self->{pfamMappingSet}};
  my $pfam_acc = undef;

  # print STDERR "mapping pfam name and accession...\n" ;

  foreach my $line (@pfamMappingSet) {
    chomp ($line);
    if ($line =~ /.+$pfam_name$/) {
      $line =~ /\s*(\w+)\s+.+/;
      $pfam_acc = $1;
      last;
    }
  }

  return $pfam_acc;
}

sub mapPfam2Ipro {
  my $self = shift;
  my ($pfam_acc) = @_;
  my @iproMappingSet = @{$self->{iproMappingSet}};
  my $ipro_acc  = "undefined";
  my $ipro_name = "undefined";

  # print STDERR "mapping pfam 2 interpro accession numbers... - lines: " . @iproMappingSet . "\n";

  foreach my $line (@iproMappingSet) {
    chomp ($line);
    if ($line =~ /.+$pfam_acc.+/) {
      $line =~ /.+(IPR\d+)\s+(.+)/;
      $ipro_acc  = $1;
      $ipro_name = $2;
      last;
    }
  }

  return ($ipro_acc, $ipro_name);
}

sub mapIpro2GO {
  my $self = shift;
  my ($ipro_acc) = @_;
  my @goMappingSet = @{$self->{goMappingSet}};
  my @go_list = ();
  
  # print STDERR "mapping interpro 2 GO accession numbers...\n";

  foreach my $line (@goMappingSet) {
    chomp ($line);
    if ($line =~ /.+$ipro_acc.+/) {
      $line =~ /.+GO:(.+)\s\;\s(GO:\d+)/;
      my $go_term = $1;
      my $go_acc  = $2;

      my %go = (
		go_term => $go_term,
		go_acc  => $go_acc,
	       );
      push (@go_list, \%go);
    }
  }

  return @go_list;
}


sub parse {
  my $self = shift;
  my %results;
  my @results = ();

  if (-d $self->{estwiseOutput}) {
    my @files = getFiles ($self->{estwiseOutput});
    foreach my $file (@files) {
      my @results_tmp = $self->parseOneFile ($file);
      push (@results, @results_tmp);
    }
  }
  else {
    @results = $self->parseOneFile;
  }

  # generate the %results hashtable

  print STDERR "generating estwise hashtable...\n";

  foreach my $result (@results) {
    # print STDERR "estwise hit: " . Dumper ($result) . "\n";
    if (defined $result->getESTAccession) {
      my $est = $result->getESTAccession;
      my @resultsPerEST = ();
      if (defined ($results{$est})) {
	my @resultsPerEST = @{$results{$est}};
	push (@resultsPerEST, $result);
	$results{$est} = \@resultsPerEST;
      }
      else {
	push (@resultsPerEST, $result);
	$results{$est} = \@resultsPerEST;
      }
    }
  }
  return %results;
}


sub parseOneFile {
  my $self = shift;
  my $file = "";
  if (@_) {
    $file = shift;
  }
  else {
    $file = $self->{estwiseOutput};
  }
  my @results = ();

  print STDERR "parsing file, $file...\n";

  open ESTWISE, "< $file" or die "can't open file, $!\n";

  # activate filter initially
  # no longer filter !!!
  #my $filter = 0;

  while (<ESTWISE>) {
    my $line = $_;
    
    # parsing line
    $line =~ /^Protein\s+(\w+\-*\w*)\s+DNA\s+(\S+)\s+(\S+)\s+(\d+\.\d+)/;
    
    my $est_acc   = $3;
    my $strand    = $2;
    my $pfam_name = $1;
    my $pfam_acc  = undef;
    
    if ($debug) {
      print STDERR "Found Pfam name: $pfam_name\n";
    }
    
    $pfam_acc  = $self->mapPfam_name_accession ($pfam_name);
    my $bits      = $4;    
    my $ipro_acc  = "undefined";
    my $ipro_name = "undefined";
    my @go_accs   = ();
    
    if (defined ($pfam_acc)) {
      ($ipro_acc, $ipro_name)  = $self->mapPfam2Ipro ($pfam_acc);
      @go_accs = $self->mapIpro2GO ($ipro_acc);
    }
    else {
      $pfam_acc = "undefined";
    }

    my $estwise_hit = Bio::PSU::Utils::EstwiseHit->new (
							est_accession => $est_acc,
							strand => $strand,
							pfam_name => $pfam_name,
							pfam_accession => $pfam_acc,
							ipro_accession => $ipro_acc,
							ipro_name => $ipro_name,
							go_accessions => \@go_accs,
							bits => $bits,
						       );
    
    push (@results, $estwise_hit);
  }

  print STDERR "parsing done.\n";

  close ESTWISE;
  return @results;
}

sub sort_by_EST {
  my $self = shift;
  my (@results) = @_;
  my @new_results = ();
  
  while (@results) {
    my ($index, $estwise_hit) = getMinEST (@results);
    push (@new_results, $estwise_hit);
    splice (@results, $index, 1);
  }

  return @new_results;
}

sub getMinEST {
  my (@results) = @_;
  my $index = 0;
  my $estwise_hit = undef;
  my $est_acc_min = undef;

  my $i = 0;
  while ($i<@results) {
    my $estwise_hit_tmp = $results[$i];
    if (not defined $estwise_hit) {
      $estwise_hit = $estwise_hit_tmp;
      $est_acc_min = $estwise_hit_tmp->getESTAccession;
    }
    else {
      my $est_acc_tmp = $estwise_hit_tmp->getESTAccession;
      if (($est_acc_tmp cmp $est_acc_min) == -1) {
	$estwise_hit = $estwise_hit_tmp;
	$est_acc_min = $est_acc_tmp;
	$index = $i;
      }
    }
    $i++;
  }

  return ($index, $estwise_hit);
}

sub getFiles {
  my ($directory) = @_;
  my @files = ();
  
  opendir (THISDIR, $directory) or die "can not read this directory, $!\n";
  @files = map {$directory."/".$_} grep {($_ =~ /.*\.res$/)} readdir THISDIR;
  closedir THISDIR;
  
  return @files;
}
