package Bio::PSU::Utils::BlastParsing;

use strict;
use Data::Dumper;
use Carp;

# Bio modules
# PSU
use Bio::PSU::SearchFactory;
use Bio::PSU::IO::BufferFH;

# Bioperl
# use Bio::SearchIO;

# eval "require Bio::PSU::Utils::BlastHit";
use Bio::PSU::Utils::BlastHit;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name => value pairs (odd number supplied)";
  }

  my $self = {
	      blastOutput  => undef,
	      max_hits => 20,
	      @_,
	     };
  bless ($self, $class);
  return $self;
}

#########
# set the blast output
# it can be a file or a directory containing a set of blast output files
##

sub setBlastOutput {
  my $self = shift;
  my ($blastOutput) = @_;
  $self->{blastOutput} = $blastOutput;
}

sub getBlastOutput {
  my $self = shift;
  return $self->{blastOutput};
}

sub setMaxHits {
  my $self = shift;
  my ($max_hits) = @_;
  $self->{max_hits} = $max_hits;
}

sub getMaxHits {
  my $self = shift;
  return $self->{max_hits};
}

sub parse {
  my $self = shift;
  my %results;
  my @results = ();

  if (-d $self->{blastOutput}) {
    my @files = getFiles ($self->{blastOutput});
    foreach my $file (@files) {
      my @results_tmp = $self->parseOneFile ($file);
      push (@results, @results_tmp);
    }
  }
  else {
    @results = $self->parseOneFile;
  }

  # generate the hashtable

  foreach my $result (@results) {
    # print STDERR "parsing new result...\n";

    ##########
    # PSU
    ##
    my $q_acc = $result->getQueryAccession;

    ##########
    # Bioperl
    ##

    # name or accession
    # my $q_acc = $result->query_accession;

    # get rid of the description
    my @tmp = split (/\s/, $q_acc);
    $q_acc = $tmp[0];
    my @resultsPerEST = ();
    if (defined ($results{$q_acc})) {
      my @resultsPerEST = @{$results{$q_acc}};
      push (@resultsPerEST, $result);
      $results{$q_acc} = \@resultsPerEST;
    }
    else {
      push (@resultsPerEST, $result);
      $results{$q_acc} = \@resultsPerEST;
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
    $file = $self->{blastOutput};
  }
  my $max_hits = $self->{max_hits};
  my @hits = ();

  print STDERR "parsing a new file, $file...\n";

  ##############
  # PSU
  ##

  my $blast = Bio::PSU::SearchFactory->make (
					     -file      => $file,
					     -program   => 'blast'
					    );
  
  while (my $result = $blast->next_result) {
    my $q_acc = $result->q_id;

    # print STDERR "Processing blastx results for q_acc, $q_acc...\n";
    
    my $search_type = $result->type;
    my $i = 0;
    while ((my $hit_psu=$result->next_hit) && ($i<$max_hits)) {
      my $s_acc  = $hit_psu->s_id;
      my $s_desc = $hit_psu->s_desc;
      my $s_length = $hit_psu->s_len;
      my $gene_name = undef;

      # parsing the subject accession number retruned by the PSU modules.
      
      # print STDERR "hiting s_acc, $s_acc\n";

      if ($s_acc =~ /\|/) {
	my @tmp = split ('\|', $s_acc);
	
	# print STDERR "tmp: @tmp\n";
	# re. Droso, there are 4 accession numbers - keep the second one
	# because this is the flybase accession number
	$s_acc = $tmp[1];
	
	# the first one is the gene name
	$gene_name = $tmp[0];
      }

      # Just report the best HSP...
      # Apparently it is alread rporting the best HSP first...

      # while (my $hsp = $hit_psu->next_hsp) {
      if (my $hsp = $hit_psu->next_hsp) {
	my $score    = $hsp->score;
	my $percent  = $hsp->percent;
	my $evalue   = $hsp->expect;
	my $q_strand = $hsp->q_strand;
	my ($q_start, $q_end);
	if ($q_strand == 1) {
	  $q_start  = $hsp->q_begin;
	  $q_end    = $hsp->q_end;
	}
	else {
	  # Reverse start end stop to make sure stop > start !!
	  $q_start  = $hsp->q_end;
	  $q_end    = $hsp->q_begin;
	}
	my $q_frame  = $hsp->q_frame;
	my $s_strand = $hsp->s_strand;
        my $s_start  = $hsp->s_begin;
	my $s_end    = $hsp->s_end;

	my $hit = Bio::PSU::Utils::BlastHit->new (
						  search_type => $search_type,
						  query_accession => $q_acc,
						  q_strand => $q_strand,
						  q_start  => $q_start,
						  q_end    => $q_end,
						  q_frame  => $q_frame,
						  s_strand => $s_strand,
						  s_start  => $s_start,
						  s_end    => $s_end,
						  subject_accession => $s_acc,
						  subject_description => $s_desc,
						  subject_gene_name => $gene_name,
						  subject_length => $s_length,
						  score => $score,
						  percent => $percent,
						  evalue => $evalue,
						 );
	push (@hits, $hit);
      }
      $i++;
    }
  }
  ####

  ###############
  # Bioperl
  ##

  # my $searchio = new Bio::SearchIO (-format => 'blast',
  #				    -file   => $file
  #				   );
  #  while ( my $result = $searchio->next_result() ) {
  #    push (@results, $result);
  #  }
  ###

  # @hits = sort_by_query_accession (@hits);
  
  return @hits;
}

sub sort_by_query_accession {
  my (@hits) = @_;
  my @new_hits = ();

  while (@hits) {
    my ($index, $hit) = getMinEST (@hits);
    push (@new_hits, $hit);
    splice (@hits, $index, 1);
  }

  return @new_hits;
}

sub getMinEST {
  my (@hits) = @_;
  my $index = 0;
  my $hit = undef;
  my $query_acc_min = undef;

  my $i = 0;
  while ($i<@hits) {
    my $hit_tmp = $hits[$i];
    if (not defined $hit) {
      $hit = $hit_tmp;
      $query_acc_min = $hit_tmp->getQueryAccession;
    }
    else {
      my $query_acc_tmp = $hit_tmp->getQueryAccession;
      if (($query_acc_tmp cmp $query_acc_min) == -1) {
	$hit = $hit_tmp;
	$query_acc_min = $query_acc_tmp;
	$index = $i;
      }
    }
    $i++;
  }

  return ($index, $hit);
}

sub getFiles {
  my ($directory) = @_;
  my @files = ();
  
  opendir (THISDIR, $directory) or die "can not read this directory, $!\n";
  @files = map {$directory."/".$_} grep {($_ =~ /.*\.res$/)} readdir THISDIR;
  closedir THISDIR;
  
  # print STDERR "files: @files\n";

  return @files;
}

1;
