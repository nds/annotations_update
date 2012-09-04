package Bio::PSU::Utils::BlastHit;

use strict;
use Carp;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name=>value pairs (odd number supplied)";
  }

  my $self     = {
		  search_type => undef,
		  query_accession => undef,
		  q_start  => undef,
		  q_end    => undef,
		  q_strand => undef,
		  q_frame  => undef,
                  s_start  => undef,
                  s_end    => undef,
		  s_strand => undef,
		  subject_accession => undef,
		  subject_description => undef,
		  subject_gene_name => undef,
		  subject_length => undef,
		  score => undef,
		  percent => undef,
		  evalue => undef,
		  @_,
		 };
  bless ($self, $class);
  return $self;
}

sub getSearchType {
  my $self = shift;
  return $self->{search_type};
}

sub getQueryAccession {
  my $self = shift;
  return $self->{query_accession};
}

sub getQueryStart {
  my $self = shift;
  return $self->{q_start};
}

sub getQueryEnd {
  my $self = shift;
  return $self->{q_end};
}

sub getQueryStrand {
  my $self = shift;
  return $self->{q_strand};
}

sub getQueryFrame {
  my $self = shift;
  return $self->{q_frame};
}

sub getSubjectAccession {
  my $self = shift;
  return $self->{subject_accession};
}

sub getSubjectStart {
  my $self = shift;
  return $self->{s_start};
}

sub getSubjectEnd {
  my $self = shift;
  return $self->{s_end};
}

sub getSubjectStrand {
  my $self = shift;
  return $self->{s_strand};
}

sub getSubjectDescription {
  my $self = shift;
  return $self->{subject_description};
}

sub getSubjectGeneName {
  my $self = shift;
  return $self->{subject_gene_name};
}

sub getSubjectLength {
  my $self = shift;
  return $self->{subject_length};
}

sub getScore {
  my $self = shift;
  return $self->{score};
}

sub getPercent {
  my $self = shift;
  return $self->{percent};
}

sub getEValue {
  my $self = shift;
  return $self->{evalue};
}

1;
