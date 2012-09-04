package Bio::PSU::Utils::GeneAssociationEntry;

use strict;
use Carp;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name=>value pairs (odd number supplied)";
  }

  my $self = {
	      db => undef,
	      db_object_id => undef,
	      db_object_symbol => undef,
	      not => undef,
	      go_id => undef,
	      db_ref => undef,
	      evidence => undef,
	      with_or_from => undef,
	      aspect => undef,
	      db_object_name => undef,
	      synonym => undef,
	      db_object_type => undef,
	      taxon => undef,
	      date => undef,
	      assigned_by => undef,
	      @_,
	     };
  bless ($self, $class);
  return $self;
}

sub setDb {
  my $self = shift;
  my ($db) = @_;
  $self->{db} = $db;
}

sub getDb {
  my $self = shift;
  return $self->{db};
}

sub getDbObjectId {
  my $self = shift;
  return $self->{db_object_id};
}

sub setDbObjectId {
  my $self = shift;
  my ($db_object_id) = @_;
  $self->{db_object_id} = $db_object_id;
}

sub getDbObjectSymbol {
  my $self = shift;
  return $self->{db_object_symbol};
}

sub setDbObjectSymbol {
  my $self = shift;
  my ($db_object_symbol) = @_;
  $self->{db_object_symbol} = $db_object_symbol;
}

sub setNot {
  my $self = shift;
  my ($not) = @_;
  $self->{not} = $not;
}

sub getNot {
  my $self = shift;
  return $self->{not};
}

sub setGOid {
  my $self = shift;
  my ($go_id) = @_;
  $self->{go_id} = $go_id;
}

sub getGOid {
  my $self = shift;
  return $self->{go_id};
}

sub setDbReference {
  my $self = shift;
  my ($db_ref) = @_;
  $self->{db_ref} = $db_ref;
}

sub getDbReference {
  my $self = shift;
  return $self->{db_ref};
}

sub setEvidence {
  my $self = shift;
  my ($evidence) = @_;
  $self->{evidence} = $evidence;
}

sub getEvidence {
  my $self = shift;
  return $self->{evidence};
}

sub getWithOrFrom {
  my $self = shift;
  return $self->{with_or_from};
}

sub setWithOrFrom {
  my $self = shift;
  my ($with_or_from) = @_;
  $self->{with_or_from} = $with_or_from;
}

sub setAspect {
  my $self = shift;
  my ($aspect) = @_;
  $self->{aspect} = $aspect;
}

sub getAspect {
  my $self = shift;
  return $self->{aspect};
}

sub setDbObjectName {
  my $self = shift;
  my ($db_object_name) = @_;
  $self->{db_object_name} = $db_object_name;
}

sub getDbObjectName {
  my $self = shift;
  return $self->{db_object_name};
}

sub setSynonym {
  my $self = shift;
  my ($synonym) = @_;
  $self->{synonym} = $synonym;
}

sub getSynonym {
  my $self = shift;
  return $self->{synonym};
}

sub setDbObjectType {
  my $self = shift;
  my ($db_object_type) = @_;
  $self->{db_object_type} = $db_object_type;
}

sub getDbObjectType {
  my $self = shift;
  return $self->{db_object_type};
}

sub setTaxon {
  my $self = shift;
  my ($taxon) = @_;
  $self->{taxon} = $taxon;
}

sub getTaxon {
  my $self = shift;
  return $self->{taxon};
}

sub setDate {
  my $self = shift;
  my ($date) = @_;
  $self->{date} = $date;
}

sub getDate {
  my $self = shift;
  return $self->{date};
}

sub setAssignedBy {
  my $self = shift;
  my ($assigned_by) = @_;
  $self->{assigned_by} = $assigned_by;
}

sub getAssignedBy {
  my $self = shift;
  return $self->{assigned_by};
}

1;