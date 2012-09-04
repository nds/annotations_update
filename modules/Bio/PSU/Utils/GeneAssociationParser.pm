package Bio::PSU::Utils::GeneAssociationParser;

use strict;
use Data::Dumper;
use Carp;

# eval "require Bio::PSU::Utils::GeneAssociationEntry";
use Bio::PSU::Utils::GeneAssociationEntry;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name => value pairs (odd number supplied)";
  }

  my $self = {
	      geneAssociationFile => undef,
	      debug => 0,
	      @_,
	     };
  
  bless ($self, $class);
  return $self;
}

#########
# set the estwise output
# it can be a file or a directory containing a set of estwise output files
##

sub setGeneAssociationFile {
  my $self = shift;
  my ($ga_file) = @_;
  $self->{geneAssociationFile} = $ga_file;
}

sub getGeneAssociationFile {
  my $self = shift;
  return $self->{geneAssociationFile};
}


sub parse {
  my $self = shift;
  my $debug = $self->{debug};
  my @entries = ();
  my %entries;
  
  if (-d $self->{geneAssociationFile}) {
    my @files = getFiles ($self->{geneAssociationFile});
    foreach my $file (@files) {
      my @entries_tmp = $self->parseOneFile ($file);
      push (@entries, @entries_tmp);
    }
  }
  else {
    @entries = $self->parseOneFile;
  }
  
  if ($debug) {
    print STDERR "ga entries: " . Dumper (@entries) . "\n";
  }
  
  # generate the entries hashtable
  
  if ($debug) {
    print STDERR "generating gene association entries hashtable...\n";
  }

  foreach my $entry (@entries) {
    if ($debug) {
      print STDERR "db_object_id: " . $entry->getDbObjectId . "\n";
    }
    if (defined $entry->getDbObjectId) {
      my $db_object_id = $entry->getDbObjectId;
      
      if ($debug) {
	print STDERR "db_object_id: $db_object_id\n";
      }
      
      my @entriesPerGene = ();
      if (defined ($entries{$db_object_id})) {
	my @entriesPerGene = @{$entries{$db_object_id}};
	push (@entriesPerGene, $entry);
	$entries{$db_object_id} = \@entriesPerGene;
      }
      else {
	push (@entriesPerGene, $entry);
	$entries{$db_object_id} = \@entriesPerGene;
      }
    }
  }
  
  return %entries;
}


sub parseOneFile {
  my $self = shift;
  my $debug = $self->{debug};
  my $file = "";
  if (@_) {
    $file = shift;
  }
  else {
    $file = $self->{geneAssociationFile};
  }
  my @entries = ();

  if ($debug) {
    print STDERR "parsing ga file, $file...\n";
  }

  open GA, "< $file" or die "can't open file, $!\n";

  while (<GA>) {
    my $line = $_;
    
    # parsing line
    $line =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.*)/;

    my $db = $1;
    my $db_object_id = $2;
    my $db_object_symbol = $3;
    my $not = $4;
    my $go_id = $5;
    my $db_ref = $6;
    my $evidence = $7;
    my $with_or_from = $8;
    my $aspect = $9;
    my $db_object_name = $10;
    my $synonym = $11;
    my $db_object_type = $12;
    my $taxon = $13;
    my $date = $14;
    my $assigned_by = $15;

    if ($debug) {
      print STDERR "line: $line\n";
      print STDERR "\t* db: $db\n";
      print STDERR "\t* db_object_id: $db_object_id\n";
      print STDERR "\t* db_object_symbol: $db_object_symbol\n";
      print STDERR "\t* not: $not\n";
      print STDERR "\t* go_id: $go_id\n";
      print STDERR "\t* db_ref: $db_ref\n";
      print STDERR "\t* evidence: $evidence\n";
      print STDERR "\t* with_or_from: $with_or_from\n";
      print STDERR "\t* aspect: $aspect\n";
      print STDERR "\t* db_object_name: $db_object_name\n";
      print STDERR "\t* synonym: $synonym\n";
      print STDERR "\t* $db_object_type: $db_object_type\n";
      print STDERR "\t* taxon: $taxon\n";
      print STDERR "\t* date: $date\n";
    }
    
    my $ga_entry = Bio::PSU::Utils::GeneAssociationEntry->new (
							       db => $db,
							       db_object_id => $db_object_id,
							       db_object_symbol => $db_object_symbol,
							       not => $not,
							       go_id => $go_id,
							       db_ref => $db_ref,
							       evidence => $evidence,
							       with_or_from => $with_or_from,
							       aspect => $aspect,
							       db_object_name => $db_object_name,
							       synonym => $synonym,
							       db_object_type => $db_object_type,
							       taxon => $taxon,
							       date => $date,
							       assigned_by => $assigned_by,
							      );
    
    push (@entries, $ga_entry);
  }

  if ($debug) {
    print STDERR "parsing done.\n";
  }

  close GA;
  return @entries;
}

1;