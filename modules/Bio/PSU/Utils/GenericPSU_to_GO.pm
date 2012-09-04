package Bio::PSU::Utils::GenericPSU_to_GO;

use strict;
use Data::Dumper;
use Carp;

use Bio::PSU::Utils::FeatureIdHelper::BioPSUHelper;

my $debug = 1;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name => value pairs (odd number supplied)";
  }

  my $self = {
	      @_,
	     };
  
  bless ($self, $class);
  return $self;
}

sub embl_parsing {
  my $self = shift;
  my ($seq, $tigr_format, $db_name, $taxon_id, $assigned_by, $regex) = @_;
  
  my $id_helper =
    Bio::PSU::Utils::FeatureIdHelper::BioPSUHelper->new(-selection_regex => $regex,
							-verbose => 0);
  
  my $time;
  
  {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
    
    $year += 1900;
    
    $time = sprintf "%s%2.2d%2.2d", $year, $mon, $mday;
  }
  
  if (defined $seq) {
    my @features = ();
    
    for my $feature ($seq->features) {
      next unless (($feature->key () eq "CDS") || ($feature->key () eq "mRNA"));
      
      if ((not $tigr_format) && $feature->qexists ('pseudo')) {
	if ($debug) {
	  print STDERR "pseudogene!\n";
	}
	next;
      }

      # processing systematic id

      my $systematic_id = undef;
      
      eval {
	$systematic_id = $id_helper->get_systematic_id($feature);
      };
      
      if ($@) {
	my $loc = join ",", map {
	  $_->{start} . ".." . $_->{end}
	} $feature->ranges ();
	
	print STDERR "Error in CDS at: $loc\n";
	print STDERR "A Systematic Id is missing !!\n";
	next;
      }

      if (defined $systematic_id) {
	if ($debug) {
	  print STDERR "\nprocessing gene, $systematic_id...\n";
	}
      }
      
      # processing standard name

      my $primary_name = undef;

      eval {
	$primary_name = $id_helper->get_primary_name($feature);
      };

      if (not defined $primary_name) {
	if ($debug) {
	  print STDERR "no gene primary name at: " . get_range_string($feature) . "\n";
	  print STDERR "the primary name will be given as the systematic name\n";
	}
	$primary_name = $systematic_id;
      }

      # processing product

      my $product = getProduct ($feature);

      # Return an array of GO hashes
      my ($GO_hashes) = getGOTerms ($feature);

      # processing synonyms

      my $synonym = getSynonyms ($feature);
	
      if ($debug) {
	print STDERR "\nprocessing GO qualifiers...\n\n";
      }

      # processing GO qualifiers ...

      if ($debug) {
	print STDERR "go qualifiers: " . @$GO_hashes . "\n";
      }

      foreach my $GO_qualifier (@$GO_hashes) {
	my $ontology      = $GO_qualifier->{aspect};
	my $ontology_id   = mapOntologyAspect ($ontology);
	my $go_id         = $GO_qualifier->{id};
	my $go_term       = $GO_qualifier->{term};
	my $evidence_code = $GO_qualifier->{ev_code};
	my $dbxref        = $GO_qualifier->{dbxref};
	my $with_or_from  = "";
	defined $GO_qualifier->{with} and $with_or_from = $GO_qualifier->{with};
	my $date;
	if (defined $GO_qualifier->{date}) {
          $date           = $GO_qualifier->{date};
	}
	else {
	  # check time is defined !
	  $date = $time;
	}
	
	##################
	# NOT information
	##################

	# to modify accordingly...

	my $not_flag = 0;
	
	if ($go_id =~ s/^NOT //) {
	  $not_flag = 1;
	}
	
	my $db_object_id     = $systematic_id;
	my $db_object_symbol = undef;
	if (defined $primary_name) {
	  $db_object_symbol = $primary_name;
	}
	else {
	  $db_object_symbol = $systematic_id;
	}
	
	print "$db_name\t$db_object_id\t$db_object_symbol\t\t$go_id\t$dbxref\t$evidence_code\t$with_or_from\t$ontology_id\t$product\t$synonym\tgene\ttaxon:$taxon_id\t$date\t$assigned_by\n";
      }
    }
  } else {
    if ($debug) {
      die "End of file while reading tab file\n";
    }
  }

}

#########
# END
##

sub get_range_string
{
  my ($feature) = shift;
  my @ranges = map {
    $_->{start} . ".." . $_->{end}
  } $feature->ranges ();

  return join ",", @ranges;
}


sub getSystematicId {
  my ($feature) = @_;
  
  my $systematic_id = undef;
  if ($feature->systematic_id ()) {
    $systematic_id = $feature->systematic_id ();
  }
  elsif ($feature->temporary_systematic_id ()) {
    $systematic_id = $feature->temporary_systematic_id ();
  }
  else {
    my @genes = $feature->gene ();
    $systematic_id = $genes[0];
  }

  return $systematic_id;
}

sub getPrimaryName {
  my ($feature) = @_;

  my $primary_name = undef;
  
  if ($feature->primary_name) {
    $primary_name = $feature->primary_name;
  }
  else {
    my @genes = $feature->gene ();
    if (@genes>1) {
      $primary_name = $genes[1];
    }
  }

  return $primary_name;
}


sub getProduct {
  my ($feature) = @_;

  my $product = "";
  
  if ($feature->qexists ('product')) {
    my @products = $feature->qvalues ('product');
    
    # concatenate all the product qualifiers
    # $product = join "|", @products;
    
    # NOW just pick the first one !!!
    $product = $products[0];
  }
  else {
    if ($debug) {
      print STDERR "no /product at: " . get_range_string ($feature);
    }
  }
  
  return $product;
  
}

sub getSynonyms {
  my ($feature) = @_;
  my $synonym = "";

  if ($feature->qexists ('synonym')) {
    my @synonyms = $feature->qvalues('synonym');
    $synonym = join "|", @synonyms;
  }
  else {
    if ($debug) {
      print STDERR "no /synonym at: " . get_range_string ($feature);
    }
  }

  return $synonym;

}

sub getGOTerms {
  my ($feature) = @_;
  
  my @GO_hashes = ();
  
  if (($feature->qexists ("GO_component")) || ($feature->qexists ("GO_function")) || ($feature->qexists ("GO_process"))) {
 
    # old syntax

    if ($feature->qexists ("GO_component")) {
      my @values = $feature->qvalues ('GO_component');
      foreach my $go_component (@values) {
	my %GO_hash;
	if ($go_component =~ /^(GO:\d+)\s*\(([^\)]+)\);\s*([^;]+);\s*([^;]+);/) {
	  $GO_hash{aspect}  = "component";
	  $GO_hash{id}      = $1;
	  $GO_hash{term}    = $2;
	  $GO_hash{ev_code} = $3;
	  my $db_xref       = $4;

	  if ($db_xref =~ /\s*([^\(]+)\(.+\)/) {

	    print STDERR "getting rid of the source and keep just the dbxref\n";

	    $db_xref = $1;

	    print STDERR "dbxref: $db_xref\n";

	  }

	  $GO_hash{dbxref}  = $4;
	}
	else {
	  print STDERR "ERROR: can't parse GO component qualifier, $go_component!\n";
	}
      }
    }

    if ($feature->qexists ("GO_process")) {
      my @values = $feature->qvalues ('GO_process');
      foreach my $go_process (@values) {
	my %GO_hash;
	if ($go_process =~ /^(GO:\d+)\s*\(([^\)]+)\);\s*([^;]+);\s*([^;]+);/) {
	  $GO_hash{aspect}  = "process";
	  $GO_hash{id}      = $1;
	  $GO_hash{term}    = $2;
	  $GO_hash{ev_code} = $3;
	  my $db_xref       = $4;

	  if ($db_xref =~ /\s*([^\(]+)\(.+\)/) {

	    print STDERR "getting rid of the source and keep just the dbxref\n";

	    $db_xref = $1;

	    print STDERR "dbxref: $db_xref\n";

	  }

	  $GO_hash{dbxref}  = $4;
	}
	else {
	  print STDERR "ERROR: can't parse GO process qualifier, $go_process!\n";
	}
      }
    }

    if ($feature->qexists ("GO_function")) {
      my @values = $feature->qvalues ('GO_function');
      foreach my $go_function (@values) {
	my %GO_hash;
	if ($go_function =~ /^(GO:\d+)\s*\(([^\)]+)\);\s*([^;]+);\s*([^;]+);/) {
	  $GO_hash{aspect}  = "function";
	  $GO_hash{id}      = $1;
	  $GO_hash{term}    = $2;
	  $GO_hash{ev_code} = $3;
	  my $db_xref       = $4;

	  if ($db_xref =~ /\s*([^\(]+)\(.+\)/) {

	    print STDERR "getting rid of the source and keep just the dbxref\n";

	    $db_xref = $1;

	    print STDERR "dbxref: $db_xref\n";

	  }

	  $GO_hash{dbxref}  = $4;
	}
	else {
	  print STDERR "ERROR: can't parse GO function qualifier, $go_function!\n";
	}
      }
    }

  }
  elsif ($feature->qexists ("GO")) {
    my @GO_qualifiers = $feature->qvalues ("GO");
    
    foreach my $GO_qualifier (@GO_qualifiers) {

      my %GO_hash;

      if ($GO_qualifier =~ /aspect=([^;]+);/) {
        $GO_hash{aspect} = $1;
      }
      if ($GO_qualifier =~ /GOid=([^;]+);/) {
        $GO_hash{id} = $1;
      }
      if ($GO_qualifier =~ /term=([^;]+);/) {
        $GO_hash{term} = $1;
      }
      if ($GO_qualifier =~ /evidence=([^;]+);/) {
        $GO_hash{ev_code} = $1;
      }
      if ($GO_qualifier =~ /db_xref=([^;]+);/) {
        $GO_hash{dbxref} = $1;
      }
      if ($GO_qualifier =~ /with=([^;]+);/) {
        $GO_hash{with} = $1;
      }
      if ($GO_qualifier =~ /date=([^;]+);/) {
        $GO_hash{date} = $1;
      }
      
      push (@GO_hashes, \%GO_hash);
    }
  }
  else {
    # print STDERR "no GO qualifiers...\n";
  }
  
  return (\@GO_hashes);
  
}

sub mapOntologyAspect {
  my ($ontology) = @_;
  my $ontology_id;
  
  # here is the mapping between "component" and "C" etc.
  
  if ($ontology =~ /component/i) {
    $ontology_id = "C";
  }
  elsif ($ontology =~ /function/i) {
    $ontology_id = "F";
  }
  elsif ($ontology =~ /process/i) {
    $ontology_id = "P";
  }
  else {
    print STDERR "ERROR: don't know about this ontology, $ontology!!\n";
  }
  
  return $ontology_id;
}

1;
