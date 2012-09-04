
=head1 NAME - Bio::EnsEMBL::DBSQL::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql'
        );

    $gene_adaptor = $db->get_GeneAdaptor();

    $gene = $gene_adaptor()->fetch_by_stable_id($stable_id);

    $slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end('X', 1, 10000);

=head1 DESCRIPTION


Formerly this class provided database connectivity and a means to retrieve
object adaptors.  This class is now provided for convenience and backwards
compatibility, but delegates its connection responsibilities to the
DBConnection class (no longer inherited from) and its object adaptor
retrieval to the static Bio::EnsEMBL::Registry.


=head1 CONTACT

Post questions to the EnsEMBL development list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

#package Bio::EnsEMBL::DBSQL::DBAdaptor;
package Bio::PSU::Chado;

use vars qw(@ISA $AUTOLOAD);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::SeqRegionCache;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

#use Bio::EnsEMBL::Utils::Exception qw(warning throw  deprecate stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;
my $reg = "Bio::EnsEMBL::Registry";

#CHADO stuff
use DBI;
use Bio::GMOD::Config;
use Bio::GMOD::DB::Config;

=head2 new

  Arg [-DNADB]: (optional) Bio::EnsEMBL::DBSQL::DBAdaptor DNADB 
               All sequence, assembly, contig information etc, will be
               retrieved from this database instead.
  Arg [..]   : Other args are passed to superclass
               Bio::EnsEMBL::DBSQL::DBConnection
  Example    : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    -user   => 'root',
						    -dbname => 'pog',
						    -host   => 'caldy',
						    -driver => 'mysql' );
  Exmaple2   : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                    -species => 'Homo_sapiens',
                                                    -group   => 'core'
						    -user   => 'root',
						    -dbname => 'pog',
						    -host   => 'caldy',
						    -driver => 'mysql');
  Description: Constructor for DBAdaptor.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub new2 {
	my ( $class, @args ) = @_;

	my $self = {};
	bless $self, $class;

	my ( $species, $group, $con, $dnadb ) =
	  rearrange( [qw(SPECIES GROUP DBCONN DNADB)], @args );

	if ( defined($con) ) {
		$self->dbc($con);
	}
	else {
		$self->dbc( new Bio::EnsEMBL::DBSQL::DBConnection(@args) );
	}

	if ( defined($species) ) {
		$self->species($species);
	}

	if ( defined($group) ) {
		$self->group($group);
	}

	$self = Bio::EnsEMBL::Utils::ConfigRegistry::gen_load($self);

	if ( defined $dnadb ) {
		$self->dnadb($dnadb);
	}

	return $self;
}

sub new {
	my ( $class, @args ) = @_;

	my $self = {};
	bless $self, $class;

	my ($db) = rearrange( [qw(DB)], @args );

	my $gmod_conf =
	  $ENV{'GMOD_ROOT'}
	  ? Bio::GMOD::Config->new( $ENV{'GMOD_ROOT'} )
	  : Bio::GMOD::Config->new();

	my $db_conf = Bio::GMOD::DB::Config->new( $gmod_conf, $db );
	my $dbh     = $db_conf->dbh;

	if ( defined $dbh ) {
		$self->dbh($dbh);
	}

	#return $dbh or die "Unable to connect to database\n";
	return $self;
}

sub get_seq {
	my ( $self, $name ) = @_;
	my $sth = $self->dbh->prepare(
		"select residues
                             from feature
                             where uniquename = ?"
	);
	$sth->execute($name);

	my ($residues) = $sth->fetchrow_array();

	my $query = Bio::Seq->new(
		-display_id => $name,
		-seq        => $residues
	);
	return $query;

}

sub create_fasta_file {
	my ( $self, $org, $type ) = @_;
	my $file = "/tmp/$type.$org.fasta";

	my $sth = $self->dbh->prepare(
		"SELECT uniquename, residues
                             FROM feature
                             WHERE type_id = ( SELECT cvterm_id 
                                               FROM cvterm
                                               WHERE name = ? and cv_id = ( SELECT cv_id 
                                                                            FROM cv 
                                                                            WHERE name='Sequence Ontology' )
                                             )
                             AND organism_id = ( SELECT organism_id 
                                                 FROM organism
                                                 WHERE common_name = ?)"
	);
	$sth->execute( $type, $org );

	my $seqio = Bio::SeqIO->new(
		-format => 'fasta',
		-file   => ">$file"
	);

	foreach my $row ( @{ $sth->fetchall_arrayref() } ) {
		my $seq = Bio::Seq->new(
			-display_id => $$row[0],
			-seq        => $$row[1]
		);
		$seqio->write_seq($seq);
	}
	return $file;
}

sub get_type_id {
	my ( $self, $type ) = @_;

	my $sth = $self->dbh->prepare(
		"SELECT cvterm_id
		FROM cvterm
		WHERE name = ?
		AND cv_id = ( SELECT cv_id
			FROM cv
			WHERE name='Sequence Ontology')
		"
	);
	$sth->execute($type);
	print STDERR "Getting type_id for $type\n";
	my ($id) = $sth->fetchrow_array();
	return $id;
}

sub get_org_id {
	my ( $self, $org ) = @_;

	my $sth = $self->dbh->prepare(
		"SELECT organism_id
		FROM organism
		WHERE common_name = ?
		"
	);
	$sth->execute($org);
	my ($id) = $sth->fetchrow_array();
	return $id;
}

#
# deal with 'set_' and 'get_' methods
#
sub AUTOLOAD {
	my ( $self, $newval ) = @_;

	my ($key) = ( $AUTOLOAD =~ /.*::(\w+)/ );    # We don't want Feature::fmin

	if ( defined $newval ) {
		$self->{$key} = $newval;
	}
	else {

		#$DB::single = 1;
		if ( defined( $self->{$key} ) ) {
			return $self->{$key};
		}
		else {
			return undef;
		}
	}
}


sub DESTROY { }    # required due to AUTOLOAD

1;
