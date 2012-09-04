package Bio::PSU::MultiFastaUtil;
use strict;
use Carp;
use Bio::SeqIO;

sub new {

	my $class = shift;
	my $self = {};
	my $reference_file = shift;
	bless ($self, $class);
	$reference_file or croak "Reference file not provided";
	$self->reference_file($reference_file);
	my $offsets = $self->_read_file($reference_file);
	$self->_offsets($offsets);
	return $self;
}


sub _read_file {
	my $self = shift;
	my $reference_file = shift;
	my $total_length = 0;
	my %offsets;
	my $in = Bio::SeqIO->new(-file=>$reference_file, -format=>'fasta');
	
	while (my $seq = $in->next_seq){
		my $id = $seq->id;
		my $length = length($seq->seq);
		$offsets{$id} = $total_length;
		$total_length += $length;
	}
	
	return \%offsets;
}

sub reference_file {
	my $self = shift;
	$self->{reference_file} = shift if @_;
	return $self->{reference_file};
}

sub _offsets {
	my $self = shift;
	$self->{offsets} = shift if @_;
	return $self->{offsets};
}

sub get_offset {
	my $self = shift;
	my $id = shift;
	$self->_offsets and return $self->_offsets->{$id};
	return undef;
}

sub get_ids {
	my $self = shift;
	$self->_offsets and return keys %{$self->_offsets};
}	return {};

1;
