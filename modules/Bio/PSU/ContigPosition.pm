package Bio::PSU::ContigPosition;
use strict;
use Carp;
use Bio::Location::Simple;


=head1 NAME
Bio::PSU::ContigPosition - Perl module to represent a Contig positioned on a reference sequence

=head1  SYNOPSIS

my $cp = Bio::PSU::ContigPosition->new( id=>$id, method=>$method, sequence=>$sequence, start=>$start, end=>$end, strand=>$strand));

or 

my $cp = Bio::PSU::ContigPosition->new( id=>$id, method=>$method, sequence=>$sequence, $bio_location_simple));

=head1 DESCRIPTION

Object to contain the sequence, ID, and placement method (e.g. blast alignment) of a contig and it's position relative to a reference sequence.

Has methods to write crunch and feature table representations of itself.

=head1 AUTHOR

Nick Peters (np1)

=head1 METHODS

=cut 



=head2 new

 Title   : new
 Usage   : my $cp = Bio::PSU::ContigPosition->new( id=>$id, method=>$method, sequence=>$sequence, start=>$start, end=>$end, strand=>$strand));
 Function: Create a new Bio::PSU::ContigPosition object
 Returns : <Bio::PSU::ContigPosition>
 Args    : Contig id, method, sequence  and eith a Bio::Location::Simple object or start, end, strand

=cut 


sub new {

	my $class = shift;
	my $self = {};
	@_ % 2 and croak "Uneven parameters list";
	my (%params) = @_;
	$params{id} or croak "Must provide a sequence ID";
	$params{method} or croak "Must provide a location method";
	$params{sequence} or croak "Must provide a sequence";
	bless ($self, $class);
	$self->id($params{id});
	$self->method($params{method});
	if ($params{location}){
		$self->location($params{location});
	}
	elsif ( $params{start} && $params{end} && $params{strand}){
		$self->location($params{start}, $params{end}, $params{strand});
	}	
	else{
		croak " Must provide either (location), (start, end, strand), or an alignment location information";
	}
	length $params{sequence} == $self->location->length or croak "Sequence length does not match location length";
	$self->sequence($params{sequence});
	
	
	
	return $self;
	
}

=head2 id

 Title   : id
 Usage   : my $id = $c_position->id or $c_position->id('myID');
 Function: Getter/Setter for the ContigPosition ID
 Returns : string 
 Args    : None or ID to set

=cut 

sub id {
	my $self = shift;
	my $id = shift;
	$id and $self->{id} = $id;
	return $self->{id};
}

=head2 method

 Title   : method
 Usage   : my $method = $c_position->method or $c_position->method('mymethod');
 Function: Getter/Setter for the ContigPosition method
 Returns : string 
 Args    : None or method to set

=cut 

sub method {
	my $self = shift;
	my $id = shift;
	$id and $self->{method} = $id;
	return $self->{method};
}


=head2 location

 Title   : location
 Usage   : my $location = $c_position->location or $c_position->location($bio_location) or $c_position->location($start, $stop, $strand);
 Function: Getter/Setter for the ContigPosition location
 Returns : <Bio::Location::Simple> 
 Args    : None or location to set

=cut 


sub location {
	my $self = shift;

	if($_[0] && ref($_[0])){
		if ($_[0]->isa('Bio::Location::Simple')){
			$_[0]->to_FTstring;
			$self->{location} = $_[0];
		}
		else{
			croak " Must provide either location or start, end, strand information";	
		}
	}
	elsif (@_ == 3) {
		my ($start, $end, $strand) = @_;
		$self->{location} = Bio::Location::Simple->new(-start => $start, -end => $end, -strand => $strand);
	}
	
	return $self->{location};
}


=head2 reference_location

 Title   : reference_location
 Usage   : my $reference_location = $c_position->reference_location or $c_position->reference_location($bio_reference_location) or $c_position->reference_location($start, $stop, $strand);
 Function: Getter/Setter for the ContigPosition reference_location
 Returns : <Bio::reference_location::Simple> 
 Args    : None or reference_location to set

=cut 

sub reference_location {
	my $self = shift;

	if($_[0] && ref($_[0])){
		if ($_[0]->isa('Bio::Location::Simple')){
			$_[0]->to_FTstring;
			$self->{location} = $_[0];
		}
		else{
			croak " Must provide either location or start, end, strand information";	
		}
	}
	elsif (@_ == 3) {
		my ($start, $end, $strand) = @_;
		$self->{ref_location} = Bio::Location::Simple->new(-start => $start, -end => $end, -strand => $strand);
	}
	
	return $self->{ref_location};
}

=head2 sequence

 Title   : sequence
 Usage   : my $sequence = $c_position->sequence or $c_position->sequence('AAATGCTTTTT');
 Function: Getter/Setter for the ContigPosition sequence
 Returns : string 
 Args    : None or sequence to set

=cut 

sub sequence {

	my $self = shift;
	my $sequence = shift;
	if ($sequence){
		$sequence =~ /[^atgcxn]/i and croak "Sequence must only be DNA (and X or N)";
		$self->{sequence} = $sequence;
	}
	
	return $self->{sequence};
}


=head2 complementary_sequence

 Title   : complementary_sequence
 Usage   : my $complementary_sequence = $c_position->complementary_sequence;
 Function: Getter for the ContigPosition complementary_sequence
 Returns : string 
 Args    : None 

=cut

sub complementary_sequence {

	my $self = shift;
	my $sequence  = $self->{sequence};
	$sequence =~ tr/ATGCatgc/TACGtacg/;
	return reverse($sequence);
}




=head2 start

 Title   : start
 Usage   : my $start = $c_position->start;
 Function: Getter for the ContigPosition start
 Returns : Number 
 Args    : None 

=cut 

sub start {
	my $self = shift;
	return $self->location->start;
}


=head2 end

 Title   : end
 Usage   : my $end = $c_position->end;
 Function: Getter for the ContigPosition end
 Returns : Number 
 Args    : None 

=cut 

sub end {
	my $self = shift;
	return $self->location->end;
}

=head2 strand

 Title   : strand
 Usage   : my $strand = $c_position->strand;
 Function: Getter for the ContigPosition strand
 Returns : Number 
 Args    : None 

=cut 

sub strand {
	my $self = shift;
	return $self->location->strand;
}


=head2 to_FTstring

 Title   : to_FTstring
 Usage   : my $c_position->to_FTstring($offset, $colour);
 Function: Returns the feature table representation of the pseudomolecule to match the sequence produced by the sequence method.
 Returns : String
 Args    : Offset to increase coordinates by, colour code to use in file

=cut 

sub to_FTstring {
	my $self = shift;
	my $offset = shift;
	my $colour = shift;
	$colour or $colour = 1;
	my $location; 
	my $loc_string;
	if ($offset){
		my $start = $self->location->start + $offset;
		my $end = $self->location->end + $offset;
		my $strand = $self->location->strand;
		$location = Bio::Location::Simple->new (-start=>$start,-end=>$end,-strand=>$strand);
		$loc_string = $location->to_FTstring;
	}else{
	 	$loc_string = $self->location->to_FTstring;
	}
	
	my $id = $self->id;
	my $method = $self->method;
	
return <<EOF;
FT   CONTIG          $loc_string
FT                   /systematic_id="$id"
FT                   /method="$method"
FT                   /colour="$colour"
EOF

}


=head2 to_crunch

 Title   : to_crunch
 Usage   : my $c_position->to_crunch();
 Function: Returns the crunch representation of the ContigPosition mapped to the reference
 Returns : String
 Args    : offset to increase pseudomolecul coords by, reference sequence length, reference offset coords

=cut 

sub to_crunch {

	my $self = shift;
	my $offset = shift;
	my $reference_length = shift;
	my $reference_offset = shift;
	$reference_offset or $reference_offset = 0;
	$reference_length or croak "Must provide reference length";
	$offset or $offset = 0;
	my $ref_location = $self->reference_location;
	$ref_location or return undef;
	my $start = $self->start + $offset;
	my $end = $self->end + $offset;
	$self->strand or ($start, $end) = ($end, $start);
	my $ref_start = $ref_location->start;
	$ref_start < 0 and $ref_start = 1;
	$ref_start += $reference_offset;
	my $ref_end = $ref_location->end;
	$ref_end > $reference_length and $ref_end = $reference_length;
	$ref_end += $reference_offset;
	
	return join (' ', (1, , 100, $start, $end, $self->id, $ref_start, $ref_end, 'unknown', 'NONE'));
	
}



1;
