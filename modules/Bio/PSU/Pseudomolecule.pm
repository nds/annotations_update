package Bio::PSU::Pseudomolecule;
use strict;
use Carp;
use Bio::SeqIO;
use ContigPosition;

=head1 NAME
Bio::PSU::Pseudomolecule - Perl module to represent a pseudomolecule

=head1  SYNOPSIS

my $ps = Bio::PSU::Pseudomolecule->new(id=>$id, offset=>$file_offset, reference_length=>$length, reference_offset=> $reference_offset);


=head1 DESCRIPTION

A Pseudomolecule object holds multiple ContigPosition objects and has methods for producing the pseudomolecule sequence,  feature tables and crunch.

Contigs can be added to the pseudomolecule by giving it the details of an alignment using the add_contig_by_alignment method.

A single object should represent the ordering against a single reference molecule. Only one contig placement can exist for a specific contig ID

=head1 AUTHOR

Nick Peters (np1)

=head1 METHODS

=cut 


our $DEBUG = 0;


=head2 new

 Title   : new
 Usage   : my $ps = Bio::PSU::Pseudomolecule->new(id=>$id, offset=>$file_offset, reference_length=>$length, reference_offset=> $reference_offset);
 Function: Parse the next result from the data stream
 Returns : L<Bio::PSU::Pseudomolecule>
 Args    : Minimum of id, reference_length

=cut 

sub new {

	my $class = shift;
	my $self = {};
	@_ % 2 and croak "Uneven parameters list";
	my (%parameters) = @_;
	$parameters{id} or croak "Must provide a sequence ID";
	$parameters{reference_length} or croak "Must provide a reference sequence length";
	bless ($self, $class);
	$self->id($parameters{id});
	$self->reference_length($parameters{reference_length});
	$self->reference_offset(0);
	$parameters{reference_offset} and $self->reference_offset($parameters{reference_offset});
	$parameters{offset} ? $self->offset($parameters{offset}) :  $self->offset(0);
	return $self;
	
}


=head2 id

 Title   : id
 Usage   : my $id = $pseudomolecule->id or $pseudomolecule->id('myID');
 Function: Getter/Setter for the pseudomolecule ID
 Returns : string 
 Args    : None or ID to set

=cut 


sub id {
	my $self = shift;
	my $id = shift;
	$id and $self->{id} = $id;
	return $self->{id};
}

=head2 reference_length

 Title   : reference_length
 Usage   : my $id = $pseudomolecule->reference_length or $pseudomolecule->reference_length(1000);
 Function: Getter/Setter for the pseudomolecule reference_length
 Returns : number 
 Args    : None or reference_length to set

=cut 

sub reference_length {
	my $self = shift;
	my $id = shift;
	$id and $self->{reference_length} = $id;
	return $self->{reference_length};
}

=head2 reference_offset

 Title   : reference_offset
 Usage   : my $id = $pseudomolecule->reference_offset or $pseudomolecule->reference_offset(1000);
 Function: Getter/Setter for the pseudomolecule reference_offset - this is the position of the reference molecule in a multi fasta file
 Returns : number 
 Args    : None or reference_offset to set

=cut 

sub reference_offset {
	my $self = shift;
	my $id = shift;
	$id and $self->{reference_offset} = $id;
	return $self->{reference_offset};
}

=head2 add_contig

 Title   : add_contig
 Usage   : $seudomolecule->new($id, $method, $sequence, $start, $end, $strand);
 Function: Parse the next result from the data stream
 Returns : undef
 Args    : Contig id, placement method (e.g. blast), contig sequence,  start, end, strand

=cut 

sub add_contig {
	my $self = shift;
	my $id = shift;
	my $method = shift;
	my $sequence = shift;
	my $start = shift;
	my $end = shift;
	my $strand = shift;
	$self->add_CongigPosition( Bio::PSU::ContigPosition->new ( id=>$id, method=>$method, sequence=>$sequence, start=>$start, end=>$end, strand=>$strand));
}


=head2 offset

 Title   : offset
 Usage   : my $id = $pseudomolecule->offset or $pseudomolecule->offset(1000);
 Function: Getter/Setter for the pseudomolecule offset - this is the position of the pseudomolecule in a multi fasta file
 Returns : number 
 Args    : None or offset to set

=cut 

sub offset {
	my $self = shift;
	my $offset = shift;
	defined $offset and $self->{offset} = $offset;
	return $self->{offset};
}

=head2 add_contig_by_alignment

 Title   : add_contig_by_alignment
 Usage   : $pseudomolecule->add_contig_by_alignment($id, $method, $sequence, $query_start, $query_end, $reference_start, $reference_end, $strand);
 Function: Add a new ContigPosition based on an alignment.   This returns the ContigPosition and a list of subsumed ContigPositions.  If the contig position being added
 		lies within a currently existing contig position there will be a single ID returned matching the added ID
 Returns : <Bio::PSU::ContigPosition>,  Array reference to a list of the IDs of subsumed ContigPositions
 Args    : id, method, sequence, query start, query end, reference start, reference end, strand

=cut 

sub add_contig_by_alignment {

	my $self = shift;
	my $id = shift;
	my $method = shift;
	my $sequence = shift;
	my $query_start = shift;
	my $query_end = shift;
	my $reference_start = shift;
	my $reference_end = shift;
	my $strand = shift;
	
	
	#need to work out the real bounds of the contig
	
	#$query_start > $query_end and ($query_start,$query_end) = ($query_end, $query_start); 

	$DEBUG and print STDERR "IN:\nREFSTART: $reference_start	REFEND: $reference_end\n";
	$DEBUG and print STDERR "IN:\nQUERYSTART: $query_start	QUERYEND: $query_end\n";
	
	if ($strand == 1){
		$reference_start = $reference_start - ($query_start - 1);
	}
	else{
		$reference_start = $reference_start - (CORE::length($sequence) - $query_end);
	}
	
	$reference_end = $reference_start + (CORE::length($sequence) -1 );
	
	$DEBUG and print STDERR "OUT:\nREFSTART: $reference_start	REFEND: $reference_end\n";
	$DEBUG and print STDERR "OUT:\nQUERYSTART: $query_start	QUERYEND: $query_end\n";
	
	
	#print "$reference_end  $reference_start " . CORE::length($sequence) . "\n";
	
	my ($cp, $subsumed) = $self->add_ContigPosition( Bio::PSU::ContigPosition->new ( id=>$id, method=>$method, sequence=>$sequence, start=>$reference_start, end=>$reference_end, strand=>$strand));
	

	$cp or return undef;
	
	$cp->reference_location ($reference_start, $reference_end, 1);
	
	return ($cp, $subsumed);
}

=head2 add_ContigPosition

 Title   : add_ContigPosition
 Usage   : $pseudomolecule->add_ContigPosition($contig_position);
 Function: Add a new ContigPosition.   This returns the ContigPosition and a list of subsumed ContigPositions.  If the contig position being added
 		lies within a currently existing contig position there will be a single ID returned matching the added ID
 Returns : <Bio::PSU::ContigPosition>,  Array reference to a list of the IDs of subsumed ContigPositions
 Args    : <Bio::PSU::ContigPosition>

=cut 

sub add_ContigPosition {

	my $self = shift;
	my $position = shift;
	$position->isa('Bio::PSU::ContigPosition') or croak "argument to addContig position must be a ContigPosition";
	my $id = $position->id();
	#Check whether it subsumes or is subsumed by another ContigPosition
	
	my @subsumed = $self->_subsumes($position);
	
	
	#if our CP to be added is inside another don't add
	if (@subsumed == 1 && $subsumed[0] eq $id){
		return (undef, \@subsumed);
	}
	foreach  my $id_sub (@subsumed) {
		$self->remove_ContigPosition($id_sub);
	}
	
	$self->{position}->{$id} = $position;
	return ($position, \@subsumed);
}

=head2 get_ContigPosition

 Title   : get_ContigPosition
 Usage   : my $cp = $pseudomolecule->get_ContigPosition($id);
 Function: Retrieve a ContigPosition by its ID
 Returns : L<Bio::PSU::ContigPosition>
 Args    : contig ID

=cut 

sub get_ContigPosition {

	my $self = shift;
	my $id = shift;
	return $self->{position}->{$id};
}

=head2 get_ContigPosition_ids

 Title   : get_ContigPosition_ids
 Usage   : my @ids = $pseudomolecule->get_ContigPosition_ids();
 Function: Retrieve a list of ContigPosition IDs
 Returns : Array of IDs
 Args    : None

=cut

sub get_ContigPosition_ids {

	my $self = shift;
	my @ids = sort  { $self->get_ContigPosition($a)->start <=> $self->get_ContigPosition($b)->start }  (keys %{$self->{position}});
	return @ids;
}


=head2 relative_start

 Title   : relative_start
 Usage   : my $start  = $pseudomolecule->relative_start();
 Function: Returns the start of the lowest starting ContigPosition - NB this may be less than 0
 Returns : Number
 Args    : None

=cut 

sub relative_start {

	my $self = shift;
	my @ids = $self->get_ContigPosition_ids;
	return $ids[0]->start;
}

=head2 relative_end

 Title   : relative_end
 Usage   : my $start  = $pseudomolecule->relative_end();
 Function: Returns the end of the highest ending ContigPosition 
 Returns : Number
 Args    : None

=cut 

sub relative_end {

	my $self = shift;
	my @ids = sort  { $self->get_ContigPosition($a)->end <=> $self->get_ContigPosition($b)->end }  (keys %{$self->{position}});
	return $ids[-1]->end;
}



=head2 remove_ContigPosition

 Title   : remove_ContigPosition
 Usage   : $pseudomolecule->remove_ContigPosition($id);
 Function: Remove a ContigPosition by its ID
 Returns : undef
 Args    : contig ID

=cut 


sub remove_ContigPosition {

	my $self = shift;
	my $id = shift;
	delete $self->{position}->{$id};
}

=head2 sequence

 Title   : sequence
 Usage   : my $ps = Pseudomolecule->sequence();
 Function: Returns the sequence of the pseudomolecule, padded with Ns for internal gaps
 Returns : String
 Args    : None

=cut 

sub sequence {
	my $self = shift;
	
	$DEBUG and print STDERR "Writing sequence...\n";
	
	#sort by start position - need to put everything in an array
	
	my @ids = $self->get_ContigPosition_ids;
	
	my $output_sequence = "";
	
	my $previous_end = 0;
	
	
	foreach my $id (@ids){
		$DEBUG and print STDERR "ID $id\n";
		my $position = $self->get_ContigPosition($id);
		my $start = $position->start;
		my $end = $position->end;	
		my $strand = $position->strand;
		my $sequence;
				
		#Pad with Ns if needed
		if ($previous_end && $start > $previous_end + 1){
		
			$output_sequence .= 'N' x ($start - $previous_end - 1);
		}
		
		if ($strand == 1){
			$sequence = $position->sequence;
			
		}	
		else {
			$sequence = $position->complementary_sequence;
		}
		
		$output_sequence .= $sequence;
		
		$previous_end = $end;
		
		$DEBUG and print STDERR "finished contig\n";
	}
	
	
	#write end section if 
	#if ($previous_end < $self->length){
	#
	#	$output_sequence .= 'N' x  ($self->length - $previous_end);
	#}
	
	return $output_sequence;
}

=head2 to_FTstring

 Title   : to_FTstring
 Usage   : my $pseudomolecule->to_FTstring();
 Function: Returns the feature table representation of the pseudomolecule to match the sequence produced by the sequence method
 Returns : String
 Args    : None

=cut 

sub to_FTstring {

	my $self = shift;
	my $ft_string = $self->_build_molecule('ft');
	
}

=head2 to_crunch

 Title   : to_crunch
 Usage   : my $pseudomolecule->to_crunch();
 Function: Returns the crunch representation of the pseudomolecule/reference sequence to match the sequence produced by the sequence method
 Returns : String
 Args    : None

=cut 

sub to_crunch {

	my $self = shift;
	my $ft_string = $self->_build_molecule('crunch');	
}

=head2 _build_molecule

 Title   : _build_molecule
 Usage   : $pseudomolecule->_build_molecule('ft');
 Function: Private method to build a molecule in the same way the sequence method does.  Handles offsets.  Used by the to_crunch and to_FTString methods
 Returns : String
 Args    : output type (currently 'ft' or 'crunch')

=cut 

sub _build_molecule {

	my $self = shift;
	my $output_type = shift;
	$output_type or croak "Must provide an output type";
	my $return_string = '';
	
	my $processed = 0;
	#we pass artemis colour codes through to the FT method
	my @colours = (3,4);
	
	my @ids = $self->get_ContigPosition_ids;
	my $lowest_cp  = $self->get_ContigPosition($ids[0]);
	
	my $offset = 0;
	my $file_offset = 0;
	$self->offset and $file_offset = $self->offset;
	
	#allow for the fact that we might get negative values - create an offset to make it start at 1 and everything else relative
	if ($lowest_cp->start < 0){
	
		$offset = abs($lowest_cp->start)   + 1;
	}
	#or we may want a negative offset cause the first contig isn't at 1
	elsif ($lowest_cp->start > 1){
	
		$offset = 0 - $lowest_cp->start + 1;
	}
	
	$DEBUG and print STDERR "OFFSET $offset\n";
	
	my $previous_end = 0;
	
	foreach my $id ( @ids ){
	
		#need to increase offset when things overlap 
		my $cp = $self->get_ContigPosition($id);
		
		$DEBUG and print STDERR join (' ',  ($id, $cp->start, $cp->end, $cp->strand)) . "\n"; 
		
		if ($previous_end &&  $cp->start <= $previous_end){
		
			$offset += (($previous_end - $cp->start) + 1);
		}

		if( $output_type eq 'crunch'){
		
			$return_string .= $cp->to_crunch(($offset + $file_offset), $self->reference_length, $self->reference_offset) . "\n";
		}
	       elsif ($output_type eq 'ft'){
	       		my  $colour = $colours[$processed % 2];
			$return_string .= $cp->to_FTstring($offset + $file_offset, $colour );
		}
		else{
			croak "I don't know how to write $output_type";
		}
		
		
		$previous_end = $cp->end;
		$processed++;
	}
	
	return $return_string;
}


=head2 overlapping_location

 Title   : overlapping_location
 Usage   : $pseudomolecule->overlapping_location($location);
 Function: Pase the next result from the data stream
 Returns : Array ref to  a list of ContigPositions overlapping the Location object
 Args    : <Bio::Location::Simple>

=cut  

sub overlapping_location {

	my $self = shift;
	my $location = shift;
	my $overlapping_array = [];
	
	
	
	foreach my $id ( $self->get_ContigPosition_ids ){
	
		my $cp = $self->get_ContigPosition($id);
		if ($location->overlaps($cp->location)){
		
			push ( @{$overlapping_array},  $cp);
		}
	}
	
	return $overlapping_array;
	
}


=head2 _subsumes

 Title   : _subsumes
 Usage   : $pseudomolecule->_subsumes($contig_position);
 Function: Private method.  Checks to see if a ContigPosition to be added subsumes or is subsumed by any existing ContigPositions
 Returns : Array of ids of subsumed ContigPositions
 Args    : <Bio::PSU::ContigPosition>

=cut 

sub _subsumes {

	my $self = shift;
	my $cp = shift;
	my @subsumed_ids = ();

	foreach my $id ( $self->get_ContigPosition_ids ){
	
		my $current_cp = $self->get_ContigPosition($id);
		#our added CP contains a pre-existing CP
		if ($cp->location->contains($current_cp->location)){
			
			$DEBUG and print STDERR $cp->id .  " overlaps $id\n";
			push (@subsumed_ids, $id);
		}
		#a pre-existing one subsumes our added CP
		elsif ($current_cp->location->contains($cp->location)){
			
			$DEBUG and print STDERR "Added CP is inside an existing CP\n"; 
			return ($cp->id);
		}
		
	}
	
	
	
	
	return @subsumed_ids;
	
}


=head2 length

 Title   : length
 Usage   : my $length =  $pseudomolecule->length();
 Function: Returns the length of the pseudomolecule
 Returns : Number
 Args    : None

=cut 

sub length {

	my $self = shift;
	return length ( $self->sequence );
}




1;
