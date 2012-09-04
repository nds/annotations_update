#!/bin/perl -w

# Copyright (C) 2005 Genome Research Limited. All Rights Reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=head1 Name

FastFasta.pm: some methods for accessing sequences and information about
              sequences within a fasta file quickly using little memory

=head1 Author

Sendu Bala (email: bix@sendu.me.uk)

=cut

package Bio::PSU::FastFasta;

use strict;
use warnings;
use Carp;
use vars qw($VERSION);

$VERSION = 1.00;


# constructor
sub new {
    my $class = shift;
    
    my $self = { @_ };
    
    defined $self->{-infile} or defined $self->{-outfile} or die "-infile or -outfile must be defined\n";
    defined $self->{-infile} and my $in_file = $self->{-infile};
    defined $self->{-outfile} and my $out_file = $self->{-outfile};
    $in_file and $out_file and ($in_file eq $out_file) and die "-infile cannot be the same as -outfile\n";
    
    # open the file
    my ($in_fh, $out_fh);
    if ($in_file) {
        open($in_fh, $in_file) or croak "can't read input file '$in_file'";
        $self->{in_fh} = $in_fh;
    }
    if ($out_file) {
        open($out_fh, ">$out_file") or croak "can't write output file '$out_file'";
        $self->{out_fh} = $out_fh;
    }
    
    unless ($self->{-fast} || ! $in_file) {
        my ($by_index, $by_id, $lengths_by_index, $lengths_by_id, $smallest, $largest, $total, $descriptions_by_id) = _starts($in_fh);
        $self->{by_index} = $by_index;
        $self->{by_id} = $by_id;
        $self->{descriptions} = $descriptions_by_id;
        $self->{lengths_by_index} = $lengths_by_index;
        $self->{lengths_by_id} = $lengths_by_id;
        $self->{smallest} = $smallest;
        $self->{largest} = $largest;
        $self->{total} = $total;
        $self->{average} = sprintf("%.2f", ($total / @{$by_index}));
    }
    
    return bless($self, $class);
}

# indexes the start position and lengths of all the seqeuences in the file by id and order
sub _starts {
    my ($fh) = @_;
    
    # even though trying to be memory efficent, Tie::IxHash is too slow
    my @by_index;
    my %by_id;
    my %descriptions_by_id;
    my %lengths_by_index;
    my %lengths_by_id;
    
    my $smallest = 99999999999999999999;
    my $largest = 0;
    my $total = 0;
    
    local $/ = "\n>";
    while (! eof($fh)) {
        my $start = tell($fh);
        
        my $seq = <$fh>;
        
        $seq =~ s/^>?(\S+)(.*)$//m;
        my $id = $1;
        $descriptions_by_id{$id} = $2;
        $seq =~ s/\n//g;
        $seq =~ s/>$//;
        my $length = length($seq);
        
        push(@by_index, $start);
        $by_id{$id} = $#by_index;
        $lengths_by_index{$#by_index} = length($seq);
        exists $lengths_by_id{$id} and warn "warning: multiple sequences in the file share the same id '$id', only the length of the last can be accessed\n";
        $lengths_by_id{$id} = $length;
        
        if ($length < $smallest) {
            $smallest = $length;
        }
        if ($length > $largest) {
            $largest = $length;
        }
        $total += $length;
    }
    
    seek($fh, $by_index[0], 0);
    
    return (\@by_index, \%by_id, \%lengths_by_index, \%lengths_by_id, $smallest, $largest, $total, \%descriptions_by_id);
}

# given a filehandle returns (id, sequence) for the sequence at the current seek position
sub _get_seq {
    my ($fh) = @_;
    my ($id, $seq);
    
    # setting $/ here every time slows us down significantly, but not much more than setting $/ globally in fast mode and detecting fast mode here
    local $/ = "\n>";
    
    $seq = <$fh>;
    $seq or return ();
    $seq =~ s/^>?(\S+).*$//m;
    $id = $1;
    $seq =~ s/\n//g;
    $seq =~ s/>$//;
    
    return ($id, $seq);
}

# Get the next (id, sequence) from the fasta file
sub next_seq {
    my $self = shift;
    my $fh = $self->{in_fh};
    
    return _get_seq($fh);
}

# Get the number of sequences in the file
sub num_of_seqs {
    my $self = shift;
    
    return scalar(@{$self->{by_index}});
}

# Get the total bases in the file
sub total {
    my $self = shift;
    
    return $self->{total};
}

# Get the smallest sequence size
sub smallest {
    my $self = shift;
    
    return $self->{smallest};
}

# Get the largest sequence size
sub largest {
    my $self = shift;
    
    return $self->{largest};
}

# Get the average sequence size
sub average {
    my $self = shift;
    
    return $self->{average};
}

# get the description of a sequence (everything after the id on an id line)
sub desc {
    my ($self, $id) = @_;
    $id or return;
    
    return $self->{descriptions}->{$id} if exists $self->{descriptions}->{$id};
}

# Given (many) id(s) or number(s) (to ask for the nth sequence in the file), returns a reference to an id->sequence hash
sub get_seqs {
    my ($self, $request, @requests) = @_;
    $request or return;
    
    # id or number?
    my $id_mode;
    if ($request =~ /^\d+$/) {
        exists $self->{by_id}->{$request} ? ($id_mode = 1) : (exists $self->{by_index}->[$request - 1] or return undef);
    }
    else {
        exists $self->{by_id}->{$request} ? ($id_mode = 1) : return undef;
    }
    
    unshift(@requests, $request);
    
    # store the current position within the file (incase the user is currently looping through the file with next_seq)
    my $fh = $self->{in_fh};
    my $old_position = tell($fh);
    
    # get the sequence(s)
    my %seqs;
    foreach my $request (@requests) {
        $id_mode ? seek($fh, $self->{by_index}->[$self->{by_id}->{$request}], 0) : seek($fh, $self->{by_index}->[$request - 1], 0);
        my ($id, $seq) = _get_seq($fh);
        $seqs{$id} = $seq;
    }
    
    # restore the old position within the file
    seek($fh, $old_position, 0);
    
    return \%seqs;
}

# Get the length of a sequence given its id or a number (to ask for the nth sequence in the file)
sub seq_length {
    my ($self, $request) = @_;
    $request or carp "FastaFasta: you must supply an id to seq_length";
    
    exists $self->{lengths_by_id}->{$request} ? (return $self->{lengths_by_id}->{$request}) : (exists $self->{lengths_by_index}->{$request} ? return $self->{lengths_by_index}->{$request} : return undef);
}

# Get the lengths of all sequences in the file, results returned as a hash reference, hash has id => length pairs
sub seq_lengths {
    my $self = shift;
    
    return $self->{lengths_by_id};
}

# Get the id of the first sequence matching perfectly the given string of letters
sub identify {
    my ($self, $seq) = @_;
    
    # store the current position within the file (incase the user is currently looping through the file with next_seq)
    my $fh = $self->{in_fh};
    my $old_position = tell($fh);
    
    # go through the file from the beginning and look for given sequence
    my $id;
    while (my ($this_id, $this_seq) = $self->next_seq()) {
        if ($this_seq eq $seq) {
            $id = $this_id;
            last;
        }
    }
    
    # restore the old position within the file
    seek($fh, $old_position, 0);
    
    $id ? return $id : return undef;
}

# Get all the ids (as an array) of the sequences containing a perfect full-length match to the given string of letters
sub find_matches {
    my ($self, $seq) = @_;
    
    # store the current position within the file (incase the user is currently looping through the file with next_seq)
    my $fh = $self->{in_fh};
    my $old_position = tell($fh);
    
    # go through the file from the beginning and look for given sequence
    my @ids;
    while (my ($this_id, $this_seq) = $self->next_seq()) {
        if (index($this_seq, $seq) > -1) {
            push(@ids, $this_id);
        }
    }
    
    # restore the old position within the file
    seek($fh, $old_position, 0);
    
    @ids > 0 ? return @ids : return ();
}

# Get the %hash{id}->{start, end} for the sequences containing a near-perfect near-full-length match to the given string of letters
sub find_partials {
    
}

1;
