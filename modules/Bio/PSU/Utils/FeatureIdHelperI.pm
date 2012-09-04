
package Bio::PSU::Utils::FeatureIdHelperI;

use strict;
use vars qw(@ISA);

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

sub selection_regex
{
    my $self = shift;
    $self->throw_not_implemented;
}

sub get_systematic_id
{
    my $self = shift;
    $self->throw_not_implemented;
}

sub get_primary_name
{
    my $self = shift;
    $self->throw_not_implemented;
}

sub get_synonym
{
    my $self = shift;
    $self->throw_not_implemented;
}

1;
