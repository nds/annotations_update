
package Bio::PSU::Utils::FeatureIdHelper::BioPerlHelper;

use strict;
use vars qw(@ISA);

use Bio::Root::Root;
use Bio::SeqFeatureI;

use Bio::PSU::Utils::FeatureIdHelperI;

@ISA = qw(Bio::Root::Root Bio::PSU::Utils::FeatureIdHelperI);

sub new
{
    my ($proto, @args) = @_;
    my $self = $proto->SUPER::new(@args);

    my ($regex) = $self->_rearrange([qw(SELECTION_REGEX)], @args);
    $self->selection_regex($regex);
    return $self;
}

sub selection_regex
{
    my ($self, $regex) = @_;

    if (defined $regex)
    {
        my $type = ref($regex);

        if ($type)
        {
            if (($type ne 'Regexp') && ($type ne 'ARRAY'))
            {
                $self->throw("Regex must be a scalar, a compiled regex or an array reference.");
            }
        }

        $self->{'selection_regex'} = $regex;
    }

    return $self->{'selection_regex'};
}

sub get_systematic_id
{
    my ($self, $obj) = @_;

    my $id = undef;

    eval
    {

     $id = $self->_select_single_value($obj,
                                       'systematic_id',
                                       'temporary_systematic_id');
    };

    if ($@)
    {
        if ($self->verbose > 0)
        {
            print STDERR "$@\n";
        }
    }

    if (! defined $id)
    {
        eval
        {
            $id = $self->_select_matching_value($obj, 'gene');
        };

        if ($@)
        {
            if ($self->verbose > 0)
            {
                print STDERR "$@\n";
            }
        }
    }

    if (! defined $id)
    {
        warn "Warning: problem with feature at " .
             $obj->location()->to_FTstring () . "\n";
        $self->throw("Failed to recover systematic_id\n");
    }

    return $id;
}

sub get_primary_name
{
    my ($self, $obj) = @_;

    my $id = undef;

    eval
    {
        $id = $self->_select_single_value($obj,
                                          'primary_name',
                                          undef);
    };

    if ($@)
    {
        if ($self->verbose > 0)
        {
            print STDERR "$@\n";
        }
    }

    if (! defined $id)
    {
        eval
        {
            $id = $self->_select_matching_value($obj, 'gene');
        };

        if ($@)
        {
            if ($self->verbose > 0)
            {
                print STDERR "$@\n";
            }
        }
    }

    if (! defined $id)
    {
        $self->throw("Failed to recover primary_name:\n$@\n");
    }

    return $id;
}

sub get_synonyms
{
    my ($self, $obj) = @_;

    my @synonyms;
    my $primary_name;

    if ($obj->isa('Bio::SeqFeatureI'))
    {
        $primary_name = $self->get_primary_name($obj);

        if ($obj->has_tag('synonym'))
        {
            @synonyms = $obj->get_tag_values('synonym');
        }
        elsif ($obj->has_tag('gene'))
        {
            @synonyms = $obj->get_tag_values('gene');
        }
    }
    else
    {
        $self->throw("Argument was not a Bio::SeqFeatureI\n");
    }

    return grep ! /^$primary_name$/, @synonyms;
}

sub _select_single_value
{
    my ($self, $obj, $primary_tag, $secondary_tag) = @_;

    my $value;

    if ($obj->isa('Bio::SeqFeatureI'))
    {
        if ($obj->has_tag($primary_tag))
        {
            my @vals = $obj->get_tag_values($primary_tag);

            if (scalar @vals > 1)
            {
                $self->throw("Feature has more than one $primary_tag tag\n");
            }

            $value = $vals[0];
        }
        elsif (defined $secondary_tag)
        {
            if ($obj->has_tag($secondary_tag))
            {
                my @vals = $obj->get_tag_values($secondary_tag);

                if (scalar @vals > 1)
                {
                    $self->throw("Feature has more than one $secondary_tag tag\n");
                }

                $value = $vals[0];
            }
            else
            {
                $self->throw("Feature has neither $primary_tag nor $secondary_tag tags\n");
            }
        }
        else
        {
            $self->throw("Feature has no $primary_tag tag\n");
        }
    }
    else
    {
        $self->throw("Argument was not a Bio::SeqFeatureI\n");
    }

    return $value;
}

sub _select_matching_value
{
    my ($self, $obj, $primary_tag) = @_;

    my $value;

    if ($obj->isa('Bio::SeqFeatureI'))
    {
        if ($obj->has_tag($primary_tag))
        {
            my @vals = $obj->get_tag_values($primary_tag);
            my $regex = $self->selection_regex;

            if (scalar @vals == 1)
            {
                $value = $vals[0];
            }
            elsif (defined $regex)
            {
                my @regexes = ($regex);

                if (ref($regex) && ref($regex) eq 'ARRAY')
                {
                    @regexes = @$regex;
                }

                foreach my $r (@regexes)
                {
                    my @matched_vals = grep /$r/, @vals;
                    next unless @matched_vals;

                    $value = $matched_vals[0];

                    if (scalar @matched_vals > 1)
                    {
                      $self->warn("Selection regex $r matches more than one value - picking first: $value\n");
                    }
                    last;
                }

                if (!defined $value) {
                    $value = $vals[0];
                    $self->warn("Tag has more than one value and selection regex doesn't match - picked first: $value\n");
                }
            }
            else
            {
              $value = $vals[0];
              $self->warn("Tag has more than one value and no selection regex was supplied - picked first: $value\n");
            }
        }
        else
        {
            $self->throw("Feature has no $primary_tag tag\n");
        }
    }
    else
    {
         $self->throw("Argument was not a Bio::SeqFeatureI\n");
    }

    return $value;
}

1;
