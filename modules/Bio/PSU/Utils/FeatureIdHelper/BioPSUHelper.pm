
package Bio::PSU::Utils::FeatureIdHelper::BioPSUHelper;

use strict;
use vars qw(@ISA);

use Bio::Root::Root;

use Bio::PSU::Feature;
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

sub _restrict_re
{
  return '^' . $_[0] . '$';
}

sub get_systematic_id
{
    my ($self, $obj) = @_;

    my $id = undef;

    eval
    {
        $id = $self->_select_single_value($obj,
                                          "systematic_id",
                                          "temporary_systematic_id",
                                          "locus_tag");
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
            $id = $self->_select_matching_value($obj, "gene");
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
        my $fth = Bio::PSU::IO::FTHandler->new;
        my $text = $fth->make_feature_block($obj);
        print STDERR "problem with this feature:\n\n$text\n";

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
                                          "primary_name",
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
            $id = $self->_select_matching_value($obj, "gene");
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
        $self->throw("Failed to recover primary_name\n");
    }

    return $id;
}

sub get_synonyms
{
    my ($self, $obj) = @_;

    my @synonyms;
    my $primary_name;

    if ($obj->isa('Bio::PSU::Feature'))
    {
        $primary_name = $self->get_primary_name($obj);

        if ($obj->qexists(_restrict_re 'synonym'))
        {
            @synonyms = $obj->qvalues(_restrict_re "synonym");
        }
        elsif ($obj->qexists(_restrict_re 'gene'))
        {
            @synonyms = $obj->qvalues(_restrict_re "gene");
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
    my ($self, $obj, $primary_tag, $secondary_tag, $tertiary_tag) = @_;

    my $value;


    if ($obj->isa('Bio::PSU::Feature'))
    {
        if ($obj->qexists(_restrict_re $primary_tag))
        {

            my @vals = $obj->qvalues(_restrict_re $primary_tag);

            $value = $vals[0];
            if (scalar @vals > 1)
            {
                $self->warn("Feature at " . $obj->start() .
                            " has more than one $primary_tag tag - picking first: $value\n");
            }

        }
        elsif (defined $secondary_tag)
        {
            if ($obj->qexists(_restrict_re $secondary_tag))
            {
                my @vals = $obj->qvalues(_restrict_re $secondary_tag);

                $value = $vals[0];
                if (scalar @vals > 1)
                {
                    $self->warn("Feature at " . $obj->start() .
                                " has more than one $secondary_tag tag - picking first: $value\n");
                }
            }
	    elsif (defined $tertiary_tag)
	    {
        	if ($obj->qexists(_restrict_re $tertiary_tag))
		{
         		my @vals = $obj->qvalues(_restrict_re $tertiary_tag);

	        	$value = $vals[0];
			if (scalar @vals > 1)
 		        {
			    $self->warn("Feature at " . $obj->start() .
					" has more than one $tertiary_tag tag - picking first: $value\n");
			}
		    }
		else
		{
		    $self->throw("Feature at " . $obj->start() .
                             " has neither $primary_tag nor $secondary_tag nor $tertiary_tag tags\n");
		}
	    }
	}
	else
	{
	    $self->throw("Feature at " . $obj->start() .
                         " has no $primary_tag tag\n");
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

    if ($obj->isa('Bio::PSU::Feature'))
    {
        if ($obj->qexists(_restrict_re $primary_tag))
        {
            my @vals = $obj->qvalues(_restrict_re $primary_tag);
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
                    $self->warn("More than one qualifier and selection regex doesn't match - picked first: $value\n");
                }
            }
            else
            {
                $value = $vals[0];
                $self->warn("More than one qualifier and no selection regex was supplied - picked first: $value\n");
            }
        }
        else
        {
            $self->throw("Feature at " . $obj->start() .
                         " has no $primary_tag tag\n");
        }
    }
    else
    {
         $self->throw("Argument was not a Bio::SeqFeatureI\n");
    }

    return $value;
}
