package Bio::PSU::PSUReposUtils;
#
# Hopefully most of the PSURepository code will end up in this module.
#
# At the time of writing, this module is placed here when live;
#    /nfs/pathsoft/prod/genlib/perl/src
# Not ideal, but it means no new path as to be created.
#
# Author: Paul Mooney (pjm@sanger.ac.uk)
#

use Data::Dumper;

my $INSTANCE = 'utlp';
my $USER     = 'psu_repos';
my $PASSWORD = 'psu_repos';
my $DRIVER   = 'Oracle';


################################################################################
#
# new
#
################################################################################
sub new {
    my ($type, $debug) = @_;

    # Message Of The Day can go here
    #print "\n", "x"x80, "\n";
    #print "  The PSU Repository system will be unavailable on Tuesday 7th September\n";
    #print "  between 12:00 and 15:00 for maintenance of the database\n";
    #print "x"x80, "\n\n";


    my $self = {
        'sql_get_via_nickname' => 
            'SELECT NICK_NAME, PATH_DATA_DIR, PROJECT_NAME '.
            'FROM   PATH_RESPOSITORY '.
            'WHERE  NICK_NAME = ?',

        'sql_get_all' => 
            'SELECT NICK_NAME, PATH_DATA_DIR, PROJECT_NAME '.
            'FROM   PATH_RESPOSITORY ',

        'sql_get_via_project_name' => 
            'SELECT NICK_NAME, PATH_DATA_DIR, PROJECT_NAME '.
            'FROM   PATH_RESPOSITORY '.
            'WHERE  PROJECT_NAME = ?',

        'sql_get_via_path_data_dir' => 
            'SELECT NICK_NAME, PATH_DATA_DIR, PROJECT_NAME '.
            'FROM   PATH_RESPOSITORY '.
            'WHERE  PATH_DATA_DIR = ?',

        'sql_insert_entry' => 
            'INSERT INTO PATH_RESPOSITORY '.
            ' (NICK_NAME, PATH_DATA_DIR, PROJECT_NAME) '.
            'VALUES (?, ?, ?)',

        'PATH_DATA_DIR' => "/nfs/pathdata"
    };

    if (defined($debug)) {
        $self->{'debug'} = 1;
    }
    else {
        $self->{'debug'} = 0;
    }

    bless($self, $type);
    return $self;
}

################################################################################
#
# connect_to_DB
#
# This needs to be made less Oracle specific
#
################################################################################
sub connect_to_DB {
    my ($self, $instance, $user, $password, $driver) = @_;

    if (defined $self->{'dbh'}) {
        return $self->{'dbh'};
    }

    $instance = $instance || $INSTANCE;
    $user     = $user     || $USER;
    $password = $password || $PASSWORD;
    $driver   = $driver   || $DRIVER;

    my $dbh = DBI->connect("dbi:Oracle:$instance", 
                           $user,
                           $password,
                           { RaiseError => 1, AutoCommit => 0 },
                           )
        or die "ERROR:$DBI::errstr\n";

    if( !defined $dbh ) {
        die "Could not connect to database";
    }

    $self->{'dbh'} = $dbh;

    return $dbh;
}

################################################################################
#
# disconnect
#
# delete all statement handles and disconnect
#
################################################################################
sub disconnect {
    my ($self) = @_;

    foreach my $sth (keys %{$self}) {
        if ($sth =~ /^sth/) {
            delete $self->{$sth};
        }
    }

    $self->{'dbh'}->disconnect() or die $self->{'dbh'}->errstr;
    delete $self->{'dbh'};
}

################################################################################
#
# commit 
#
################################################################################
sub commit {
    my ($self) = @_;

    $self->{'dbh'}->commit() or die $self->{'dbh'}->errstr;
}

################################################################################
#
# rollback
#
################################################################################
sub rollback {
    my ($self) = @_;

    $self->{'dbh'}->rollback() or die $self->{'dbh'}->errstr;
}


################################################################################
#
# get_data_for_nickname
#
# Executes prepared SQL for $nick_name and returns the row from the database
#
################################################################################
sub get_data_for_nickname {
    my ($self, $nick_name) = @_;

    unless (defined $self->{'dbh'}) {
        $self->connect_to_DB();
    }

    unless (defined $self->{'sth_get_via_nickname'}) {
        $self->{'sth_get_via_nickname'} = 
            $self->{'dbh'}->prepare($self->{'sql_get_via_nickname'});
    }

    if ($self->{'debug'}) {
        print STDERR "get_data_for_nickname(): sth_get_via_nickname = ",
            $self->{'sth_get_via_nickname'},"\n";
        print STDERR "get_data_for_nickname(): nick_name = $nick_name\n";
    }

    $self->{'sth_get_via_nickname'}->execute($nick_name)
        or die $dbh->errstr;

    my $hash = $self->{'sth_get_via_nickname'}->fetchrow_hashref();

    print STDERR "\$hash = ", Data::Dumper::Dumper($hash)
        if $self->{'debug'};

    return $hash;
}

################################################################################
#
# get_all_data
#
# Executes prepared SQL to get all rows from PATH_RESPOSITORY and the hash
# returned is keyed by NICK_NAME
#
################################################################################
sub get_all_data {
    my ($self) = @_;

    unless (defined $self->{'dbh'}) {
        $self->connect_to_DB();
    }

    unless (defined $self->{'sth_get_all_via_nickname'}) {
        $self->{'sth_get_all_via_nickname'} = 
            $self->{'dbh'}->prepare($self->{'sql_get_all'});
    }

    $self->{'sth_get_all_via_nickname'}->execute()
        or die $dbh->errstr;

    my $hash = $self->{'sth_get_all_via_nickname'}->fetchall_hashref('NICK_NAME');

    return $hash;
}

################################################################################
#
# get_data_for_project_name
#
# Executes prepared SQL to get the row from PATH_RESPOSITORY for the project
# with the 
#
################################################################################
sub get_data_for_project_name {
    my ($self, $project_name) = @_;

    unless (defined $self->{'dbh'}) {
        $self->connect_to_DB();
    }

    unless (defined $self->{'sth_get_via_project_name'}) {
        $self->{'sth_get_via_project_name'} = 
            $self->{'dbh'}->prepare($self->{'sql_get_via_project_name'});
    }

    $self->{'sth_get_via_project_name'}->execute($project_name)
        or die $dbh->errstr;

    my $hash = $self->{'sth_get_via_project_name'}->fetchrow_hashref();

    return $hash;
}

################################################################################
#
# get_data_for_path_data_dir
#
# Executes prepared SQL for $nick_name and returns the row from the database
#
################################################################################
sub get_data_for_path_data_dir {
    my ($self, $path_data_dir) = @_;

    unless (defined $self->{'dbh'}) {
        $self->connect_to_DB();
    }

    unless (defined $self->{'sth_get_via_path_data_dir'}) {
        $self->{'sth_get_via_path_data_dir'} = 
            $self->{'dbh'}->prepare($self->{'sql_get_via_path_data_dir'});
    }

    $self->{'sth_get_via_path_data_dir'}->execute($path_data_dir)
        or die $dbh->errstr;

    my $hash = $self->{'sth_get_via_path_data_dir'}->fetchrow_hashref();

    return $hash;
}

################################################################################
#
# get_path_data_base_dir
#
################################################################################
sub get_path_data_base_dir {
    my ($self) = @_;

    return $self->{'PATH_DATA_DIR'};
}



################################################################################
#
# get_path_data_dir
#
################################################################################
sub create_entry {
    my ($self, $nick_name, $path_data_dir, $project_name) = @_;

    print STDERR "create_entry()\n";

    unless (defined $self->{'dbh'}) {
        $self->connect_to_DB();
    }

    unless (defined $self->{'sth_insert_entry'}) {
        $self->{'sth_insert_entry'} = 
            $self->{'dbh'}->prepare($self->{'sql_insert_entry'});
    }

    $self->{'sth_insert_entry'}->execute($nick_name,
                                         $path_data_dir,
                                         $project_name) or die $dbh->errstr;
}

################################################################################
#
# create_directories_for
#
# For a given path to a file, it checks the directory structure specified
# exists. If it does not is creates it.
#
# It can cope with these;
#         filename_on_its_own
#         ./filename
#         /tmp/filename <--- SPECAIL CASE
#         annotation/Trypanosoma/brucei/chr02/Tb927_02_mV2.embl
##

sub create_directories_for {
    my ($self, $new_file) = @_;
    
    my @dirs = split('/', $new_file);
    pop @dirs;
    print "dirs = ", join(", ", @dirs), "\n" if $self->{'debug'};


    ##
    # SPECIAL CASE: Fix path that has slash as first char i.e. /tmp/thing
    ##
    if ($new_file =~ m|^/|) {
        shift @dirs;
        $dirs[0] = "/".$dirs[0];
    }

    ##
    # Create each dir as need be
    ##
    my $prev_dirs = "";

    foreach my $dir (@dirs) {
        next if ($dir eq '.');

        my $whole_dir = "${prev_dirs}${dir}";

        unless (-d $whole_dir) {
            print "Making directory: $whole_dir\n" if $self->{'debug'};
            mkdir "$whole_dir" or die "ERROR: Unable to create directory $whole_dir";
        }
        else {
            print "$whole_dir already exists\n" if $self->{'debug'};
        }

        $prev_dirs .= "${dir}/";
    }

}



################################################################################
#
# Auto disconnect from the DB when this object is destroyed
#
################################################################################

sub DESTROY {
    my $self = shift;
    die "Wrong object type\n" unless ref($self);
    $self->disconnect();
    return;
}


1;
