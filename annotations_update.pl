#!/usr/local/bin/perl -w

#
# $Id: annotations_update.pl,v 1.48 2008-03-27 09:47:43 adh Exp $
#

############################################
############################################
#
# MAKE SURE EMBL FILES HAVE A SEQUENCE IDS
#
############################################
############################################

# Constraint variables declaration before using them
use strict;
use DaemonsLog::Log;
##########################################

# Bio modules
# PSU
use Bio::PSU::IO::Blast::Result;
use Bio::PSU::IO::Blast::Search;
use Bio::PSU::SeqFactory;
use Bio::PSU::Range;
use Bio::PSU::Feature;
use Bio::PSU::Utils::FeatureIdHelper::BioPSUHelper;

use File::Temp qw/tempdir/;
use Getopt::Long;
##########################################

# Benchmark Module
use Benchmark;
##################################################################

use Data::Dumper;
##################################################################

##########################################
##
# Initialisation
##
##########################################

my $t1 = Benchmark->new ();

my $debug = 0;

my $_warning_colour = 10;
my $_error_colour   = 2;

my @mapped_qualifiers_to_keep = ('gene', 'systematic_id', 'curation', 'psu_db_xref', 'label', 'synonym');
my @newseq_qualifiers_to_keep = ('all');

#########################################
##
# Use of NCBI Blast 2.2.6, with customized parameters, seems much better than WuBlast !!!
##
#########################################

my $blast_software = "wublast";
$blast_software    = "ncbiblast";

# Blast database indexing command
my $_indexingdb = undef;
# Blast search command
my $_blastn = undef;

if ($blast_software eq "wublast") {
  
  $_indexingdb = "/usr/local/pubseq/bin/wu_blast/pressdb";
  $_indexingdb = "pressdb";
  
  # WU blast2n last version - good and faster than blastn 1.4
  $_blastn   = "/usr/local/pubseq/bin/wu_blast/blastn";
}
else {
  $_indexingdb = "formatdb -p f -i ";
  $_blastn = "blastall -p blastn ";
}

my $oldseqfile    = '';
my $oldseqformat  = '';
my $newseqfile    = '';
my $newseqformat  = '';
my $tabfile       = '';
my $taboutput     = 'mappedAnnotations.embl';
my $speedup       = 0;
my $filtering     = 90;
my $merging       = 0;
my $extension     = 500;
my $conflictScore = 300;
my $mapped_qualifiers_to_keep_file = $ENV{HOME} . '/.mapped_qualifiers_to_keep';
my $newseq_qualifiers_to_keep_file = $ENV{HOME} . '/.newseq_qualifiers_to_keep';

if (!(&GetOptions('oldseq=s' => \$oldseqfile, "newseq=s" => \$newseqfile, "oldseq_if=s" => \$oldseqformat, "newseq_if=s" => \$newseqformat, "tab=s" => \$tabfile, "output=s" => \$taboutput, "speedup" => \$speedup, "filtering=s" => \$filtering, "merging" => \$merging, "extension=s" => \$extension, "oq=s" => \$mapped_qualifiers_to_keep_file, "nq=s" => \$newseq_qualifiers_to_keep_file))) {
  print STDERR "Error: problem processing command-line options\n";
  display_help ();
  die("Error: problem processing command-line options\n");
}

if ($extension != 0) {
  print STDERR "\nConflict definition: there is a conflict when two hits have a similar score - not much than $conflictScore difference\n";
}

if ($extension != 0) {
  print STDERR "In case of conlicts (ie several potential mapping locations), the features to map will be extended of $extension bp upstream and downstream for better mapping!!\n\n";
}

if ($oldseqfile eq '') {
  print STDERR "\noldseq parameter missing!\n\n";
  display_help ();
  exit 1;
}
elsif ($newseqfile eq '') {
  print STDERR "\nnewseq parameter missing!\n\n";
  display_help ();
  exit 1;
}
if ($oldseqformat eq '') {
  print STDERR "\noldseq_if parameter missing!\n\n";
  display_help ();
  exit 1;
}
elsif ($newseqformat eq '') {
  print STDERR "\nnewseq_if parameter missing!\n\n";
  display_help ();
  exit 1;
}
elsif ($newseqformat eq 'fasta') {
  # newseq is a fasta sequence then we don't need any merging process
  $merging = 0;
}

my $regex;
my $id_helper =
  Bio::PSU::Utils::FeatureIdHelper::BioPSUHelper->new(-selection_regex => $regex, -verbose => $debug);

if ($speedup) {
  print STDERR "\nSpeedup mode activated...\n";
  print STDERR "bear in mind that the results may not be consistent or the mapping may not be fully done!\n\n";
}

if ($merging) {
  print STDERR "\nWARNING: If the EMBL file, $newseqfile, the features are mapped to, contains annotations, they will be merged into a unique annotation file, $taboutput\n\n";
  
  # process the files
  
  if (not -f $mapped_qualifiers_to_keep_file) {
    print STDERR "No file for mapped qualifiers to keep specified, default list is: " . join (',', @mapped_qualifiers_to_keep) . "\n";
  }
  else {
    open OLDFILE, "<$mapped_qualifiers_to_keep_file" or die "can't open file, $mapped_qualifiers_to_keep_file!\n";
    my @list_tmp = ();
    while (<OLDFILE>) {
      chomp ($_);
      push (@list_tmp, $_);
    }
    close OLDFILE;
    @mapped_qualifiers_to_keep = @list_tmp;
  }
  
  if (not -f $newseq_qualifiers_to_keep_file) {
    print STDERR "No file for new qualifiers to keep specified, default list is: " . join (',', @newseq_qualifiers_to_keep) . "\n";
  }
  else {
    open NEWFILE, "<$newseq_qualifiers_to_keep_file" or die "can't open file, $newseq_qualifiers_to_keep_file!\n";
    my @list_tmp = ();
    while (<NEWFILE>) {
      chomp ($_);
      push (@list_tmp, $_);
    }
    close NEWFILE;
    @newseq_qualifiers_to_keep = @list_tmp;
  }
  
}
else {
  print STDERR "no merging process.\n";
}

print STDERR "Loading the sequence and feature objects from files...\n";

###############################
# Get the old sequence object
# the annotations are within this object or within a tab file
##

my $old_s1 = Bio::PSU::SeqFactory->make(
                                        -file   => "<$oldseqfile",
                                        -format => $oldseqformat
                                       );

my $oldseq_obj = $old_s1->next_seq;

if (not defined $oldseq_obj) {
  print STDERR "can\'t get the old sequence object from file, $oldseqfile!\n";
  print STDERR "Are you sure you specified the right format, $oldseqformat?\n";
  
  exit 1;
}

# check if there are annotations to map !!

if ((! $oldseq_obj->has_features) and ($tabfile eq '')) {
  print STDERR "\nNo features to map !!\n\n";
  display_help ();
  exit 1;
}

# my $blast_directory = "blast_search";
#my $blast_directory = tempdir( 
#			      'BLASTtmpdirXXXXXX',
#			      DIR => File::Spec->curdir,
#			      CLEANUP => 1 
#			     );

my $blast_directory = tempdir(
                              'BLASTtmpdirXXXXXX',
                               DIR => "/tmp/",
                               CLEANUP => 1
                             );

if ($debug) {
  print STDERR "blast directory name: $blast_directory\n";
}

###########################
# Get the new sequence objects
##

my ($subjectDb_file, %new_seq_objects) = getNewSequenceObjects ($newseqfile, $newseqformat, $blast_directory);

my @report_features = ();
my @failed_features = ();

# the tab object

my $tab_seq_obj;
if (($tabfile ne '')) {
  my $tab_obj = Bio::PSU::SeqFactory->make(-file   => "<$tabfile",
                                           -format => 'embl'
                                          );
  $tab_seq_obj = $tab_obj->next_seq;  
}
else {
  $tab_seq_obj = $oldseq_obj;
}

#####################################################
# Mapping ...
#####################################################

if ((my $nb_features = $tab_seq_obj->has_features) > 0) {
  print STDERR "$nb_features features are being mapped ...\n";
  my $i = 1;
  my $percent_memorized = -1;
  
  foreach my $feature ($tab_seq_obj->features) {
    
    if (($feature->key eq "source") || ($feature->key eq "gap")) {
      # don't map this one
      next;
    }
    
    if ($feature->key eq "CDS") {
      my $query_gene_name = getGeneName ($feature);
      print STDERR "\nmapping new feature, $query_gene_name...\n";
    }
    else {
      print STDERR "\nmapping new " . $feature->key . " feature...\n";
    }
    
    my $merged = 0;
    my $merged_gene_name = undef;
    
    my ($mapped_feature, $best_hsps_ref) = mapping ($feature, $oldseq_obj, \%new_seq_objects, $speedup, $filtering, $subjectDb_file, $blast_directory, $extension, $conflictScore);
    
    # my $seq_id = $best_hsp_hash->{seq_id};
    
    if (defined ($mapped_feature)) {
      
      if ($debug) {
        print STDERR "Mapping done,\n";
      }
      my $best_hsp_hash  = $best_hsps_ref->{1};
      my $seq_id         = $best_hsp_hash->{seq_id};
      
      my $partial_mapping = 0;
      my @oldranges  = $feature->ranges;
      my @newranges  = $mapped_feature->ranges;
      if (@newranges < @oldranges) {
	print STDERR "partial mapping done\n";
	$partial_mapping = 1;
      }
      
      my $seq_obj = $new_seq_objects{$seq_id};
      
      if ($merging && !$partial_mapping) {
	my @newseq_features = $seq_obj->features;
        my ($overlapping_features) = getOverlappingFeatures ($mapped_feature, @newseq_features);
        if (defined $overlapping_features) {
          my @indexes = keys (%$overlapping_features);
          if (@indexes > 1) {
            print STDERR "Found more than one overlapping feature - don\'t process any merging!\n";
            $seq_obj->features ($mapped_feature);
          }
          else {
	    print STDERR "Found one overlapping feature\n";
	    $merged = 1;
            my $index = $indexes[0];
            my $overlaped_feature = $overlapping_features->{$index};
	    # get the overlap feature gene name if CDS
	    if ($overlaped_feature->key eq "CDS") {
	      $merged_gene_name = getGeneName ($overlaped_feature);
	    }
	    
            my $merged_feature = mergeFeatures ($mapped_feature, $overlaped_feature, \@mapped_qualifiers_to_keep, \@newseq_qualifiers_to_keep);
            
	    # Replace the overlapping one bt the merged one
	    splice (@newseq_features, $index, 1, $merged_feature);
	    
            print STDERR "Merging done!\n";
	    
	    my $new_obj = $seq_obj->clone (-features => \@newseq_features);
	    $seq_obj = $new_obj;
	    
	    # this way the mapped one is attached to the new sequence and it's possible to get the sequence string - useful for the report - see further down
	    $mapped_feature = $merged_feature;
	    
          }
	}
	else {
	  print STDERR "no overlapping feature found!\n";
	  $seq_obj->features ($mapped_feature);
	}
      }
      else {
        $seq_obj->features ($mapped_feature);
      }
      
      $new_seq_objects{$seq_id} = $seq_obj;
      
      # Feature mapping Report
      
      my $report_feature = Bio::PSU::Feature->new (
						   -key => $feature->key,
						   -ranges => [$feature->ranges]
                                                  );
      
      # Get the gene name(s)
      
      my $id = getGeneName ($feature);
      if (defined $id) {
	$report_feature->qadd ('systematic_id', $id);
      }
      
      my $note = "mapped to sequence, " . $best_hsp_hash->{seq_id} . ", on strand, " . $mapped_feature->strand . ", at position ";
      foreach my $range ($mapped_feature->ranges) {
        $note .= $range->start . ".." . $range->end . ", ";
      }
      $report_feature->qadd ('note', $note);
      $note = "mapped to sequence, " . $best_hsp_hash->{seq_id} . ", with confidence, percentage_id = " . $best_hsp_hash->{percent_id} . ", score = " . $best_hsp_hash->{score} . ", evalue = " . $best_hsp_hash->{evalue};
      $report_feature->qadd ('note', $note);
      
      my $tomapped_f_length = 0;
      foreach my $range ($feature->ranges) {
        $tomapped_f_length += 1 + $range->end - $range->start;
      }
      
      if ((length ($mapped_feature->str)) != $tomapped_f_length) {
        $note = "WARNING - mapped length is " . length ($mapped_feature->str) . ", expecting $tomapped_f_length";
        $report_feature->qadd ('exception', $note);
	$report_feature->qadd ('colour', $_warning_colour);
      }
      
      if ($merged) {
	if (defined $merged_gene_name) {
	  $note = "MERGED with one of the features, whose gene name is, $merged_gene_name, attached to the new sequence";
	}
	else {
	  $note = "MERGED with one of the features attached to the new sequence";
	}
	$report_feature->qadd ('note', $note);
      }
      
      if ($partial_mapping) {
	$note = "WARNING - Partial mapping done - the mapped feature represents the region where the mapped feature is expected";
	$report_feature->qadd ('note', $note);
	$report_feature->qadd ('colour', $_warning_colour);
      }
      
      push (@report_features, $report_feature);
      
    }
    else {
      
      my $report_feature = Bio::PSU::Feature->new (
						   -key => $feature->key,
						   -ranges => [$feature->ranges]
                                                  );
      
      # Get the gene name
      
      my $id = getGeneName ($feature);
      
      if (defined $id) {
        $report_feature->qadd ('systematic_id', $id);
      }
      
      # Leave a note about why it failed (ie no significant hits or more several ones!!)
      
      if (not defined $best_hsps_ref) {
	my $note = "ERROR - Mapping failed - no significant hit found";
	$report_feature->qadd ('note', $note);
	$report_feature->qadd ('colour', $_error_colour);
      }
      else {
	my $note = "ERROR - Mapping failed - several significant hits found";
	$report_feature->qadd ('note', $note);
	foreach my $index (keys (%$best_hsps_ref)) {
	  my $best_hsp_hash = $best_hsps_ref->{$index};
	 
	  #$note = "potential mapping to sequence, " . $best_hsp_hash->{seq_id} . ", with confidence, percentage_id = " . $best_hsp_hash->{percent_id} . ", score = " . $best_hsp_hash->{score} . ", evalue = " . $best_hsp_hash->{evalue};
	  #$report_feature->qadd ('note', $note);

	  $note = "potential mapping to sequence, " . $best_hsp_hash->{seq_id} . ", on strand, " . $best_hsp_hash->{strand} . ", at position " . $best_hsp_hash->{start} . ".." . $best_hsp_hash->{end};
	  $report_feature->qadd ('note', $note);
	  $note = "mapping with confidence, percentage_id = " . $best_hsp_hash->{percent_id} . ", score = " . $best_hsp_hash->{score} . ", evalue = " . $best_hsp_hash->{evalue};
	  $report_feature->qadd ('note', $note);
	  $report_feature->qadd ('colour', $_error_colour);

	}	  

      }
      
      push (@report_features, $report_feature);
    }
    
    my $percent = int ($i*100/$nb_features);
    if ((($percent%5) == 0) && ($percent != $percent_memorized)) {
      print STDERR "\n$percent% done.\n\n";
      $percent_memorized = $percent;
    }
    $i++;
  }
}

#####################################################
# Mapping done
# copy the new set of features into a new tab file
# copy the failed features into the failure tab file
#####################################################

#############
# Mapped Features
#############

print STDERR "copying the mapped features into $taboutput file...\n";

my $new_tab_obj = Bio::PSU::SeqFactory->make(-file   => ">$taboutput",
					     -format => 'embl');
foreach my $seq_id (keys (%new_seq_objects)) {
  $new_tab_obj->write_seq($new_seq_objects{$seq_id});
}

#############
# Report Features
#############

my $report_file = $oldseq_obj->id . ".report.tab";

print STDERR "report file: $report_file\n";

print STDERR "copying the report features into $report_file...\n";

my $report_seq_obj = Bio::PSU::Seq->new (-id   => $tab_seq_obj->id,
			 	         -desc => "Reporting mapping of features from sequence, " . $tab_seq_obj->id,
				         -type => $tab_seq_obj->type);
$report_seq_obj->features (@report_features);
my $report_tab_obj = Bio::PSU::SeqFactory->make(-file   => ">$report_file",
						-format => 'embl');
$report_tab_obj->write_seq($report_seq_obj);


#############
# failed features - NO MORE !!
#############

#my $failed_tab_seq_obj = Bio::PSU::Seq->new (-id   => $tab_seq_obj->id,
#					     -desc => $tab_seq_obj->desc,
#					     -type => $tab_seq_obj->type);
#$failed_tab_seq_obj->features (@failed_features);

# copy this seq obj into a new file

#my $failed_tab_obj = Bio::PSU::SeqFactory->make(-file   => ">$failedfeatfile",
#						-format => 'embl');
#$failed_tab_obj->write_seq($failed_tab_seq_obj);

my $t2 = Benchmark->new ();
print  STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";

###################################################################
# The end
###################################################################


###################################################################
# Methods
###################################################################

# display the help message

sub display_help {
  print "annotations_update.pl --oldseq [old sequence file] --oldseq_if [old sequence file format (Fasta or EMBL)] --newseq [new sequence file] --newseq_if [new sequence file format (Fasta or EMBL)] --tab [annotation file] --output [embl file with the mapped annotations and the sequence] --filtering [blast percentage identity] --speedup --merging --extension [number of bp the features to map will be extended of] --oq [old sequence qualifiers to map over] --nq [new sequence qualifiers to keep]\n\n";
  print "\t--tab the annotation file (a tab file - with or without the sequence) - COMPULSORY\n";
  print "\t--oldseq the previous sequence file (in fasta or EMBL format)         - COMPULSORY\n";
  print "\t--oldseq_if the previous sequence file format (fasta or EMBL format)  - COMPULSORY\n";
  print "\t--newseq the new sequence file (in embl or fasta format)              - COMPULSORY\n";
  print "\t--newseq_if the new sequence file format  (embl or fasta format)      - COMPULSORY\n";
  print "\t--output the output file name                                         - DEFAULT is mappedAnnotations.embl\n";
  print "\t--speedup blast search speedup by softmasking the query - speedup is not activated per default\n";
  print "\t--filtering give the blast percentage identity limit below the feature won't be mapped  - DEFAULT is 90%\n";
  print "\t--merging merge the mapped features with overlapping features already attached to the new sequence\n";
  print "\t--extension The number of bp the features are extended of (upstream and downstream) for a more accurate mapping - DEFAULT is 2000 bp\n";
  print "\t--oq Specify the file location giving the set of old sequence qualifiers to map over    - DEFAULT is HOME/.mapped_qualifiers_to_keep\n";
  print "\t--nq Specify the file location giving the set of new sequence qualifiers to keep        - DEFAULT is HOME/.newseq_qualifiers_to_keep\n";
}


###################################################################
# get min and max from a set of ranges
# the purpose is to blast against a big range to make sure the mapping will be correct
# so if the gene contains several exons, 
# the query sequence will be the all gene sequence - exons + introns.
###################################################################

sub create_range {
  my (@ranges) = @_;
  
  my $strand;
  my $start = 10**12;
  my $stop  = 0;
  
  foreach my $range (@ranges) {
    if ($range->start < $start) {
      $start = $range->start;
    }
    if ($range->end > $stop) {
      $stop = $range->end;
    }
    $strand = $range->strand;
  }
  my $new_range = Bio::PSU::Range->new(-start       => $start,
				       -end         => $stop,
				       -strand      => $strand);
  return $new_range;
}


sub setBlastOptions {
  my ($speedup, $is_exon_mapping, $blast_directory) = @_;
  my $blastOptions = {
		      'speedup'         => $speedup,
		      'exon_mapping'    => $is_exon_mapping,
		      'blast_directory' => $blast_directory,
		     };
  
  return $blastOptions;
}

####################################################################
# run a blast search and return the results as a result PSU object
# use ONLY WuBLAST 1.4 !!
##################################################################

sub doBlast {
  my ($queryFile, $subjectFile, $blastOptions, $reprocessing) = @_;
  
  my $blast_dir       = $blastOptions->{blast_directory};
  my $speedup         = $blastOptions->{speedup};
  my $is_exon_mapping = $blastOptions->{exon_mapping};
  
  my $blast_options = "";
  
  if ($_blastn =~ /blastall/) {
    
    # Filtering deactivated => "-F F"
    
    # two modes

    # less stringent
    $blast_options .= "-a 3 -K 10 -q -1 ";

    # filtering and more stringent - not compatible with extension option !!
    # $blast_options .= "-a 3 -K 10 ";

    if ($is_exon_mapping) {
      # no filtering, short words
      $blast_options = $blast_options . "-W 7 -F F";
    }
    elsif ($speedup || not $reprocessing) {
      # softmasking and big words
      $blast_options = $blast_options . "-W 20 -F \"m D\"";
    }
    else {
      # no speedup and reprocessing => more sensitive

      # no masking and big words
      # more hits so more sensitive but much slower
      $blast_options = $blast_options . "-W 20 -F F";
    }
  }
  else {
    $blast_options .= "-warnings -hspmax 20 -cpus 3";
    if ($is_exon_mapping) {
      # small words
      $blast_options = $blast_options . "W=7";
    }
    elsif ($speedup) {
      # big words and softmasking with both seq and dust
      $blast_options = $blast_options . "W=20 -span wordmask=seg wordmask=dust";
    }
    else {
      # big words, no masking
      $blast_options = $blast_options . "W=20";
    }
  }
  
  if ($debug) {
    print STDERR "\tQuery File: $queryFile\n";
    print STDERR "\tSubject File: $subjectFile\n";
    print STDERR "\tblast_options: $blast_options\n";
  }
  
  print STDERR "blast search running ...\n";
  
  # generate a buffer - BUG !!!!!!!!!!!!!
  # my @blast_output = qx/$_blastn $subjectFile $queryFile $blast_options/;
  
  # generate a file
  my $blast_output_name = "$blast_directory/blast.res";

  my @args = ();
  if ($_blastn =~ /blastall/) {
    push (@args, "$_blastn -d $subjectFile -i $queryFile $blast_options > $blast_output_name");
  }
  else {
    push (@args, "$_blastn $subjectFile $queryFile $blast_options > $blast_output_name");
  }
  system (@args);
  
  if ($debug) {
    print STDERR "search done.\n";
    # print STDERR "Dumping blast output array: " . Dumper (@blast_output) . "\n";
  }
  
  #my $bfh    = Bio::PSU::IO::BufferFH->new(-buffer => \@blast_output);
  #my $search = Bio::PSU::IO::Blast::Search->new(-bfh => $bfh);
  my $bfh_2  = Bio::PSU::IO::BufferFH->new(-file => $blast_output_name);
  my $search = Bio::PSU::IO::Blast::Search->new(-bfh => $bfh_2);
  
  my $result = $search->next_result;
  
  if ($debug && not defined $result) {
    print STDERR "no search result!\n";
  }
  else {
    # print STDERR "returning a result object reference...\n";
  }
  
  return $result;
}

sub processResults {
  my ($filtering, $result, $extension, $conflictScore) = @_;
  
  my $has_been_reversed = undef;
  my $mapped_range  = undef;
  
  my $best_hsps_ref = getBestHSP ($filtering, $result, $conflictScore);
  if (defined ($best_hsps_ref)) {
    if (keys (%$best_hsps_ref) == 1) {
      my $best_hsp_hash = $best_hsps_ref->{1};
      my $hsp = $best_hsp_hash->{hsp};
      ($mapped_range, $has_been_reversed) = getRangeBestHSP ($hsp, $extension);
    }
    else {
      # conflict
      # do not do any mapping, report the potential mapping positions though
      if ($debug) {
        print STDERR "conflict\n";
      }
    }
  }
  
  return ($best_hsps_ref, $mapped_range, $has_been_reversed);
}  

####################################################################
# select the best HSP and return it
####################################################################

# HSP : High Scoring Pair

sub getBestHSP {
  my ($hit_filtering, $result, $conflictScore) = @_;
  
  my %best_hsps;
  my $i = 1;
  my $score = 0;
  my $percent_id = 0;

  # $conflictScore
  # a conflict is when the query hits two HSPs wih the same score
  # now there is a conflict when their percentages of identity are close one to another.
  
  while (defined (my $hit = $result->next_hit)) {
    while (defined (my $hsp = $hit->next_hsp)) {
      if (($hsp->score-$conflictScore) > $score) {
      # if ($hsp->score > $score) {

	$score      = $hsp->score;
	$percent_id = $hsp->percent;

	# (re)initializing the best hsp hash

	print STDERR "new hsp - reinitialisation\n";
	print STDERR "\tscore: "   . $hsp->score . "\n";
	print STDERR "\tpercent: " . $hsp->percent . "\n";

	my %best_hsps_tmp;
	$i = 1;
	
	my $best_hsp_hash = {
			     'seq_id' => $hit->s_id,
			     'percent_id' => $hsp->percent,
			     'score'  => $hsp->score,
			     'evalue' => $hsp->expect,
			     'start'  => $hsp->s_begin,
			     'end'    => $hsp->s_end,
			     'strand' => $hsp->s_strand,
			     'hsp'    => $hsp,
			    };
	$best_hsps_tmp{$i} = $best_hsp_hash;
	%best_hsps = %best_hsps_tmp;
	$i++;
      }
      elsif ($hsp->score >= ($score-$conflictScore)) {

	# we consider this case as a conflict (because two hits with the same score or two hits with very close score) - in that case it's difficult to decide which one to choose.

	# should we change the score and percent_id we take for reference ???

	print STDERR "add an extra hsp\n";
	print STDERR "\tscore: "   . $hsp->score . "\n";
	print STDERR "\tpercent: " . $hsp->percent . "\n";

	my $best_hsp_hash = {
			     'seq_id' => $hit->s_id,
			     'percent_id' => $hsp->percent,
			     'score'  => $hsp->score,
			     'evalue' => $hsp->expect,
			     'start'  => $hsp->s_begin,
			     'end'    => $hsp->s_end,
			     'strand' => $hsp->s_strand,
			     'hsp'    => $hsp,
			    };
	
	$best_hsps{$i} = $best_hsp_hash;
	$i++;
      }
    }
  }
  
  # processing done!
  if (keys (%best_hsps) == 0) {
    print STDERR "\tno hits found\n";
    return undef;
  }
  elsif (keys (%best_hsps) == 1) {
    
    my $best_hsp_hash = $best_hsps{1};
    
    print STDERR "\thit found whose seq id is: " . $best_hsp_hash->{seq_id} . "\n";
    print STDERR "\thit found whose score is : " . $best_hsp_hash->{score} . "\n";
    print STDERR "\thit found whose percentage identity is : " . $best_hsp_hash->{percent_id} . "\n";
    
    if ($best_hsp_hash->{percent_id} >= $hit_filtering) {
      return \%best_hsps;
    }
    else {
      print STDERR "\tno significant hit found\n";
      return undef;
    }
  }
  else {
    # conflict - more than one best hit !!
    return \%best_hsps;
  }
}


####################################################################
# create a new range object from a HSP object
####################################################################

sub getRangeBestHSP {
  my ($best_hsp, $extension) = @_;
  
  my $start  = $best_hsp->s_begin;
  my $end    = $best_hsp->s_end;
  my $strand = $best_hsp->q_strand;
  
  if ($debug) {
    # print STDERR "dumping best_hsp: " . Dumper ($best_hsp) . "\n";
  }
  
  my $has_been_reversed = 0;
  
  if ($start > $end) {
    
    # NCBI Blast output reports reversed coordinates if the query hits the reversed strand of the subject sequence.
    
    if ($debug) {
      print STDERR "NCBI Blast search, start and end corrdinates to reverse !! strand: $strand\n";
    }
    
    my $tmp = $start;
    $start  = $end;
    $end    = $tmp;
    
    $has_been_reversed = 1;
    
  }
  
  # Extension recorrection
  # Bear in mind that this recorrection may not be right !!
  $start += $extension;
  $end   -= $extension;

  my $range = Bio::PSU::Range->new(-start       => $start,
				   -end         => $end,
				   -strand      => $strand);
  # always on the forward strand
  # if not, means the new sequence is the reverse and complement of the old sequence
  # return this information
  
  if ($strand == -1) {
    $has_been_reversed = 1;
  }
  
  return ($range, $has_been_reversed);
}


sub setSubRangeMappedToAllSequence {
  my ($sub_range, $sub_mapped_range, $strand) = @_;
  
  my $range_mapped_to_all_sequence = Bio::PSU::Range->new (
							   -start => $sub_range->start + $sub_mapped_range->start - 1,
							   -end => $sub_range->start + $sub_mapped_range->end - 1,
							   -strand => $strand,
							  );
  
  return $range_mapped_to_all_sequence;
}


sub createFeature {
  my ($feature, @new_ranges) = @_;
  my $new_feature = Bio::PSU::Feature->new (
					    -key    => $feature->key,
					    -ranges => \@new_ranges,
					   );
  foreach my $qualifier ($feature->qnames) {
    $new_feature->qadd ($qualifier, $feature->qvalues ($qualifier));
  }
  
  return $new_feature;
}


sub check_gene_accurracy {
  my ($feature, $newseq_obj) = @_;
  my @notes = ();
  
  # would get confused if cds with exons on different strands
  
  my $i = 1;
  my $nb_ranges = $feature->ranges+0;
  if ($feature->has_ranges) {
    foreach my $range ($feature->ranges) {
      
      my $sequence = $newseq_obj->str;
      
      # check on the appropriate strand
      
      # feature on the forward strand
      
      if ($range->strand == 1) {
	# check if the first exon starts with a start codon
	if ($debug) {
	  print STDERR "checking on the forward strand...\n";
	}
	if ($i == 1) {
	  my $first_codon = substr ($sequence, $range->start-1, 3);
	  if ($debug) {
	    print STDERR "first codon : $first_codon - (should be atg)\n";
	  }
	  if (not (lc($first_codon) =~ /atg/)) {
	    my $note = "first exon doesn't start with a start codon";
	    push (@notes, $note);
	  }
	}
	# check if the last exon ends with a stop codon
	if ($i == $nb_ranges) {
	  my $stop_codon = substr ($sequence, $range->end-3, 3);
	  if ($debug) {
	    print STDERR "stop codon : $stop_codon - (should be taa or tag or tga)\n";
	  }
	  if (not (lc($stop_codon) =~ /taa|tag|tga/)) {
	    my $note = "last exon doesn't end with a stop codon";
	    push (@notes, $note);
	  }
	}
      }
      # feature on the reverse strand
      else {
	# check if the first exon starts with a start codon
	if ($debug) {
	  print STDERR "checking on the reverse strand...\n";
	}
	if ($i == $nb_ranges) {
	  my $first_codon = substr ($sequence, $range->end-3, 3);
	  if ($debug) {
	    print STDERR "first codon : $first_codon - (should display cat)\n";
	  }
	  if (not (lc($first_codon) =~ /cat/)) {
	    my $note = "first exon doesn't start with a start codon";
	    push (@notes, $note);
	  }
	}
	# check if the last exon ends with a stop codon
	if ($i == 1) {
	  my $stop_codon = substr ($sequence, $range->start-1, 3);
	  if ($debug) {
	    print STDERR "stop codon : $stop_codon - (should display tta or cta or tca)\n";
	  }
	  if (not (lc ($stop_codon) =~ /tta|cta|tca/)) {
	    my $note = "last exon doesn't end with a stop codon";
	    push (@notes, $note);
	  }
	}
      }
      $i++;
    }
  }  
  
  if ((@notes+0) > 0) {
    $feature->qadd('note', @notes);
  }
  
  return $feature;
}

# take as input a feature and map it into the new sequence
# return a feature with a new range object
# if feature without range, return the given feature without any processing
# if failed, return undef

sub mapping {
  my ($feature, $oldseq_obj, $new_seq_objects_ref, $speedup, $filtering, $subjectDb_file, $blast_directory, $extension, $conflictScore) = @_;
  my %new_seq_objects = %$new_seq_objects_ref;
  my $is_exon_mapping = 0;
  
  if ($feature->has_ranges) {
    my @ranges = $feature->ranges;
    my $nb_ranges = @ranges;
    
    # mapping of the feature
    
    # create a range which will be converted into the query sequence, by taking the min and the max corrdinates
    # min will the start position
    # max will be the end position
    my $range  = create_range (@ranges);
    
    # Set the QueryFile - no need to set the subject file
    my $querySequenceObject = getSubSequenceObject ($range, $oldseq_obj, 0);
    
    if (length ($querySequenceObject->str) < 20) {
      print STDERR "the query sequence is shorter than the word length, W=20.";
      print STDERR "Can not process the mapping\n";
      return (undef, undef);
    }
    
    my $queryFile       = setQueryFile ($querySequenceObject, $blast_directory);
    my $blastOptions    = setBlastOptions ($speedup, $is_exon_mapping, $blast_directory);
    
    my $result = doBlast ($queryFile, $subjectDb_file, $blastOptions, 0);
    
    my ($best_hsps_ref, $mapped_range, $has_been_reversed) = processResults ($filtering, $result, 0, $conflictScore);
    
    # my $seq_id = $best_hsp_hash->{seq_id};
    
    if ($debug && defined $has_been_reversed) {
      print STDERR "has_been_reversed? $has_been_reversed\n";
    }
    
    if (not defined $mapped_range) {
      
      # two failure modes
      # 1/ no significant hits
      # 2/ several significant hits
      
       # speedup action mode here !!
       # if speedup, no reprocessing

      if ($debug) {
	print STDERR "failed - mapped_range not defined\n";
      }
      if ($speedup && not defined $best_hsps_ref) {
	print STDERR "because no significant hits\n";
	return (undef, undef);
      }
      else {
	print STDERR "no hit or several significant hits - reprocessing the mapping of a larger region containing the feature...\n";

	# Redo the blast with the feature extension
	
	$querySequenceObject = getSubSequenceObject ($range, $oldseq_obj, $extension);

	$queryFile = setQueryFile ($querySequenceObject, $blast_directory);
	$result = doBlast ($queryFile, $subjectDb_file, $blastOptions, 1);
    
	my ($extended_best_hsps_ref, $extended_mapped_range, $extended_has_been_reversed) = processResults ($filtering, $result, $extension, $conflictScore);

	if (keys (%$extended_best_hsps_ref) > 1) {

	  print STDERR "still a conflit :-(\n";
	  
	  # still conflict - return the previous run
	  return (undef, $best_hsps_ref);
	}
	elsif (keys (%$extended_best_hsps_ref) == 1) {

	  print STDERR "conflict resolved :-)\n";

	  # conflict resolved !!
	  # return ($extended_mapped_range, $extended_best_hsps_ref);
	  $best_hsps_ref = $extended_best_hsps_ref;
	  $mapped_range = $extended_mapped_range;
	  $has_been_reversed = $extended_has_been_reversed;
	  
	}
	else {
	  # no significant hits or hit below filtering %id
	  return (undef, undef);
	}
      }
    }

    if (not defined $mapped_range) {
      print STDERR "mapping failed - also the extension - should not be here - don't know why we get here !!!!!\n";
      return (undef, undef);
    }

    if ($nb_ranges == 1) {
      
      my $best_hsp_hash = $best_hsps_ref->{1};
      my $seq_id = $best_hsp_hash->{seq_id};
      
      # if only one range, the mapping for this feature is done
      if ($debug) {
	print STDERR "only one range - no submapping\n";
      }
      if (! $has_been_reversed) {
	$mapped_range->strand($feature->strand);
      }
      else {
	$mapped_range->strand(- $feature->strand);
      }
      my @new_ranges = ($mapped_range);
      
      # create mapped_feature
      
      my $mapped_feature = createFeature ($feature, @new_ranges);
      
      # check if it makes sense
      # ie if fisrt exon has a start codon
      # and the last exon the stop codon
      # if not accurate, add a note message
      
      # deactivated !!
      # should be reported in the report file

      #my $newseq_obj = $new_seq_objects{$seq_id};
      #if ($mapped_feature->key =~ /cds/i) {
	# $mapped_feature = check_gene_accurracy ($mapped_feature, $newseq_obj);
      #}
      
      return ($mapped_feature, $best_hsps_ref);
    }
    else {
      print STDERR "submapping...\n";
      if ($debug) {
	print STDERR "working on the following range (on the new sequence) for this feature: start, end " . $mapped_range->start . ", " . $mapped_range->end . "\n"; 
      }
      
      my $best_hsp_hash = $best_hsps_ref->{1};
      my $seq_id = $best_hsp_hash->{seq_id};
      
      $is_exon_mapping = 1;
      my $blastOptions = setBlastOptions ($speedup, $is_exon_mapping, $blast_directory);
      
      my $newseq_obj = $new_seq_objects{$seq_id};
      my $subjectSequenceObject = getSubSequenceObject ($mapped_range, $newseq_obj, 0);
      my $subjectFile = setSubjectFile ($subjectSequenceObject, $blast_directory);      
      
      # set the strand - each sub range have a unique one
      # presupposed that they are all in the same strand !!!!
      
      my $strand = 0;
      my @ranges_tmp = $feature->ranges;
      my $range1 = $ranges_tmp[0];
      if (! $has_been_reversed) {
	$strand = $range1->strand;
      }
      else {
	$strand = - $range1->strand;
      }
      
      # foreach subfeature for the current feature (e.g.) a set of exons, map into the new seq.
      # if one of the subfeature can not be mapped, copy the current feature into a failure tab file.
      
      my @mapped_ranges = ();
      my $index = 0;
      foreach my $sub_range ($feature->ranges) {
	if ($debug) {
	  print STDERR "mapping new sub feature mapping...\n";
	}
	
	my $querySequenceObject = getSubSequenceObject ($sub_range, $oldseq_obj, 0);
	
	if (length ($querySequenceObject->str) < 7) {
	  print STDERR "the query sequence is shorter than the word length, W=7.\n";
	  print STDERR "Can not process the sub mapping\n";
	  
	  # Report the mapped region though
	  
	  if (! $has_been_reversed) {
	    $mapped_range->strand($ranges[0]->strand);
	  }
	  else {
	    $mapped_range->strand(- $ranges[0]->strand);
	  }
	  my @new_ranges = ($mapped_range);
	  
	  # create mapped_feature
	  
	  my $mapped_feature = createFeature ($feature, @new_ranges);
	  
	  # Add a message informing than the feature mapping needs to be refines - just reporting the region
	  # give the number of exons
	  
	  # ...
	  
	  # check if it makes sense
	  # ie if fisrt exon has a start codon
	  # and the last exon the stop codon
	  # if not accurate, add a note message
	  
	  # deactivated !!
	  # should be reported in the report file

	  #my $newseq_obj = $new_seq_objects{$seq_id};
	  #if ($mapped_feature->key =~ /cds/i) {
	    # $mapped_feature = check_gene_accurracy ($mapped_feature, $newseq_obj);
	  #}
	  
	  return ($mapped_feature, $best_hsps_ref);
	  
	}
	
	$queryFile = setQueryFile ($querySequenceObject, $blast_directory);
	
	my $result = doBlast ($queryFile, $subjectFile, $blastOptions, 0);
	
	if (not defined $result) {
	  print STDERR "Blast search didn't work !!!\n";
	  return (undef, undef);
	}
	
	# as it's submapping, no extension and conflict processings

	# my ($sub_mapped_range, $sub_has_been_reversed, $sub_seq_id) = processResults ($filtering, $result);
	my ($sub_best_hsps_ref, $sub_mapped_range, $sub_has_been_reversed) = processResults ($filtering, $result, 0, 0);
	
	if (not defined $sub_mapped_range) {
	  if ($debug) {
	    print STDERR "sub mapping has failed - couldn't map the exon $index!!\n";
	  }
	  return (undef, undef);
	}
	else {
	  # Get the range mapped to the all sequence
	  
	  my $sub_mapped_range_to_all_sequence = setSubRangeMappedToAllSequence ($mapped_range, $sub_mapped_range, $strand);
	  push (@mapped_ranges, $sub_mapped_range_to_all_sequence);
	}
	
	$index++;
	
      }
      
      # Sub range mapping done
      # Create the new feature
      
      my $mapped_feature = createFeature ($feature, @mapped_ranges);
      
      if ($debug) {
	print STDERR "Feature ranges: " . Dumper (@mapped_ranges) . "\n";
      }
      
      # check if it makes sense
      # ie if fisrt exon has a start codon
      # and the last exon the stop codon
      # if not accurate, add a note message
      
      # deactivated !!
      # should be reported in the report file
      
      #if ($mapped_feature->key =~ /cds/i) {
	# $mapped_feature = check_gene_accurracy ($mapped_feature, $newseq_obj);
      #}
      
      return ($mapped_feature, $best_hsps_ref);
      
    }
  }
  else {
    # no range, return without any processing the given feature
    print STDERR "Feature without range !!!!\n";
    return (undef, undef);
  }
}

######################################################################
# set the subsequence object from a sequence object, giving a range object
# take always the forward sequence, because the mapping is always done on the forward strand
######################################################################

sub getSubSequenceObject {
  my ($range, $sequenceObject, $extension) = @_;
  
  my $extended_start = $range->start-$extension;
  if ($extended_start < 0) {
    $extended_start = 1;
  }
  
  my $extended_end = $range->end+$extension;
  if ($extended_end > length ($sequenceObject->str)) {
    $extended_end = length ($sequenceObject->str);
  }

  my $subSequenceObject = $sequenceObject->subseq ($extended_start, $extended_end, 0);
  $subSequenceObject->id ($sequenceObject->id . ".subSeq");
  return $subSequenceObject;
}

sub setQueryFile {
  my ($querySequenceObject, $blast_directory) = @_;
  my $queryFileName = "$blast_directory/query.fasta";
  
  my $outFasta = Bio::PSU::SeqFactory->make(
					    -file => ">$queryFileName",
					    -format => 'fasta',
					   );
  $outFasta->write_seq ($querySequenceObject);
  
  return $queryFileName;
}


sub setSubjectFile {
  my ($subjectSequenceObject, $blast_directory) = @_;
  my $subjectFileName = "$blast_directory/subject.fasta";
  
  my $outFasta = Bio::PSU::SeqFactory->make(
					    -file => ">$subjectFileName",
					    -format => 'fasta',
					   );
  $outFasta->write_seq ($subjectSequenceObject);
  
  # index the database
  qx "$_indexingdb $subjectFileName 2>/dev/null";
  
  return $subjectFileName;
}

sub getNewSequenceObjects {
  my ($newseqfile, $newseqformat, $blast_directory) = @_;
  my %sequenceObjects;
  
  my $subjectDb_file  = $newseqfile;
  my $newseqfile_name;
  # Check if the path is included in the file name, if yes, just keep the file name
  if ($newseqfile =~ /\//) {
    my @dirs = split ('/', $newseqfile);
    $newseqfile_name = $dirs[@dirs-1];
  }
  
  my $outSeqFactory  = undef;
  
  if (not $newseqformat =~ /fasta/i) {
    # generate a fasta database
    
    $subjectDb_file = "$blast_directory/subjectDb.fa";
    
    $outSeqFactory = Bio::PSU::SeqFactory->make(
						-file   => ">$subjectDb_file",
						-format => "fasta"
					       );
  }
  
  my $inSeqFactory = Bio::PSU::SeqFactory->make(
						-file   => "<$newseqfile",
						-format => "$newseqformat"
					       );
  
  while (my $sequence = $inSeqFactory->next_seq) {
    
    if (defined $sequence->id) {
      $sequenceObjects{$sequence->id} = $sequence;
    }
    else {
      print STDERR "ERROR the subject sequences MUST have a sequence ID!!\n";
      exit 1;
    }
    
    if (not $newseqformat =~ /fasta/i) {
      $outSeqFactory->write_seq ($sequence);
    }
  }
  
  if (keys (%sequenceObjects) == 0) {
    print STDERR "can\'t get the new sequence objects from file, $newseqfile!\n";
    print STDERR "Are you sure you specified the right format, $newseqformat?\n";
    exit 1;
  }
  
  if ($newseqformat =~ /fasta/i) {
    if (defined $newseqfile_name) {
      qx/ln -fs $newseqfile $blast_directory\/$newseqfile_name/;
      $subjectDb_file = "$blast_directory/$newseqfile_name";
    }
    else {
      my $curdir = qx/pwd/;
      chomp $curdir;
      
      qx/ln -fs $curdir\/$newseqfile $blast_directory\/$subjectDb_file/;
      $subjectDb_file = "$blast_directory/$subjectDb_file";
    }
  }
  qx "$_indexingdb $subjectDb_file 2>/dev/null";
  
  if ($debug) {
    print STDERR "Subject Database, $subjectDb_file done!!\n";
  }
  
  return ($subjectDb_file, %sequenceObjects);
}


sub getOverlappingFeatures {
  my ($mapped_feature, @newseq_features) = @_;
  my %overlaping_features;
  
  my $i = 0;
  while ($i < @newseq_features) {
    my $newseq_feature = $newseq_features[$i];
    
    # if full overlap
    
    #    if (fullOverlap ($mapped_feature, $newseq_feature)) {

    # overlap

    if ($mapped_feature->overlaps ($newseq_feature)) {

      if ($mapped_feature->key eq $newseq_feature->key) {
	$overlaping_features{$i} = $newseq_feature;
      }
    }
    $i++;
  }
  
  # return a hash whose key is the array index and the value is the feature
  
  if (keys (%overlaping_features) > 0) {
    return \%overlaping_features;
  }
  else {
    return undef;
  }
}

####################################################

## in a function 

# check if the mapped_feature fully overlaps one of the newseq_features
# and if it is the same key !

# @mapped_qualifiers_to_keep array list to set
# @newseq_qualifiers_to_keep array list to set

# if yes, then merge the qualifiers
# if overlaps but not fully report the new gene name

# if ($mapped_feature->overlaps (@newseq_features)) {
# found overlaps !
# }
###########################################


sub mergeFeatures {
  my ($mapped_feature, $overlaped_feature, $mapped_qualifiers_to_keep_ref, $newseq_qualifiers_to_keep_ref) = @_;
  my $merged_feature = Bio::PSU::Feature->new (
					       -key    => $mapped_feature->key,
					       -ranges => [$mapped_feature->ranges]
					      );

  foreach my $qualifier_key ($mapped_feature->qnames) {
    if (@$mapped_qualifiers_to_keep_ref > 0) {
      if ($mapped_qualifiers_to_keep_ref->[0] eq 'all') {
	$merged_feature->qadd ($qualifier_key, $mapped_feature->qvalues('^' . $qualifier_key . '$'));
      }
      else {
        my $exists = is_in ($mapped_qualifiers_to_keep_ref, $qualifier_key);
        if ($exists) {
          $merged_feature->qadd ($qualifier_key, $mapped_feature->qvalues('^' . $qualifier_key . '$'));
        }
      }
    }
  }
  
  foreach my $qualifier_key ($overlaped_feature->qnames) {
    if (@$newseq_qualifiers_to_keep_ref > 0) {
      if ($newseq_qualifiers_to_keep_ref->[0] eq 'all') {
        $merged_feature->qadd ($qualifier_key, $overlaped_feature->qvalues('^' . $qualifier_key . '$'));
      }
      else {
        my $exists = is_in ($newseq_qualifiers_to_keep_ref, $qualifier_key);
        if ($exists) {
	  $merged_feature->qadd ($qualifier_key, $overlaped_feature->qvalues('^' . $qualifier_key . '$'));
        }
      }
    }
  }
  
  return $merged_feature;
}


sub is_in { 
  my ($mapped_qualifiers_to_keep_ref, $qualifier_key) = @_;
  my $exists = 0;
  
  foreach my $qualifier (@$mapped_qualifiers_to_keep_ref) {
    if ($qualifier eq $qualifier_key) {
      $exists = 1;
    }
  }
  
  return $exists;
}

# Get the PSU Gene Name 
# if the Feature is a CDS

sub getGeneName {
  my ($feature) = @_;
  my $id;

  eval {
    $id = $id_helper->get_systematic_id($feature);
  };

  if (!defined $id) {
    warn "feature at " . $feature->start . ".." . $feature->end .
      " has no gene name\n";
  }
  
  return $id;
}


sub fullOverlap {
  my ($mapped_feature, $newseq_feature) = @_;
  
  if ($mapped_feature->strand != $newseq_feature->strand) {
    return 0;
  }
  else {
    my @old_ranges    = $mapped_feature->ranges;
    my @newseq_ranges = $newseq_feature->ranges;
    if (@old_ranges != @newseq_ranges) {
      return 0;
    }
    else {
      my $i = 0;
      while ($i < @old_ranges) {
	my $old_range    = $old_ranges[$i];
	my $newseq_range = $newseq_ranges[$i];
		if (($old_range->start != $newseq_range->start) || ($old_range->end != $newseq_range->end)) {
	  return 0;
	}
	$i++;
      }
      return 1;
    }
  }

}

