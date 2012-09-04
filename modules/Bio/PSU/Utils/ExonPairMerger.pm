package Bio::PSU::Utils::ExonPairMerger;

###################################
# version with transcripts exon pairs processing
# output :
#   * GFF  => geneid processing
#   * Embl => Artemis 
###################################

use strict;

use Cwd;
##########################################

use Data::Dumper;

use Carp;

# PSU modules for EMBL writing
use Bio::PSU::SeqFactory;
use Bio::PSU::Seq;
use Bio::PSU::Feature;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name=>value pairs (odd number supplied)";
  }

  my $self     = {
		  input               => undef,
		  curdir              => undef,
                  max_intron_size     => 400,
		  verbose             => 0,
		  sembl               => undef,
		  @_,
		 };

  my $sembl = Bio::PSU::SeqFactory->make(
					 -file   => $self->{input},
					);
  $self->{sembl} = $sembl;

  bless ($self, $class);
  return $self;
}

sub getCurdir {
  my $self = shift;
  return $self->{curdir};
}

sub setCurdir {
  my $self = shift;
  my ($curdir) = @_;
  
  $self->{curdir} = $curdir;
}

# @return the predicted gene features after merging 

sub next_sequence {
  my $self = shift;
  my $seq_str = "";
  if (@_ > 0) {
    $seq_str = shift;
  }
  my $max_intron_size = $self->{max_intron_size};
  my $_verbose        = $self->{verbose};
  my $sembl           = $self->{sembl};

  my ($sequence_name, $genes, $transcripts) = getExonerateFeaturesFromEMBL ($sembl, $_verbose);

  if (defined $sequence_name) {
    my @new_features = ();
    my @new_genes = processGenes ($_verbose, @$genes);
    my @new_transcripts = processTranscripts ($_verbose, @$transcripts);
    
    # process them to do a more complete merging !!
    # ie even if they don't overlap but genes predicted from the same ESTs => merge the genes !
    
    #@new_genes = processGeneMerging ($max_intron_size, @new_genes);
    
    push (@new_features, @new_genes);
    push (@new_features, @new_transcripts);
    
    if ($_verbose) {
      print STDERR "nb generated features: ". @new_features . "\n";
    }

    # ??
    #my @exonFeaturesPerGene = generateExonFeaturesPerEntity (@new_genes);
    #$self->{exonFeaturesPerGene} = \@exonFeaturesPerGene;
    
    return ($sequence_name, \@new_features);
  }
  else {
    return (undef, undef);
  }
}

##################################################################
sub getExonerateFeaturesFromEMBL {
  my ($sembl, $_verbose) = @_;

  my $seqname = undef;
  my @genes = ();
  my @transcripts = ();

  if (my $seqobj = $sembl->next_seq) {
    $seqname = $seqobj->id;

    if ($_verbose) {
      print STDERR "parsing sequence, $seqname...\n";
    }

    foreach my $feature ($seqobj->features) {
      if ($feature->key eq "gene") {
	push (@genes, $feature);
      }
      elsif ($feature->key eq "mRNA") {
	push (@transcripts, $feature);
      }
    }

    return ($seqname, \@genes, \@transcripts);

  }
  else {
    return (undef, undef, undef);
  }

}


sub processGenes {
  my ($_verbose, @genes) = @_;

  if ($_verbose) {
    print STDERR "gene list for overlapping merging process: " . @genes . "\n";
  }

  # Gene Features merging into a unique gene model

  # loop on the genes minus those who have been merged
  
  my @new_genes = ();
  
  while (@genes) {
    my $gene1 = splice (@genes,0,1);
    my $i = 0;
    my $merging_number = 1;
    my $nb_loops = 0;
    my $nb_genes = @genes;
    
    # loop on @genes minus gene1
    
    while ($i < $nb_genes) {
      # remove the first gene of the gene array and process it
      my $gene2  = splice (@genes,0,1);
      if (can_do_complete_merge ($gene1, $gene2)) {
	
	$gene1 = merge ($gene1, $gene2, "gene");
	
	# how many merges have been done so far on this gene
	$merging_number++;
	
	$i++;
      }
      else {
	# reprocess gene2 later => add it at the end of the array
	push (@genes, $gene2);
	# incremente $i if a gene has been already processed, ie added at the end of the array, ie when the set of genes has been already processed once
	# otherwise infinite loop
	if ($nb_loops >= ($nb_genes)) {
	  $i++;
	}
      }
      $nb_loops++;
    }
    
    # ponderate the score and the percent_id by the number of merged genes
  
    if ($gene1->qexists ('score')) {
      my @a_tmp = $gene1->qvalues ('score');
      my $score = $a_tmp[0];
      $score = $score/$merging_number;
      $gene1->qremove ('score');
      $gene1->qadd ('score', $score);
    }
    if ($gene1->qexists ('percent_id')) {
      my @a_tmp = $gene1->qvalues ('percent_id');
      my $percent_id = $a_tmp[0];
      $percent_id = $percent_id/$merging_number;
      $gene1->qremove ('percent_id');
      $gene1->qadd ('percent_id', $percent_id);
    }
    if ($gene1->qexists ('ESTs')) {
      my @values = $gene1->qvalues ('ESTs');
      my $name1  = $values[0];
      if ($_verbose) {    
        # print STDERR "\n\tadd new gene - $name1 - into the genes list\n\n";
      }
    }

    push (@new_genes, $gene1);
  }

  if ($_verbose) {
    print STDERR "nb genes after processing: " . @new_genes . "\n";
  }

  return @new_genes;
} 

 
sub processTranscripts {
  my ($_verbose, @transcripts) = @_;

  if ($_verbose) {
    print STDERR "transcripts list for overlapping merging process: " . @transcripts . "\n";
  }

  # Transcript Features merging by exon pairs processing
  
  # loop on the transcripts minus those who have been merged
  
  my @new_transcripts = ();
  
  while (@transcripts) {
    my $transcript1 = splice (@transcripts,0,1);
    my $i = 0;
    my $merging_number = 1;
    my $nb_loops = 0;
    my $nb_transcripts = @transcripts;
    
    # loop on @transcripts minus transcript1
    
    while ($i < $nb_transcripts) {
      # remove the first transcript of the transcript array and process it
      my $transcript2  = splice (@transcripts,0,1);
      if (can_merge ($transcript1, $transcript2)) {
	
	$transcript1 = merge ($transcript1, $transcript2, "mRNA");
	
	# how many merges have been done so far on this transcript
	$merging_number++;
	
	$i++;
      }
      else {
	# reprocess transcript2 later => add it at the end of the array
	push (@transcripts, $transcript2);
	# incremente $i if a transcript has been already processed, ie added at the end of the array, ie when the set of transcripts has been already processed once
	# otherwise infinite loop
	if ($nb_loops >= ($nb_transcripts)) {
	  $i++;
	}
      }
      $nb_loops++;
    }
    
    # ponderate the score and the percent_id by the number of merged transcripts
    
    if ($transcript1->qexists ('score')) {
      my @a_tmp = $transcript1->qvalues ('score');
      my $score = $a_tmp[0];
      $score = $score/$merging_number;
      $transcript1->qremove ('score');
      $transcript1->qadd ('score', $score);
    }
    if ($transcript1->qexists ('percent_id')) {
      my @a_tmp = $transcript1->qvalues ('percent_id');
      my $percent_id = $a_tmp[0];
      $percent_id = $percent_id/$merging_number;
      $transcript1->qremove ('percent_id');
      $transcript1->qadd ('percent_id', $percent_id);
    }
    if ($transcript1->qexists ('ESTs')) {
      my @values = $transcript1->qvalues ('ESTs');
      my $name1  = $values[0];
      if ($_verbose) {
        # print STDERR "\n\tadd new transcript - $name1 - into the transcripts list\n\n";
      }
    }
    
    push (@new_transcripts, $transcript1);
  }

  if ($_verbose) {
    print STDERR "nb transcripts after processing: " . @new_transcripts . "\n";
  }
  
  return @new_transcripts;
}


sub processGeneMerging {
  my ($max_intron_size, @genes) = @_;
  my @new_genes = ();

  print STDERR "nb genes for sharing ESTs merging processing: " . @genes . "\n\n";

  foreach my $gene (@genes) {

    print STDERR "current gene: " . Dumper ($gene->ranges) . "\n";

    if (are_sharing_exons ($max_intron_size, $gene, @new_genes)) {
      print STDERR "sharing for current gene!!!\n";
      my ($new_genes_tmp_ref, $sharing_exons_genes_ref) = get_sharing_exons_genes ($max_intron_size, $gene, @new_genes);
      my @new_genes_tmp = @$new_genes_tmp_ref;
      my @sharing_exons_genes = @$sharing_exons_genes_ref;

      print STDERR "new genes tmp list: " . @new_genes_tmp . "\n";
      print STDERR "sharing exon genes number: " . (@sharing_exons_genes+1) . "\n";

      # if @sharing size == 0, don't need to add $gene !!!
      # because nothing to merge it with !!

      push (@sharing_exons_genes, $gene);

      my $new_gene = mergeGenes (@sharing_exons_genes);
      push (@new_genes_tmp, $new_gene);
      @new_genes = @new_genes_tmp;
    }
    else {
      print STDERR "no sharing for current gene!!!\n";
      push (@new_genes, $gene);
    }
  }

  print STDERR "nb genes after processing: " . @new_genes . "\n";

  return @new_genes;
}


sub are_sharing_exons {
  my ($max_intron_size, $gene1, @genes) = @_;
  my @ests_gene1 = $gene1->qvalues ('ESTs');
  my $ests_gene1 = $ests_gene1[0];
  @ests_gene1 = split (',', $ests_gene1);
  # print STDERR "ESTs gene 1 list: @ests_gene1\n";

  foreach my $gene2 (@genes) {
    my @ests_gene2 = $gene2->qvalues ('ESTs');
    my $ests_gene2 = $ests_gene2[0];
    if (matches ($ests_gene2, @ests_gene1)) {
      if (intron_size_good ($max_intron_size, $gene1, $gene2)) {
        return 1;
      }
      else {
        print STDERR "Sharing exons but intron size too big !!!\n";
      }
    }
  }

  return 0;
}

sub get_sharing_exons_genes {
  my ($max_intron_size, $gene1, @genes) = @_;
  my @new_genes = ();
  my @sharing_exons_genes = (); 
  
  my @ests_gene1 = $gene1->qvalues ('ESTs');
  my $ests_gene1 = $ests_gene1[0];
  @ests_gene1 = split (',', $ests_gene1);

  foreach my $gene2 (@genes) {
    my @ests_gene2 = $gene2->qvalues ('ESTs');
    my $ests_gene2 = $ests_gene2[0];
    if (matches ($ests_gene2, @ests_gene1)) {
      print STDERR "share some ESTs\n";
      if (intron_size_good ($max_intron_size, $gene1, $gene2)) {
        print STDERR "intron size ok!\n";
        push (@sharing_exons_genes, $gene2);
      }
      else {
        push (@new_genes, $gene2);
      }
    }
    else {
      push (@new_genes, $gene2);
    }
  }

  print STDERR "nb genes found for merging: " . @sharing_exons_genes . "\n";

  return (\@new_genes, \@sharing_exons_genes); 
}

sub matches {
  my ($ests_gene2, @ests_gene1) = @_;
  
  foreach my $est_gene1 (@ests_gene1) {
    if ($ests_gene2 =~ /$est_gene1/) {
      return 1;
    }
  }
  
  return 0;
}


sub intron_size_good {
  my ($max_intron_size, $gene1, $gene2) = @_;
  my @ranges1 = $gene1->ranges; # only 1 exon !!
  my $range1  = $ranges1[0];
  my @ranges2 = $gene2->ranges;
  my $downstream = 1;

  # print STDERR "gene1: " . Dumper ($gene1) . "\n";
  # print STDERR "gene2: " . Dumper ($gene2) . "\n";

  while (my $range2 = shift (@ranges2)) {
    # if gene1 upstream current exon
    if ($range1->end < $range2->start) {
      my $intron_size = $range2->start - $range1->end;
      if ($intron_size <= $max_intron_size) {
        if ($downstream == 1) {
          return 1;
        }
      }
      else {
        return 0;
      }
    }
    # if gene1 downstream current exon
    elsif ($range2->end < $range1->start) {
      my $intron_size = $range1->start - $range2->end;
      if ($intron_size <= $max_intron_size) {
        $downstream = 1;
      }
      else {
        $downstream = 0;
      }
    }
    else {
      print STDERR "not supposed to come here !!!!\n";
      print STDERR "means some overlapping...\n";
    }
  }

  return $downstream;
}


sub mergeGenes {
  my (@sharing_exons_genes) = @_;
  my @new_ranges = ();
 
  print STDERR "nb genes to merge: " . @sharing_exons_genes . "\n";

  my $nb_genes = @sharing_exons_genes;
  my $gene = pop (@sharing_exons_genes);

  push (@new_ranges, $gene->ranges);
  my @ESTs_tmp = $gene->qvalues ('ESTs');
  my $ESTs = $ESTs_tmp[0];
  my @scores_tmp = $gene->qvalues ('score');
  my $score = $scores_tmp[0];
  my @percent_ids_tmp = $gene->qvalues ('percent_id');
  my $percent_id = $percent_ids_tmp[0];
  my @similarities = $gene->qvalues ('similarity');
  
  foreach $gene (@sharing_exons_genes) {
    push (@new_ranges, $gene->ranges);
    push (@similarities, $gene->qvalues ('similarity'));
    @scores_tmp = $gene->qvalues ('score');
    $score = $score + $scores_tmp[0];
    $percent_id = $percent_id + $percent_ids_tmp[0];
    # get a unique EST list !!
    @ESTs_tmp = $gene->qvalues ('ESTs');
    my $ESTs2 = $ESTs_tmp[0];
    @ESTs_tmp = split (',', $ESTs2);
    print STDERR "ESTs_tmp array: @ESTs_tmp\n";
    foreach my $est_tmp (@ESTs_tmp) {
      $est_tmp =~ s/\s//g;
      if (not $ESTs =~ /$est_tmp/) {
        $ESTs = $ESTs . ", $est_tmp";
      }
    }
  }

  $score = $score / $nb_genes;
  $percent_id = $percent_id / $nb_genes;

  # sort ranges

  @new_ranges = sortRanges (@new_ranges);

  print STDERR "nb ranges: " . @new_ranges . "\n";
  print STDERR "ranges: " . Dumper (@new_ranges) . "\n";

  my $new_gene = Bio::PSU::Feature->new (
                                         -key => 'gene',
                                         -ranges => [@new_ranges],
                                        );
  $new_gene->qadd ('ESTs', $ESTs);
  $new_gene->qadd ('similarity' , @similarities);
  $new_gene->qadd ('score', $score);
  $new_gene->qadd ('percent_id', $percent_id);

  return $new_gene;
}


sub writeEmblfeatures {
  my ($curdir, $emblOutputFile, $seq_str, @emblFeatures) = @_;

  # done => attach the generated features to a sequence object

  #my $new_seq_obj = Bio::PSU::Seq->new (
  #					-type => 'dna',
  #					-str  => $seq_str
  #				       );

  my $new_seq_obj = Bio::PSU::Seq->new (
                                        -type => 'dna',
                                       );

  print STDERR "nb features to attach: " . @emblFeatures . "\n";
  # print STDERR "Features to attach: " . Dumper (@emblFeatures) . "\n";

  $new_seq_obj->features (@emblFeatures);
  
  # copy this seq obj into a new file
  
  # my $embl_output_file_path = $curdir . '/' . $emblOutputFile;
  # if path already included !!!
  my $embl_output_file_path = $emblOutputFile;

  my $new_output_obj = undef;
  
  # always in an EMBL format
  
  $new_output_obj = Bio::PSU::SeqFactory->make(
					       -file   => ">$embl_output_file_path",
					       -format => 'embl'
					      );
  $new_output_obj->write_seq($new_seq_obj);
}

# generate a set of GFF exon features from the set of features (gene features or transcript features)

sub generateExonsFeatures {
  my $self = shift;
  my ($sequence_name, @set) = @_;
  my @exonsFeatures = ();
  my $source = $self->{searchType};
  my $score = ".";
  my $frame = ".";
  
  # sort based on feature->ranges[0]->start value
  # so assume that foreach feature, the ranges are already sorted
  
  @set = sortFeatures (@set);
    
  foreach my $feature (@set) {
    my @ranges = sortRanges ($feature->ranges);
    my @values = $feature->qvalues ('ESTs');
    my $ESTs   = $values[0];
    my $attributes = "ESTs \"$ESTs\"";
    my $i = 0;
    
    foreach my $range (@ranges) {
      my $start  = $range->start;
      my $end    = $range->end;
      my $strand = $range->strand;
      
      # convert Exon feature into First - Internal - Terminal

      my $feature_name = "";

      # if a unique exon => Terminal
      
      if (@ranges == 1) {
	$feature_name = "Terminal";
      }
      elsif ($i != 0 && $i != @ranges) {
	$feature_name = "Internal";
      }
      else {
	# more complicate because can be an unfinished sequence !!
	# and so the last or first exon can be actually an Internal exon

	if (($i==0 && $strand==1) || ($i==@ranges-1 && $strand==-1)) { 
	  $feature_name = "First";
	}
	else {
	  $feature_name = "Terminal";
	}
      }
      
      if ($strand == 1) {
	$strand = "+";
      }
      else {
	$strand = "-";
      }
      
      my $gff_feature = "$sequence_name\t$source\t$feature_name\t$start\t$end\t$score\t$strand\t$frame\t$attributes";
      push (@exonsFeatures, $gff_feature);

      # if First or Last Exon, add another feature - Internal Exon
      
      # if ($feature_name =~ /First|Terminal/) {
      # $feature_name = "Internal";
      # $gff_feature = "$sequence_name\t$source\t$feature_name\t$start\t$end\t$score\t$strand\t$frame\t$attributes";
      # push (@exonsFeatures, $gff_feature);
      # }
      
      $i++;
    }
  }

  return @exonsFeatures;
}


# generate a set of Genomewise exon features from a set of features (gene features or transcript features)
# syntax : exon start end
# classify by the strand

sub generateExonFeaturesPerEntity {
  my (@set) = @_;
  my @exonFeaturesPerEntity = ();
  
  # sort based on feature->ranges[0]->start value
  # so assume that foreach feature, the ranges are already sorted
  
  @set = sortFeatures (@set);
    
  foreach my $feature (@set) {
    my @ranges = sortRanges ($feature->ranges);
    my $strand = 0;
    my @exons = ();

    foreach my $range (@ranges) {
      my $start  = $range->start;
      my $end    = $range->end;
      $strand = $range->strand;
      my $exonFeature = "exon\t$start\t$end\n";
     
      push (@exons, $exonFeature);
    }

    my %entity = (
		  'exons' => \@exons,
		  'strand'  => $strand,
		 );
    
    push (@exonFeaturesPerEntity, \%entity);
  }

  # print STDERR "dump: " . Dumper (@exonFeaturesPerEntity) . "\n";

  return @exonFeaturesPerEntity;
}


sub writeExonsFeatures {
  my ($gffOutputFile, @exonsFeatures) = @_;

  # copy these features into a GFF output file

  open (OUTGFF, ">$gffOutputFile") or die "can't open file, $gffOutputFile\n";

  foreach my $feature (@exonsFeatures) {
    print OUTGFF "$feature\n";
  }

  close OUTGFF;
}


sub getGeneRange {
  my (@ranges) = @_;
  my $start = 1000000000000;
  my $end = -1;
  my $strand = 0;

  foreach my $range (@ranges) {
    $strand = $range->strand();
    $start  = min ($start, $range->start);
    $end    = max ($end, $range->end);
  }
  
  # print STDERR "start, end, strand: $start, $end, $strand\n";

  my $range = Bio::PSU::Range->new(-start       => $start,
		                   -end         => $end,
				   -strand      => $strand
				  );
  return $range;
}


sub getRange {
  my ($range1, $range2) = @_;

  my $strand = $range1->strand;
  my $start  = min ($range1->start, $range2->start);
  my $end    = max ($range1->end, $range2->end);

  # print STDERR "start, end, strand: $start, $end, $strand\n";

  my $range = Bio::PSU::Range->new(-start       => $start,
				   -end         => $end,
				   -strand      => $strand
				  );
  return $range;
}


sub min {
  my ($int1, $int2) = @_;
  my $min = -1;

  if ($int1 < $int2) {
    $min = $int1;
  }
  else {
    $min = $int2;
  }

  return $min;
}


sub max {
  my ($int1, $int2) = @_;
  my $max = -1;

  if ($int1 > $int2) {
    $max = $int1;
  }
  else {
    $max = $int2;
  }

  return $max;
}


sub remove {
  my ($index, @genes) = @_;
  my @new_genes = ();
  my $i = 0;

  while (@genes) {
    my $gene = pop (@genes);
    if ($i == $index) {
      push (@new_genes, @genes);
      last;
    }
    else {
      push (@new_genes, $gene);
    }
    $i++;
  }

  return @new_genes;
}


sub can_merge {
  my ($tr1, $tr2) = @_;

  # true if (
  #          tr1 and tr2 share at least one overlaping exon
  #          and no 5' gap
  #          and no 3' gap
  #         )

  if (overlap ($tr1, $tr2)
      && (not (fiveprimegap ($tr1, $tr2)))
      && (not (threeprimegap ($tr1, $tr2)))
     ) {
    # print STDERR "yes\n";
    return 1;
  }
  else { 
    # print STDERR "no\n";
    return 0; 
  }
}

# gap or not gap, do the merging. 
# Consequence: generate a unique gene structure

sub can_do_complete_merge {
  my ($tr1, $tr2) = @_;

  # true if (
  #          tr1 and tr2 share at least one overlaping exon
  #         )

  if (overlap ($tr1, $tr2)) {
    # print STDERR "yes\n";
    return 1;
  }
  else {
    # print STDERR "no\n";
    return 0;
  }
}


sub merge {
  my ($transcript1, $transcript2, $feature_name) = @_;

  my @ranges1 = $transcript1->ranges;
  my @ranges2 = $transcript2->ranges;
  my @new_ranges = getMergedRanges (\@ranges1, \@ranges2);
  
  # update score, percent_id and ESTs information
  
  # score
  my $score = undef;
  if ($transcript1->qexists ('score') && $transcript2->qexists ('score')) {
    my @a_tmp = $transcript1->qvalues ('score');
    my $score1 = $a_tmp[0];
    @a_tmp = $transcript2->qvalues ('score');
    my $score2 = $a_tmp[0];
    $score = $score1 + $score2;
  }
  # percent_id
  my $percent_id = undef;
  if ($transcript1->qexists ('percent_id') && $transcript2->qexists ('percent_id')) {
    my @a_tmp = $transcript1->qvalues ('percent_id');
    my $percent1 = $a_tmp[0];
    @a_tmp = $transcript2->qvalues ('percent_id');
    my $percent2 = $a_tmp[0];
    $percent_id = $percent1 + $percent2;
  }

  # ESTs sequences

  my $est1_str_list = undef;
  if ($transcript1->qexists ('ESTs') && $transcript2->qexists ('ESTs')) {
    my @a_tmp = $transcript1->qvalues ('ESTs');
    $est1_str_list = $a_tmp[0];
    @a_tmp = $transcript2->qvalues ('ESTs');
    my $est2_str_list = $a_tmp[0];
    my @ests2 = split (',', $est2_str_list);

    # print STDERR "ests2 list: @ests2\n";

    # if the contig, est2, is not in the est list, add it
    foreach my $est2 (@ests2) {
      if (not $est1_str_list =~ /$est2/) {
        $est1_str_list = "$est1_str_list, $est2";
      }
    }
  }

  # note qualifiers

  my $note = undef;
  if ($transcript1->qexists ('note') && $transcript2->qexists ('note')) {
    my @notes1 = $transcript1->qvalues ('note');
    my @notes2 = $transcript2->qvalues ('note');
  
    if ($notes1[0] =~ /feature predicted by exonerate:est2genome/ && $notes2[0] =~ /feature predicted by exonerate:est2genome/) {
      $note = $notes1[0];
    }
    else {
      print STDERR "WARNING: feature origins don't match!!\n";
    }
  }

  # method qualifier

  my $method = undef;
  if ($transcript1->qexists ('method') && $transcript2->qexists ('method')) {
    my @methods1 = $transcript1->qvalues ('method');
    my @methods2 = $transcript2->qvalues ('method');
  
    if ($methods1[0] =~ /exonerate:est2genome/ && $methods2[0] =~ /exonerate:est2genome/) {
      $method = $methods1[0];
    }
    else {
      print STDERR "WARNING: Feature origins don't match!!\n";
    }
  }

  # the new transcript feature
  
  my $new_transcript = Bio::PSU::Feature->new(
					      -key    => "$feature_name",
					      -ranges => [@new_ranges]
					     );
  if (defined $est1_str_list) {
    $new_transcript->qadd ('ESTs', $est1_str_list);
  }
  if (defined $percent_id) {
    $new_transcript->qadd ('percent_id', $percent_id);
  }
  if (defined $score) {
    $new_transcript->qadd ('score', $score);
  }
  if (defined $note) {
    $new_transcript->qadd ('note', $note);
  }
  if (defined $method) {
    $new_transcript->qadd ('method', $method);
  }

  return $new_transcript;
}


sub overlap {
  my ($transcript1, $transcript2) = @_;
  my @ranges1 = $transcript1->ranges;
  my @ranges2 = $transcript2->ranges;
  
  # if on different strand - they don't overlap
  my $range1 = $ranges1[0];
  my $range2 = $ranges2[0];
  if ($range1->strand != $range2->strand) {
    return 0;
  }

  # print STDERR "ranges in overlap: @ranges2\n";

  foreach my $range (@ranges1) {
    if ($range->overlaps (@ranges2)) {
      return 1;
    }
  }
  return 0;
}

# whatever any exon pair, true if there is a gap upstream from them

sub fiveprimegap {
  my ($tr1, $tr2) = @_;
  
  my @ranges1 = sortRanges ($tr1->ranges);
  my @ranges2 = sortRanges ($tr2->ranges);
  my $i1 = 0;

  # print STDERR "ranges in fiveprimegap: @ranges2\n";

  while ($i1 < @ranges1) {
    my $range1 = $ranges1[$i1];
    if ($range1->overlaps (@ranges2)) {
      my $i2 = 0;
      while ($i2 < @ranges2) {
	# print STDERR "i2, range[i2] in fiveprimegap: $i2, " . $ranges2[$i2] . "\n";
	if ($range1->overlaps ($ranges2[$i2])) {

	  my $j1 = $i1-1;
	  my $j2 = $i2-1;
	  while ($j1>=0 && $j2>=0) {
	    # print STDERR "range[j2] in fiveprimegap: " . $ranges2[$j2] . "\n";
	    if (not $ranges1[$j1]->overlaps ($ranges2[$j2])) {	      
	      return 1;
	    }
	    $j1--;
	    $j2--;
	  }
	}
	$i2++;
      }
    }
    $i1++;
  }
  # no 5' gap - actually because no exon pair or because they all overlap each other upstream any exon pair
  return 0;
}

# whatever two overlaping exons, true if there is a gap downstream from them

sub threeprimegap {
  my ($tr1, $tr2) = @_;
  
  my @ranges1 = sortRanges ($tr1->ranges);
  my @ranges2 = sortRanges ($tr2->ranges);

  my $i1 = 0;

  # print STDERR "ranges2 in threeprimegap: @ranges2\n";

  while ($i1 < @ranges1) {
    my $range1 = $ranges1[$i1];
    if ($range1->overlaps (@ranges2)) {
      my $i2 = 0;
      while ($i2 < @ranges2) {
	# print STDERR "range2[i2] in threeprimegap: " . $ranges2[$i2] . "\n";
	if ($range1->overlaps ($ranges2[$i2])) {
	  my $j1 = $i1+1;
	  my $j2 = $i2+1;
	  while ($j1<@ranges1 && $j2<@ranges2) {
	    # print STDERR "range[j2] in threeprimegap: " . $ranges2[$j2] . "\n";
	    if (not $ranges1[$j1]->overlaps ($ranges2[$j2])) {
	      return 1;
	    }
	    $j1++;
	    $j2++;
	  }
	}
	$i2++;
      }
    }
    $i1++;
  }
  # no 3' gap - actually because no exon pair or because they all overlap each other dowstream any exon pair
  return 0;
}


sub getMergedRanges {
  my ($ranges1_ref, $ranges2_ref) = @_;
  my @ranges1 = @$ranges1_ref;
  my @ranges2 = @$ranges2_ref;

  @ranges1 = sortRanges (@ranges1);
  @ranges2 = sortRanges (@ranges2);

  my $nb_ranges1 = @ranges1;
  my $nb_ranges2 = @ranges2;

  my @new_ranges = ();
  my $i1 = 0;
  my $i2 = 0;

  while ($i1<@ranges1 && $i2<@ranges2) {
    my $range1 = $ranges1[$i1];
    my $range2 = $ranges2[$i2];

    my $new_range = undef;
    
    if ($range1->overlaps($range2)) {
      # if range conflict - doesn't start and/or end at the same coordinate
      my $conflict = has_conflict (\@ranges1, \@ranges2, $i1, $i2);
      if ($conflict) {
	$new_range = getInternalRange (\@ranges1, \@ranges2, $i1, $i2);
      }
      else {
	# get min/max range
	$new_range = getRange ($range1, $range2);
      }
      $i1++;
      $i2++;
    }
    elsif ($range1->start < $range2->start) {
      $new_range = $range1;
      $i1++;
    }
    else {
      $new_range = $range2;
      $i2++;
    }

    push (@new_ranges, $new_range);
  }

  if ($i1 < $nb_ranges1) {
    push (@new_ranges, splice (@ranges1, $i1, $nb_ranges1));
  }
  else {
    push (@new_ranges, splice (@ranges2, $i2, $nb_ranges2));
  }

  # merge contiguous overlapping ranges within the same ranges

  @new_ranges = mergeOverlappingRanges (@new_ranges);

  return @new_ranges;
}


# Merge overlapping contiguous ranges within the same feature

sub mergeOverlappingRanges {
  my (@ranges) = @_;
  my $i = 0;
  my @new_ranges = ();

  while (@ranges) {
    my $range1      = splice (@ranges,0,1);
    if ($range1->overlaps (@ranges)) {
      my $range2    = splice (@ranges,0,1);
      my $start     = min ($range1->start, $range2->start);
      my $end       = max ($range1->end, $range2->end);
      my $strand    = $range1->strand;
      my $new_range = Bio::PSU::Range->new(
					   -start       => $start,
					   -end         => $end,
					   -strand      => $strand
					  );
      push (@new_ranges, $new_range);
    }
    else {
      push (@new_ranges, $range1);
    }
    $i++;
  }

  return @new_ranges;
}


sub sortRanges {
  my (@ranges) = @_;

  my @new_ranges = sort { $a->start <=> $b->start } @ranges;

  return @new_ranges;
}


sub has_conflict {
  my ($ranges1_ref, $ranges2_ref, $i1, $i2) = @_;
  my @ranges1 = @$ranges1_ref;
  my @ranges2 = @$ranges2_ref;
  my $conflict = 0;
  
  if (
      ($i1==@ranges1-1 xor $i2==@ranges2-1) || ($i1==0 xor $i2==0)
     ) {
    $conflict = 1;
  }
  
  return $conflict;
}


sub getInternalRange {
  my ($ranges1_ref, $ranges2_ref, $i1, $i2) = @_;
  my @ranges1 = @$ranges1_ref;
  my @ranges2 = @$ranges2_ref;
  
  my $range1 = $ranges1[$i1];
  my $range2 = $ranges2[$i2];

  my $start  = -1;
  my $end   = -1;
  my $strand = $range1->strand;

  # end coordinate

  if ($i1==@ranges1-1 xor $i2==@ranges2-1) {
    if ($i1!=@ranges1-1) {
      # if $range1 is internal and $range2 is not => end = $range1->end
      # else end = $range2->end
      $end = $range1->end;
    }
    else {
      $end = $range2->end;
    }
  }
  else {
    # they're both internal
    # otherwise min/max policy
    $end = max ($range1->end, $range2->end);
  }

  # start coordinate

  if ($i1==0 xor $i2==0) {
    if ($i1!=0) {
      # if $range1 is internal and $range2 is not => start = $range1->start
      # else start = $range2->start
      $start = $range1->start;
    }
    else {
      $start = $range2->start;
    }
  }
  else {
    # they're both internal
    # otherwise min/max policy
    $start = min ($range1->start, $range2->start);
  }

  my $new_range = Bio::PSU::Range->new(-start       => $start,
				       -end         => $end,
				       -strand      => $strand
				      );
  
  return $new_range;
}


sub sortFeatures {
  my (@old_features) = @_;
  my @new_features   = ();
  
  my $feature = undef;

  while (@old_features) {
    ($feature, @old_features) = minFeature (@old_features);
    push (@new_features, $feature);
  }

  return @new_features;
}

sub minFeature {
  my (@old_features) = @_;
  my $min_feature    = $old_features[0];
  my $index = 0;
  my $i = 1;

  while ($i < @old_features) {
    my $old_feature = $old_features[$i];
    my @old_ranges  = $old_feature->ranges;
    my $old_range   = $old_ranges[0];
    my @min_ranges  = $min_feature->ranges;
    my $min_range   = $min_ranges[0];
    if ($old_range->start < $min_range->start) {
      $min_feature = $old_feature;
      $index = $i;
    }
    $i++;
  }
  splice (@old_features, $index, 1);
  return ($min_feature, @old_features);
}

1;
