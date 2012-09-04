package Bio::PSU::Utils::ExonerateGFF;

use strict;
use Carp;
use Data::Dumper;

# PSU
use Bio::PSU::Feature;
use Bio::PSU::Range;

sub new {
  my $invocant = shift;
  my $class    = ref ($invocant) || $invocant;
  if(@_ % 2) {
    croak "Default options must be name=>value pairs (odd number supplied)";
  }
  
  my $self     = {
		  file      => undef,
		  fh        => undef,
		  type      => undef,
		  feature_key => undef,
		  debug       => 0,
		  add_gene_feature => 0,
		  seq_length      => undef,
		  seq_lengths_ref => undef,
		  fileindex => 0,
		  geneindex => 0,
		  @_,
		 };
  bless ($self, $class);
  
  if (not defined $self->{fh} and defined $self->{file}) {
    my $fh = IO::File->new ($self->{file});
    $self->{fh} = $fh;
  }
  
  return $self;
}

sub next_gene_prediction {
  my $self  = shift;
  my $fh    = $self->{fh};
  my $debug = $self->{debug};
  my $seq_length    = $self->{seq_length};
  my $sequence_type = $self->{type};
  my $feature_key   = $self->{feature_key};
  my $add_gene_feature = $self->{add_gene_feature};
  my $seq_lengths   = $self->{seq_lengths_ref};
  my $fileindex = $self->{fileindex};
  my $geneindex = $self->{geneindex};
  
  if ($debug) {
    print STDERR "\nnext gene parsing...\n";
    print STDERR "starting at fileindex, $fileindex\n";
  }
  
  # initialisation stuff
  
  my $filepos = tell $fh;
  if ($debug) {
    print STDERR "filepos before repositioning: $filepos\n";
  }
  
  seek ($fh, $fileindex, 0) or die "can't seek on file, " . $self->{file} . "!\n";
  $filepos = tell $fh;
  if ($debug) {
    print STDERR "filepos after repositioning: $filepos\n";
  }
  
  if (eof ($fh)) {
    
    if ($debug) {
      print STDERR "eof reached!\n";
    }
    
    # $fh->close;
    return (undef, undef);
  }
  
  my $new_gene     = 1;
  my $gene_feature = undef;
  my $gDNA_id      = "";
  my $cDNA_id      = "";
  my $gene_id      = "";
  my $algorithm    = "";
  my $cDNA_note    = "";
  my $subject_type;
  my @features     = ();
  my @mRNA_ranges  = ();
  
  # parsing ...
  
  while (my $line = $fh->getline) {
    
    # commented lines to filter ...
    if ($line =~ /^#/) {
	$fileindex++;
	next;
      }
    
    my $feature = $line;
    
    if ($debug) {
      print STDERR "\nfeature line before feature matching: $feature\n";
    }
    
    my $feature_without_plus_minus = $feature;
    # can't make the parsing working without plus or minus substituted !!
    $feature_without_plus_minus =~ s/[+-]/0/;
    
    if ($feature =~ /gene\t/) {
      
      if ($debug) {
	print STDERR "got a gene feature line\n";
      }
      
      if ($new_gene == 1) {
	
	# ReInitialisation
	$new_gene = 0;
	my @tmp = ();
	@mRNA_ranges = @tmp;
	
	# exonerate 0.6.7 and further
	
	$feature =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t\.\s(.*)/;
	
	$gDNA_id         = $1;
	$algorithm       = $2;
	my $feature_type = "gene";
	my $start        = $4;
	my $end          = $5;
	my $score        = $6;
	my $strand       = 0;
	
	if ($7 =~ /\+/) {
	  $strand = 1;
	}
	elsif ($7 =~ /-/) {
	  $strand = -1;
	}
	
	# Rest
	
	$feature_without_plus_minus =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t\.\s(.*)/;
	
	my $rest = $8;
	
	if ($debug) {
	  print STDERR "\tscore: $score\n";
	  print STDERR "\tstrand: $strand\n";
	  print STDERR "\trest: $rest\n";
	}
	
	my @fields = split (" ; ", $rest);
	# Gene ID
	# don't rely on exonerate output,
	# don't report unique IDs !!!
	$gene_id = "Gene." . $geneindex;
	$geneindex++;
	$self->{geneindex} = $geneindex;
	
	# cDNA Sequence ID
	$cDNA_id = $fields[1];
	$cDNA_id =~ s/sequence\s//;
	
	if ($algorithm =~ /est2genome/) {
	  $subject_type = "ESTs";
	}
	elsif ($algorithm =~ /protein2genome/) {
	  $subject_type = "protein";
	}

	if ($debug) {
	  print STDERR "\tgDNA name: $gDNA_id\n";
	  print STDERR "\tcDNA name: $cDNA_id\n";
	  print STDERR "\tgene id: $gene_id\n";
	}
	
	my $range = Bio::PSU::Range->new (
					  -start  => $start,
					  -end    => $end,
					  -strand => $strand
					 ); 
	$gene_feature = Bio::PSU::Feature->new (
						-key    => $feature_type,
						-ranges => [$range],
					       );
	$gene_feature->qadd ('gene', $gene_id);
	$gene_feature->qadd ('score', $score);
	my $note = "gene feature predicted by $algorithm";
	$gene_feature->qadd ('note', $note);
	$gene_feature->qadd ('method', $algorithm);
	
      }
      elsif ($new_gene == 0) {
	$new_gene = 1;

	# parsing done for this gene
	# don't move the file cursor because need to start the parsing for the next gene prediction on this current line !!!
	
	if ($debug) {
	  print STDERR "parsing done for this gene!\n";
	}
	
	return ($gDNA_id, \@features);
      }
    }
    
    elsif ($feature =~ /exon\t/) {
      
      if ($debug) {
	print STDERR "got an exon feature line\n";
      }
      
      # exonerate 0.6.7 and further
      
      $feature =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\+|-)\t\.\s(.+)/;
      
      # gDNA processed when parsing the gene feature line
      # $algorithm too
      my $feature_type = "exon";
      my $start        = $4;
      my $end          = $5;
      # no score
      my $strand       = 0;
      
      if ($7 =~ /\+/) {
	$strand = 1;
      }
      elsif ($7 =~ /-/) {
	$strand = -1;
      }
      
      # Rest => in a note qualifier of the exon feature
      
      $feature_without_plus_minus =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t\.\s(.*)/;
      my $exon_note = $8;
      
      my $range = Bio::PSU::Range->new(
				       -start       => $start,
				       -end         => $end,
				       -strand      => $strand
				      );
      my $feature = Bio::PSU::Feature->new (
					    -key    => $feature_type,
					    -ranges => [$range],
					   );
      
      $feature->qadd ('note', $exon_note);
      $feature->qadd ('gene', $gene_id);
      my $note = "exon feature predicted by $algorithm";
      $feature->qadd ('note', $note);
      
      # Adding just mRNA and gene Features !!
      # Exon Features commented
      # push (@features, $feature);
      
      # Add exon range to the mRNA_location array
      push (@mRNA_ranges, $range);
      
    }
    
    elsif ($feature =~ /splice\d\t/) {
      
      if ($debug) {
	print STDERR "got a splice site feature line\n";
      }
      
      # exonerate 0.6.7 and further
      
      $feature =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\+|-)\t\.\s(.+)/;
      
      # gDNA processed when parsing the gene feature line
      # $algorithm too
      my $feature_type = "splicesite";
      my $splicesite_type = $3;
      my $start        = $4;
      my $end          = $5;
      # no score
      my $strand       = 0;
      
      if ($7 =~ /\+/) {
	$strand = 1;
      }
      elsif ($7 =~ /-/) {
	$strand = -1;
      }
      
      # Rest
      
      $feature_without_plus_minus =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t\.\s(.*)/;

      my $rest = $8;
      
      if ($debug) {
	print STDERR "\trest: $rest\n";
      }
      
      my @fields = split (" ; ", $rest);
      # Gene ID
      # Splice Site information
      my $splicesite = $fields[1];
      $splicesite =~ s/\"//g;
      
      if ($debug) {
	print STDERR "\tgot splicesite note: $splicesite\n";
      }
      
      my $range = Bio::PSU::Range->new(
				       -start       => $start,
				       -end         => $end,
				       -strand      => $strand
				      );
      my $feature = Bio::PSU::Feature->new (
					    -key    => $feature_type,
					    -ranges => [$range],
					   );
      
      # need to add start and end of EST match coordinates somehow
      # from the similarity feature line
      
      my $note = "$splicesite";
      if ($splicesite_type =~ /splice3/) {
	$note = "3prime " . $note;
      }
      elsif ($splicesite_type =~ /splice5/) {
	$note = "5prime " . $note;
      }
      else {
	print STDERR "ERROR - can't figure out whether the splice site is 3prime or 5prime!!\n";
	print STDERR "Splice Site type: $splicesite_type\n";
      }
      $feature->qadd ('note', $note);
      # $feature->qadd ('algorithm', $algorithm);
      $note = "splice site feature predicted by $algorithm";
      $feature->qadd ('note', $note);
      $feature->qadd ('gene', $gene_id);

      # SS Feature Commented
      # push (@features, $feature);
      
      $fileindex++;
      
    }
    
    elsif ($feature =~ /similarity\t/) {
      
      if ($debug) {
	print STDERR "got a similarity feature line\n";
      }
      
      # exonerate 0.6.7 and further
      
      $feature =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\+|-)\t\.\s(.+)/;
      
      # gDNA processed when parsing the gene feature line
      # $algorithm too
      my $feature_type = "similarity";
      my $start        = $4;
      my $end          = $5;
      # score is the same one than the the gene score
      my $score        = $6;
      my $strand       = 0;
      
      if ($7 =~ /\+/) {
	$strand = 1;
      }
      elsif ($7 =~ /-/) {
	$strand = -1;
      }
      
      # Rest
      
      $feature_without_plus_minus =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t\.\s(.*)/;

      my $rest = $8;
      
      if ($debug) {
	print STDERR "\trest: $rest\n";
      }
      
      my @fields = split (" ; ", $rest);
      # Alignment ID
      $fields[0] =~ /(\d+)/;
      my $alignment_id = $1;
      # Alignment information
      my $i = 1;
      my @aligns = ();
      while ($i < @fields) {
	my $align_field = $fields[$i];
	$align_field =~ /align\s(\d+)\s(\d+)\s(\d+)/i;
	my $cDNA_start = $2;
	my $cDNA_end   = $3;
	
	# start and end might be reversed !!
	# keep them in their given order !!
	#if ($cDNA_start > $cDNA_end) {
	#  my $tmp     = $cDNA_start;
	#  $cDNA_start = $cDNA_end;
	#  $cDNA_end   = $tmp;
	#}
	
	my $align = "$cDNA_start..$cDNA_end";
	
	push (@aligns, $align);
	
	$i++;
      }
      
      my $range = Bio::PSU::Range->new(
				       -start       => $start,
				       -end         => $end,
				       -strand      => $strand
				      );
      my $feature = Bio::PSU::Feature->new (
					    -key    => $feature_type,
					    -ranges => [$range],
					   );
      $feature->qadd ('gene', $gene_id);	
      $cDNA_note = "Matched to $sequence_type sequence $cDNA_id at ";
      $cDNA_note .= "(" . $aligns[0];
      $i = 1;
      while ($i < @aligns) {
	$cDNA_note .= ", " . $aligns[$i];
	$i++;
      }
      $cDNA_note .= ")";
      $feature->qadd ('note', $cDNA_note);
      $gene_feature->qadd ('note', $cDNA_note);
      $gene_feature->qadd ("$subject_type", $cDNA_id);
      my $note = "similarity feature predicted by $algorithm";
      $feature->qadd ('note', $note);
      
      # don't need the similarity feature which contains the same information than the gene feature
      # push (@features, $feature);

      if (($algorithm =~ /est2genome/) && $add_gene_feature) {
	push (@features, $gene_feature);
      }
      
      # position the file fileindex on the similarity line, this way it always starts the parsing on the next feature, ie on the next gene (if any)
      
      if ($debug) {
	print STDERR "filepos: " . tell ($fh) . "\n";
      }
      $self->{fileindex} = tell ($fh);

      if ($debug) {
        print STDERR "adding mRNA feature in while similarity parsing...\n";
      }

      if (@mRNA_ranges > 0) {
	
	  # Create a mRNA Feature
	
	  if (not (defined $feature_key) && ($algorithm =~ /est2genome/)) {
	    $feature_key = "mRNA";
	  }
	  if (not defined $feature_key) {
	    die "ERROR, feature key is not defined, if you are using protein2genome model, you should give 'CDS' or 'BLASTCDS'!\n";
	  }

	  # process temparary mRNA ranges in case of protein2genome alignments
	  # because of some bugs in exonerate until they get fixed !!!

	  if (not defined $seq_length) {
	    $seq_length = $seq_lengths->{$gDNA_id};
	  }

	  if ($algorithm =~ /protein2genome/) {
	    @mRNA_ranges = fixRanges ($seq_length, @mRNA_ranges);
	  }
	  
	  my $mRNA_feature = Bio::PSU::Feature->new (
						     -key    => $feature_key,
	  					     -ranges => [@mRNA_ranges],
						    );
	  $mRNA_feature->qadd ('gene', $gene_id);
	  my $note = "$feature_key feature predicted by $algorithm";
	  $mRNA_feature->qadd ('note', $note);
          $mRNA_feature->qadd ('method', $algorithm);
	  $mRNA_feature->qadd ('note', $cDNA_note);
          $mRNA_feature->qadd ("$subject_type", $cDNA_id);
          $mRNA_feature->qadd ('score', $score);
          
          # if it's a protein, get the description using getz
	  # and it is BLASTCDS mode

          if (($subject_type eq "protein") && ($feature_key eq "BLASTCDS")) {
            my $product = qx(getz '[swall-ALLTEXT:$cDNA_id]' -f "des");
            if (defined $product) {
              $product =~ /DE\s+(.+)\./;
              $product = $1;
              $mRNA_feature->qadd ('product', "putative $product");
            }
          }
          
	  if ($feature_key eq "BLASTCDS") {
	    # Used as protein similarity evidence
	    # the note below id required by GFMerge
	    
	    $mRNA_feature->qadd ('note', "Exonerate similarity (percentage match=unknown) to $cDNA_id , unknown..unknown, expect unknown, score unknown");
	  }

	  push (@features, $mRNA_feature);
      }
    }
    else {
      if ($debug) {
	print STDERR "Feature line: $feature\n";
	print STDERR "the feature key is not known or line can not be parsed!!\n";
      }
    }

    $fileindex++;
    
  } # eof !!
  
  if ($debug) {
    print STDERR "eof, returning the last gene prediction...\n";
  }

  # $self->{fileindex} = $fileindex;
  # $self->{fh} = $fh;
  # $fh->close;
  return ($gDNA_id, \@features);
  
}

sub close_fh {
  my $self = shift;
  my $fh = $self->{fh};
  $fh->close;
}

sub fixRanges {
  my ($seq_length, @mRNA_ranges) = @_;
  my @new_mRNA_ranges = ();
  my $strand = $mRNA_ranges[0]->strand;

  if ($strand eq -1) {
    # map coordinates on the forward strand...
    my $i = 0;
    foreach my $range (@mRNA_ranges) {
      my $end   = $seq_length - $range->start + 1;
      my $start = $seq_length - $range->end + 1;

      my $new_range = Bio::PSU::Range->new (
					    -start  => $start,
					    -end    => $end,
					    -strand => $strand,
					   );
      push (@new_mRNA_ranges, $new_range);
      $i++;
    }
  }
  else {
    @new_mRNA_ranges = @mRNA_ranges;
    # forward strand, just shift downstream the start of the second exon
    # or all internal + ending exons ????
    if (@mRNA_ranges > 1) {
      my $range = $mRNA_ranges[1];
      my $new_range = Bio::PSU::Range->new (
					    -start  => $range->start,
					    -end    => $range->end,
					    -strand => $strand,
					   );
      splice (@new_mRNA_ranges, 1, 1, $new_range);
    }
  }

  return @new_mRNA_ranges;
}

1;
