package NHGRI::MUMmer::AlignSet;
############################################################

=head1 NAME

AlignSet.pm - A Perl module to parse a MUMmer delta file and 
          return an AlignSet object, which contains details
          about the MUMmer alignments parsed from the delta file.

=head1 DESCRIPTION

  MUMmer's delta file format is rather obtuse, so this module
  will parse the file and allow perl programmers to access
  the alignments in a more intuitive fashion.

=head1 DATE

 October, 2016

=head1 AUTHORS

 Nancy Fisher Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=cut

############################################################
use strict;
use warnings;

use Carp;
use GTB::File qw(Open);
use GTB::FASTA;

our $VERSION  = '0.01';

###########################################################

=over 4

=item B<new()>

  This method creates a AlignSet object from a delta file path.

  Input:  -delta_file - the path to a ".delta" file created 
          by MUMmer.
          -ignore_unequal - ignore mismatches in the number of
          reference and query bases in an alignment (default is
          to die).
          -storerefentrypairs - store a hash of lists of entry pairs
          for each reference entry
          -storequeryentrypairs - store a hash of lists of entry pairs
          for each query entry
          -extend_exact - option to extend alignments as far as is
          possible with exactly matching sequence
          -reference_hashref - optional reference to hash of reference
          sequences (will be used instead of reading fasta file to
          extend alignments with -extend_exact option)
          -query_hashref - optional reference to hash of query
          sequences (will be used instead of reading fasta file to
          extend alignments with -extend_exact option)
          -reference_file - option to set alternate reference file
          (rather than use reference_file present in delta file)
          -query_file - option to set alternate query file
  Output: New AlignSet object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $delta_file = $params{-delta_file};
    my $ignore_unequal = $params{-ignore_unequal};
    my $storerefentrypairs = $params{-storerefentrypairs};
    my $storequeryentrypairs = $params{-storequeryentrypairs};
    my $extend_exact = $params{-extend_exact};
    my $reference_hashref = $params{-reference_hashref};
    my $query_hashref = $params{-query_hashref};
    my $reference_file = $params{-reference_file};
    my $query_file = $params{-query_file};

    my $self = { delta_file => $delta_file,
                 ignore_unequal => $ignore_unequal,
                 storerefentrypairs => $storerefentrypairs,
                 storequeryentrypairs => $storequeryentrypairs,
                 extend_exact => $extend_exact,
                 reference_hashref => $reference_hashref,
                 query_hashref => $query_hashref,
                 reference_file => $reference_file,
                 query_file => $query_file,
                  };

    bless $self, $class;

    # for now, reads in all alignments at the outset:

    $self->_parse_delta_file();

    return $self;
}

###########################################################

=item B<_parse_delta_file()>

  This method reads through the delta file, and 
  populates the various components of the AlignSet
  object.

  The following fields of the object hash are 
  populated:

  reference_file, query_file - files aligned in this delta file
      (only used if not specified in constructor)
  entry_pairs - reference to an array of hash references, each
      of which contains the matches between an entry from the 
      reference file and an entry from the query file.  The fields
      of an entry pair's hash are:

      ref_entry, query_entry - entry names from the two files
      ref_length, query_length - lengths of the two entries
      aligns - reference to an array of hash reference for each
          align, with keys: 
            ref_start, ref_end, query_start, query_end, mismatches, 
            nonposmatches, nonalphas, cigar_string
        

  Input: None.
  Output: None.

=cut

###########################################################
sub _parse_delta_file {
    my $self  = shift;

    $self->{delta_fh} = Open($self->{delta_file});
    my $fh = $self->{delta_fh};

    my ($ref_entry, $query_entry);
    while (<$fh>) {
        if (/^(\S+)\s(\S+)$/) { # first line contains files from run
            $self->{reference_file} = $self->{reference_file} || $1;
            $self->{query_file} = $self->{query_file} || $2;
            if ($self->{extend_exact}) {
                if (!($self->{reference_hashref}) || !($self->{query_hashref})) { # need to read seqs from FASTA
                    $self->{ref_fasta_db} = GTB::FASTA->new($self->{reference_file});
                    $self->{query_fasta_db} = GTB::FASTA->new($self->{query_file});
                    $self->{ref_fasta_db}->reindex();
                    $self->{query_fasta_db}->reindex();
                }
            }
        }
        elsif (/^\>(\S+)\s(\S+)\s(\d+)\s(\d+)$/) {
            ($ref_entry, $query_entry) = ($1, $2);
            push @{$self->{entry_pairs}}, 
                 {ref_entry => $ref_entry, query_entry => $query_entry, ref_length => $3, query_length => $4};
        }
        elsif (/^(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my ($ref_start, $ref_end, $query_start, $query_end, $mismatches, $nonposmatches, $nonalphas) = 
                ($1, $2, $3, $4, $5, $6, $7);
            my $cigar_string = $self->_read_cigar_string($ref_start, $ref_end, $query_start, $query_end);
            if ($self->{extend_exact}) {
                my ($new_ref_start, $new_ref_end, $new_query_start, $new_query_end, $new_cigar) = 
                    $self->_extend_exact($ref_entry, $ref_start, $ref_end, $query_entry, $query_start, $query_end, $cigar_string);
                $cigar_string = $new_cigar;
                $ref_start = $new_ref_start;
                $ref_end = $new_ref_end;
                $query_start = $new_query_start;
                $query_end = $new_query_end;
            }
            my $comp = ($query_start <= $query_end) ? 0 : 1;
            push @{$self->{entry_pairs}->[$#{$self->{entry_pairs}}]->{aligns}}, 
                {ref_start => $ref_start, ref_end => $ref_end, query_start => $query_start,
                 query_end => $query_end, mismatches => $mismatches, nonposmatches => $nonposmatches,
                 nonalphas => $nonalphas, cigar_string => $cigar_string, ref_entry => $ref_entry,
                 query_entry => $query_entry, comp => $comp };
            if (!$ref_start || !$ref_end) {
                die "No ref start or ref end for entry pair $self->{entry_pairs}->[$#{$self->{entry_pairs}}]->{ref_entry} and $self->{entry_pairs}->[$#{$self->{entry_pairs}}]->{query_entry}\n";
            }
        }
    }

    if ($self->{storerefentrypairs} || $self->{storequeryentrypairs}) {
        foreach my $rh_entrypair (@{$self->{entry_pairs}}) {
            if ($self->{storerefentrypairs}) {
                push @{$self->{refentrypairs}->{$rh_entrypair->{ref_entry}}}, $rh_entrypair;
            }
            if ($self->{storequeryentrypairs}) {
                push @{$self->{queryentrypairs}->{$rh_entrypair->{query_entry}}}, $rh_entrypair;
            }

        }
    }
    close $self->{delta_fh};

    return;

} ## end _parse_delta_file

###########################################################

=item B<_read_cigar_string()>

  This method reads the lines of the delta file corresponding to a single 
  alignment, and returns a cigar string for that alignment (advancing the
  file to the end of the alignment's section in the file).

  Input: None, but file handle file position must be at beginning of alignment.
  Output: Cigar string (scalar string).

=cut

###########################################################
sub _read_cigar_string {
    my $self  = shift;
    my $ref_start = shift;
    my $ref_end = shift;
    my $query_start = shift;
    my $query_end = shift;

    my $comp = ($query_start <= $query_end) ? 0 : 1;
    my $cigar_string = '';
    my ($current_dels, $current_ins) = (0, 0);
    my $fh = $self->{delta_fh};
    my ($current_ref, $current_query) = ($comp) ? ($ref_start - 1, $query_start + 1) :
                                ($ref_start - 1, $query_start - 1); # will advance to each position as line is read
    while (<$fh>) {
        chomp;
        last if ($_ eq '0'); # last line

        my $dist = $_;
        if ($dist > 0) {
            my $matches = $dist - 1;
            if ($matches) { # print out indels first
                if ($current_dels) {
                    $cigar_string .= $current_dels."D";
                    $current_ref += $current_dels;
                    $current_dels = 0;
                }
                if ($current_ins) {
                    $cigar_string .= $current_ins."I";
                    $current_query = ($comp) ? $current_query - $current_ins : $current_query + $current_ins;
                    $current_ins = 0;
                }
                $cigar_string .= $matches."M";
                $current_ref += $matches;
                $current_query = ($comp) ? $current_query - $matches : $current_query + $matches;
            }
            $current_dels++;
        }
        elsif ($dist < 0) {
            my $matches = -1*$dist - 1;
            if ($matches) { # print out indels first
                if ($current_dels) {
                    $cigar_string .= $current_dels."D";
                    $current_ref += $current_dels;
                    $current_dels = 0;
                }
                if ($current_ins) {
                    $cigar_string .= $current_ins."I";
                    $current_query = ($comp) ? $current_query - $current_ins : $current_query + $current_ins;
                    $current_ins = 0;
                }
                $cigar_string .= $matches."M";
                $current_ref += $matches;
                $current_query = ($comp) ? $current_query - $matches : $current_query + $matches;
            }
            $current_ins++;
        }
    }

    if ($current_dels) {
        $cigar_string .= $current_dels."D";
        $current_ref += $current_dels;
        $current_dels = 0;
    }
    if ($current_ins) {
        $cigar_string .= $current_ins."I";
        $current_query = ($comp) ? $current_query - $current_ins : $current_query + $current_ins;
        $current_ins = 0;
    }
    # add matches to the end:
    if (((!$comp) && ($query_end - $current_query != $ref_end - $current_ref)) ||
        (($comp) && ($current_query - $query_end != $ref_end - $current_ref))) {
        if (!$self->{ignore_unequal}) {
            die "Unequal number of remaining bases in delta alignment in $self->{delta_file}!\n";
        }
    }
    else {
        my $remaining_matches = $ref_end - $current_ref;
        $cigar_string .= $remaining_matches."M";
        return $cigar_string;
    }

} ## end _read_cigar_string

###########################################################

=item B<_extend_exact()>

  This method retrieves bases from the reference and query fasta files
  to determine whether the alignment in the delta file can be extended
  exactly.

  Input: AlignSet object, then ref entry, ref start, ref end, query entry, 
         query start, and query end
  Output: New ref start, ref end, query start, and query end

=cut

###########################################################
sub _extend_exact {
    my $self  = shift;
    my $ref_entry = shift;
    my $ref_start = shift;
    my $ref_end = shift;
    my $query_entry = shift;
    my $query_start = shift;
    my $query_end = shift;
    my $cigar = shift;
    
    my ($new_ref_start, $new_ref_end, $new_query_start, $new_query_end, $new_cigar) = 
       ($ref_start, $ref_end, $query_start, $query_end, $cigar);

    my $rh_refseqs = $self->{reference_hashref};
    my $rh_queryseqs = $self->{query_hashref};

    my $ref_fasta_db = $self->{ref_fasta_db};
    my $query_fasta_db = $self->{query_fasta_db};

    my $ref_length = ($rh_refseqs) ? length($rh_refseqs->{$ref_entry}) : $ref_fasta_db->len($ref_entry);
    my $query_length = ($rh_queryseqs) ? length($rh_queryseqs->{$query_entry}) : $query_fasta_db->len($query_entry);

    my $comp = ($query_start > $query_end) ? 1 : 0;

    # extend to right:
    my $ref_next_pos = $new_ref_end + 1;
    my $query_next_pos = ($comp) ? $new_query_end - 1 : $new_query_end + 1;

    if ($ref_next_pos <= $ref_length && $query_next_pos >= 1 && $query_next_pos <= $query_length) {
        my $ref_next_base = ($rh_refseqs) ? uc(substr($rh_refseqs->{$ref_entry}, $ref_next_pos - 1, 1)) : 
                                            uc($ref_fasta_db->seq("$ref_entry:$ref_next_pos-$ref_next_pos"));
        my $query_next_base = ($rh_queryseqs) ? uc(substr($rh_queryseqs->{$query_entry}, $query_next_pos - 1, 1)) :
                                                uc($query_fasta_db->seq("$query_entry:$query_next_pos-$query_next_pos"));
        $query_next_base =~ tr/ATGC/TACG/ if ($comp);
    
        my $right_extend_bases = ($new_cigar =~ s/(\d+)M$//) ? $1 : 0;
        while ($ref_next_base eq $query_next_base) {
            $new_ref_end = $ref_next_pos;
            $new_query_end = $query_next_pos;
            $right_extend_bases++;
    
            # retrieve next bases:
            $ref_next_pos = $new_ref_end + 1;
            $query_next_pos = ($comp) ? $new_query_end - 1 : $new_query_end + 1;
            
            last if ($ref_next_pos > $ref_length || $query_next_pos < 1 || $query_next_pos > $query_length);
            $ref_next_base = ($rh_refseqs) ? uc(substr($rh_refseqs->{$ref_entry}, $ref_next_pos - 1, 1)) : 
                             uc($ref_fasta_db->seq("$ref_entry:$ref_next_pos-$ref_next_pos"));
            $query_next_base = ($rh_queryseqs) ? uc(substr($rh_queryseqs->{$query_entry}, $query_next_pos - 1, 1)) :
                               uc($query_fasta_db->seq("$query_entry:$query_next_pos-$query_next_pos"));
            $query_next_base =~ tr/ATGC/TACG/ if ($comp);
        }
        $new_cigar = $new_cigar.$right_extend_bases.'M' if ($right_extend_bases);
        #print STDERR "Rightside REF $ref_next_base does not match QUERY $query_next_base\n";
    }

    # extend to left:
    $ref_next_pos = $new_ref_start - 1;
    $query_next_pos = ($comp) ? $new_query_start + 1 : $new_query_start - 1;

    if ($ref_next_pos >= 1 && $query_next_pos >= 1 && $query_next_pos <= $query_length) {
        my $ref_next_base = ($rh_refseqs) ? uc(substr($rh_refseqs->{$ref_entry}, $ref_next_pos - 1, 1)) : 
                            uc($ref_fasta_db->seq("$ref_entry:$ref_next_pos-$ref_next_pos"));
        my $query_next_base = ($rh_queryseqs) ? uc(substr($rh_queryseqs->{$query_entry}, $query_next_pos - 1, 1)) :
                              uc($query_fasta_db->seq("$query_entry:$query_next_pos-$query_next_pos"));
        $query_next_base =~ tr/ATGC/TACG/ if ($comp);
    
        my $left_extend_bases = ($new_cigar =~ s/^(\d+)M//) ? $1 : 0;
        while ($ref_next_base eq $query_next_base) {
            $new_ref_start = $ref_next_pos;
            $new_query_start = $query_next_pos;
            $left_extend_bases++;
    
            # retrieve next bases:
            $ref_next_pos = $new_ref_start - 1;
            $query_next_pos = ($comp) ? $new_query_start + 1 : $new_query_start - 1;

            last if ($ref_next_pos < 1 || $query_next_pos < 1 || $query_next_pos > $query_length);
            $ref_next_base = ($rh_refseqs) ? uc(substr($rh_refseqs->{$ref_entry}, $ref_next_pos - 1, 1)) : 
                             uc($ref_fasta_db->seq("$ref_entry:$ref_next_pos-$ref_next_pos"));
            $query_next_base = ($rh_queryseqs) ? uc(substr($rh_queryseqs->{$query_entry}, $query_next_pos - 1, 1)) :
                             uc($query_fasta_db->seq("$query_entry:$query_next_pos-$query_next_pos"));
            $query_next_base =~ tr/ATGC/TACG/ if ($comp);
        }
        $new_cigar = $left_extend_bases.'M'.$new_cigar if ($left_extend_bases);
        #print STDERR "Leftside REF $ref_next_base does not match QUERY $query_next_base\n";
    }

    if ($new_cigar ne $cigar || $new_ref_start != $ref_start || $new_ref_end != $ref_end) {
        print STDERR "Ref $ref_entry:$ref_start-$ref_end extended to $new_ref_start-$new_ref_end\n";
        print STDERR "Query $query_entry:$query_start-$query_end extended to $new_query_start-$new_query_end\n";
    }
    else {
        #print STDERR "No extension! Ref $ref_entry:$ref_start-$ref_end (length $ref_length) Query $query_entry:$query_start-$query_end $cigar (length $query_length)\n";
    }

    return ($new_ref_start, $new_ref_end, $new_query_start, $new_query_end, $new_cigar);

} ## end _extend_exact

###########################################################

=item B<find_ref_coords_from_query_coord()>

  This method returns a list of reference coordinates corresponding to
  query coordinate specified, searching each of the alignments in the 
  object.

  Input: Query, query coordinate
  Output: Populates each alignment for this query with a "ref_matches"
      field and an entry for the query position giving the value of 
      the reference coordinate that corresponds to the query position.
      Returns 1 if the matching ref coordinate was found, 0 otherwise.

=cut

###########################################################
sub find_ref_coords_from_query_coord {
    my $self  = shift;
    my $query_pos = shift;
    my $query = shift;
    my $rh_params = shift; # not used yet

    my $match_found = 0;
    foreach my $rh_pairentry (@{$self->{entry_pairs}}) {
        next if ($query && $rh_pairentry->{query_entry} ne $query);

        foreach my $rh_align (@{$rh_pairentry->{aligns}}) {
            my $align_match_found = 0;
            my $ref_start = $rh_align->{ref_start};
            my $ref_end = $rh_align->{ref_end};
            my $query_start = $rh_align->{query_start};
            my $query_end = $rh_align->{query_end};
            my $comp = ($query_start <= $query_end) ? 0 : 1;
            #print STDERR "Searching for $query_pos in align $query:$query_start-$query_end, ref $ref_start-$ref_end\n";
            if (!$comp && ($query_pos < $query_start || $query_pos > $query_end)) {
                #print "$query_pos is not between $query_start and $query_end!\n";
                next;
            }
            elsif ($comp && ($query_pos > $query_start || $query_pos < $query_end)) {
                #print "$query_pos is not between $query_start and $query_end!\n";
                next;
            }

            my $current_ref = $ref_start - 1;
            my $current_query = ($comp) ? $query_start + 1 : $query_start - 1;

            my $cigar_string = $rh_align->{cigar_string};
            #print "Cigar string is $cigar_string!\n";

            while ($cigar_string) {
                my ($bases, $op);
                if ($cigar_string =~ s/^(\d+)([MDI])//) {
                    ($bases, $op) = ($1, $2);
                    if ($op eq 'M') {
                        if ($comp && ($current_query - $bases <= $query_pos)) {
                            $rh_align->{ref_matches}->{$query_pos} = $current_ref - ($query_pos - $current_query);
                            $match_found = 1;
                            $align_match_found = 1;
                            last;
                        }
                        elsif (!$comp && ($current_query + $bases >= $query_pos)) {
                            $rh_align->{ref_matches}->{$query_pos} = $current_ref + ($query_pos - $current_query);
                            $match_found = 1;
                            $align_match_found = 1;
                            last;
                        }
                        else {
                            $current_ref += $bases;
                            $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                        }
                    }
                    elsif ($op eq 'D') { # deleted from reference--advance ref coord
                        $current_ref += $bases;
                    }
                    elsif ($op eq 'I') { # inserted to reference--advance query coord
                        $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
    
                        if (($comp && $current_query <= $query_pos) || (!$comp && $current_query >= $query_pos)) {
                            $rh_align->{ref_matches}->{$query_pos} = $current_ref;
                            $match_found = 1;
                            $align_match_found = 1;
                            #print STDERR "Setting ref_match within insertion to $current_ref\n";
                            last;
                        }
                    }
                }
                else {
                    die "Cigar string $cigar_string is of the wrong form!\n";
                }
            }
            if (!$align_match_found) { # check to be sure we got to the end
                print STDERR "Reached query $current_query (end is $query_end), ref $current_ref (end is $ref_end) without finding match\n";
            }
            else {
                #print STDERR "Found match: $rh_align->{ref_matches}->{$query_pos}\n";
            }
        }
    }
    return $match_found;
}

###########################################################

=item B<find_query_coords_from_ref_coord()>

  This method returns a list of query coordinates corresponding to
  reference coordinate specified, searching each of the alignments in the 
  object.

  Input: Reference, reference coordinate
  Output: Populates each alignment for this reference with a "query_matches"
      field and an entry for the reference position giving the value of 
      the query coordinate that corresponds to the reference position.

=cut

###########################################################
sub find_query_coords_from_ref_coord {
    my $self  = shift;
    my $ref_pos = shift;
    my $ref = shift;
    my $rh_params = shift; # not used yet

    #print "In find_query_coords!\n";
    foreach my $rh_pairentry (@{$self->{entry_pairs}}) {
        next if ($ref && $rh_pairentry->{ref_entry} ne $ref);

        foreach my $rh_align (@{$rh_pairentry->{aligns}}) {
            my $query_start = $rh_align->{query_start};
            my $ref_start = $rh_align->{ref_start};
            my $query_end = $rh_align->{query_end};
            my $ref_end = $rh_align->{ref_end};
            my $comp = ($query_start < $query_end) ? 0 : 1;
            if (($ref_pos < $ref_start) || ($ref_pos > $ref_end)) {
                #print "$ref_pos is not between $ref_start and $ref_end!\n";
                next;
            }

            my $current_ref = $ref_start - 1;
            my $current_query = ($comp) ? $query_start + 1 : $query_start - 1;

            my $cigar_string = $rh_align->{cigar_string};

            while ($cigar_string) {
                my ($bases, $op);
                if ($cigar_string =~ s/^(\d+)([MDI])//) {
                    ($bases, $op) = ($1, $2);
                    if ($op eq 'M') {
                        if ($current_ref + $bases >= $ref_pos) {
                            $rh_align->{query_matches}->{$ref_pos} = ($comp) ? $current_query - ($ref_pos - $current_ref) : $current_query + ($ref_pos - $current_ref);
                            #print "Value is $rh_align->{query_matches}->{$ref_pos}\n";
                            last;
                        }
                        else {
                            $current_ref += $bases;
                            $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                        }
                    }
                    elsif ($op eq 'D') { # deleted from reference--advance ref coord
                        $current_ref += $bases;
                        if ($current_ref >= $ref_pos) {
                            $rh_align->{query_matches}->{$ref_pos} = $current_query;
                            last;
                        }
                    }
                    elsif ($op eq 'I') { # inserted to query--advance query coord
                        $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
    
                    }
                }
                else {
                    die "Cigar string $cigar_string is of the wrong form!\n";
                }
            }
        }
    }
}

###########################################################

=item B<find_gap_query_coords()>

  This method returns the query contig and the query coordinates of a specified
  reference entry and alignment endpoints for the aligns surrounding a GAP found
  by the MUMmer show-diff program.

  Input: AlignSet object, reference entry, reference coordinate of end of first
         alignment, reference coordinate of start of second alignment.
  Output: array containing the query contig aligned to the reference, the
         query coordinates corresponding to the first and second reference 
         alignment endpoints, and 1 or 0 depending on whether the query contig
         is complemented or uncomplemented with respect to the reference, 
         respectively.

=cut

###########################################################
sub find_gap_query_coords {
    my $self  = shift;
    my $ref_entry = shift;
    my $ref1 = shift;
    my $ref2 = shift;

    $self->find_query_coords_from_ref_coord($ref1, $ref_entry);
    $self->find_query_coords_from_ref_coord($ref2, $ref_entry);

    my %found_query_matches = ();
    foreach my $rh_entry (@{$self->{entry_pairs}}) {
        my ($querycontig, $query1, $query2, $comp);
        next if ($rh_entry->{ref_entry} ne $ref_entry);
        foreach my $rh_align (@{$rh_entry->{aligns}}) {
            if (defined($rh_align->{query_matches}->{$ref1})) { 
                my ($rstart, $rend, $qstart, $qend) = ($rh_align->{ref_start},
                    $rh_align->{ref_end}, $rh_align->{query_start}, $rh_align->{query_end});
                if ($rend == $ref1) { # found query1
                    $querycontig = $rh_entry->{query_entry};
                    $query1 = $qend;
                    $comp = ($qend >= $qstart) ? 0 : 1;
                    $found_query_matches{$querycontig} = {'query1' => $query1, 'comp1' => $comp};
                    $rh_align->{query_matches}->{$ref1} = 'NA' if (!$rh_align->{query_matches}->{$ref1});
                    if ($rh_align->{query_matches}->{$ref1} != $query1) {
                        my $wrong_match = $rh_align->{query_matches}->{$ref1};
                        die "Didn\'t find correct match $query1 for ref1 $ref1 in first align (found $wrong_match)!\n";
                    }
                }
            }
            if (defined($rh_align->{query_matches}->{$ref2})) { 
                my ($rstart, $rend, $qstart, $qend) = ($rh_align->{ref_start},
                    $rh_align->{ref_end}, $rh_align->{query_start}, $rh_align->{query_end});
                if ($rstart == $ref2) { # found query2
                    $querycontig = $rh_entry->{query_entry};
                    $query2 = $qstart;
                    $comp = ($qend >= $qstart) ? 0 : 1;
                    $found_query_matches{$querycontig} = {'query2' => $query2, 'comp2' => $comp};
                    $rh_align->{query_matches}->{$ref2} = 0 if (!$rh_align->{query_matches}->{$ref2});
                    if ($rh_align->{query_matches}->{$ref2} != $query2) {
                        my $wrong_match = $rh_align->{query_matches}->{$ref2};
                        print "Didn\'t find correct match $querycontig:$query2 for ref2 $ref2 in second align (found $wrong_match)!\n";
                    }
                }
            }
        }
    }
    # Check for gap spanning matches to $ref_entry:$ref1-$ref2
    foreach my $varcontig (keys %found_query_matches) {
        #print "Checking $varcontig for gap spanning matches to $ref_entry:$ref1-$ref2\n";
        my $query1 = $found_query_matches{$varcontig}->{'query1'};
        my $query2 = $found_query_matches{$varcontig}->{'query2'};
        my $comp1 = $found_query_matches{$varcontig}->{'comp1'};
        my $comp2 = $found_query_matches{$varcontig}->{'comp2'};
        if ($query1 && $query2 && $comp1 && $comp2 && ($comp1 == $comp2)) { # everything's here for this contig
            #print "Returning $varcontig:$query1-$query2\n";
            return ($varcontig, $query1, $query2, $comp1);
        }
    }

}

###########################################################

1;
__END__

=back
