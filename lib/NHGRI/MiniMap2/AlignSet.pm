package NHGRI::MiniMap2::AlignSet;
############################################################

=head1 NAME

AlignSet.pm - A Perl module to parse a MiniMap2 bam file and 
          return an AlignSet object, which contains details
          about the MiniMap2 alignments parsed from the bam file.

=head1 DESCRIPTION

  MiniMap2's bam file contains alignments which can be parsed
  and returned as a set of 

=head1 DATE

 October, 2018

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

  Input:  -bam_file - the path to a BAM file created by MiniMap2.
          -storerefentrypairs - store a hash of lists of entry pairs
          for each reference entry
          -storequeryentrypairs - store a hash of lists of entry pairs
          for each query entry
          -reference_file - reference FASTA file
          -query_file - query FASTA file
  Output: New AlignSet object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $bam_file = $params{-bam_file};
    my $storerefentrypairs = $params{-storerefentrypairs};
    my $storequeryentrypairs = $params{-storequeryentrypairs};
    my $reference_file = $params{-reference_file};
    my $query_file = $params{-query_file};
    my $verbose = $params{-verbose};

    my $self = { bam_file => $bam_file,
                 storerefentrypairs => $storerefentrypairs,
                 storequeryentrypairs => $storequeryentrypairs,
                 reference_file => $reference_file,
                 query_file => $query_file,
                 verbose => $verbose,

                  };

    bless $self, $class;

    if (!$self->{reference_file} || !$self->{query_file}) {
        die "Must specify reference file with -reference_file parameter and query file with -query_file parameter in $class constructor!\n";
    }

    $self->{reference_db} = GTB::FASTA->new($reference_file);
    $self->{query_db} = GTB::FASTA->new($query_file);

    # for now, reads in all alignments at the outset:

    $self->_parse_bam_file();

    return $self;
}

###########################################################

=item B<_parse_bam_file()>

  This method reads through the bam file, and populates the
  various components of the AlignSet object.

  The following fields of the object hash are 
  populated:

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
sub _parse_bam_file {
    my $self  = shift;

    my $reference_db = $self->{reference_db};
    my $query_db = $self->{query_db};

    my $command = "samtools view -F 256 $self->{bam_file} | ";

    $self->{bam_fh} = Open($command);
    my $fh = $self->{bam_fh};

    my ($ref_entry, $query_entry);
    my $readcount = 0;
    while (<$fh>) {
        chomp;
        my ($query_entry, $flag, $ref_entry, $refpos, $mapq, $cigar_string, $matename, $mpos, $isize, $seq, $qual, $opt) = split /\t/, $_;
        next if ($cigar_string eq '*'); # skip unaligned
        my $comp = ( $flag & 16 ) ? 1 : 0;
        my $ra_align_endpoints = $self->_get_ref_query_positions($ref_entry, $query_entry, $refpos, $comp, $cigar_string);

        my $ref_length = $reference_db->len($ref_entry);
        my $query_length = $query_db->len($query_entry);

        my ($entry_pair) = grep {$_->{ref_entry} eq $ref_entry && $_->{query_entry} eq $query_entry} @{$self->{entry_pairs}};
        if (!$entry_pair) {
            push @{$self->{entry_pairs}}, {ref_entry => $ref_entry, query_entry => $query_entry, ref_length => $ref_length, query_length => $query_length};
            $entry_pair = $self->{entry_pairs}->[$#{$self->{entry_pairs}}];
        }

        foreach my $ra_align (@{$ra_align_endpoints}) {
            my ($ref_start, $ref_end, $query_start, $query_end) = @{$ra_align};

            push @{$entry_pair->{aligns}}, 
                {ref_start => $ref_start, ref_end => $ref_end, query_start => $query_start,
                 query_end => $query_end, cigar_string => $cigar_string, ref_entry => $ref_entry,
                 # could eventually try to parse these values out of minimap2 alignments, but not now
                 # mismatches => $mismatches, nonposmatches => $nonposmatches, nonalphas => $nonalphas,
                 query_entry => $query_entry, comp => $comp };

            print "ALIGNMENT\t$ref_entry\t$ref_start\t$ref_end\t$query_entry\t$query_start\t$query_end\t$comp\n" if ($self->{verbose});
        }

        $readcount++;
        #last if ($readcount > 10000);
    }

    print "Value of storerefentrypairs is $self->{storerefentrypairs}\n";
    if ($self->{storerefentrypairs} || $self->{storequeryentrypairs}) {
        print "Making refentry hash!\n";
        foreach my $rh_entrypair (@{$self->{entry_pairs}}) {
            print "Entry pair!\n";
            if ($self->{storerefentrypairs}) {
                push @{$self->{refentrypairs}->{$rh_entrypair->{ref_entry}}}, $rh_entrypair;
                print "Entry chrom: $rh_entrypair->{ref_entry}\n";
            }
            if ($self->{storequeryentrypairs}) {
                push @{$self->{queryentrypairs}->{$rh_entrypair->{query_entry}}}, $rh_entrypair;
            }

        }
    }
    close $fh;

    return;

} ## end _parse_bam_file

###########################################################

=item B<_get_ref_query_positions()>

  This method determines endpoints of alignments in the
  query and reference and the query start position from
  parameters of a SAM-formatted alignment.

  Input: Object, ref entry, query entry, reference position, comp, cigar
  Output: Reference to a list of refs to arrays of (ref start, ref end,
          query start, query end)

=cut

###########################################################
sub _get_ref_query_positions {
    my $self  = shift;
    my $ref_entry = shift;
    my $query_entry = shift;
    my $ref_pos = shift;
    my $comp = shift;
    my $cigar_string = shift;

    my $reference_db = $self->{reference_db};
    my $query_db = $self->{query_db};

    my ($ref_start, $ref_end, $query_start, $query_end);

    $ref_start = $ref_pos;
    $query_start = ($comp) ? $query_db->len($query_entry) : 1;

    my @starts_ends = ();
    my ($current_ref, $current_query) = ($ref_start, $query_start);
    my $first_op = 1; # advance start for soft clipping if first operation
    while ($cigar_string) {
        my ($bases, $op);
        if ($cigar_string =~ s/^(\d+)([MDINSHP])//) {
            ($bases, $op) = ($1, $2);
            if ($op eq 'M') {
                $current_ref += $bases;
                $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
            }
            elsif ($op eq 'D') { # deleted from reference--advance ref coord
                $current_ref += $bases;
            }
            elsif ($op eq 'I') { # inserted to reference--advance query coord
                $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
            }
            elsif ($op eq 'S') { # soft clipping--advance query start coord if at beginning
                if ($first_op) {
                    $query_start = ($comp) ? $query_start - $bases : $query_start + $bases;
                    $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                }
            }
            elsif ($op eq 'H') { # hard clipping--advance query start coord if at beginning
                if ($first_op) {
                    $query_start = ($comp) ? $query_start - $bases : $query_start + $bases;
                    $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                }
            }
        }
        else {
            die "Cigar string $cigar_string is of the wrong form!\n";
        }
        $first_op = 0;
    }

    $ref_end = $current_ref - 1;
    $query_end = ($comp) ? $current_query + 1 : $current_query - 1;

    print "SAM $ref_entry/$query_entry $ref_pos comp $comp\n";
    print "YIELDS $ref_start-$ref_end $query_start-$query_end\n";

    push @starts_ends, [$ref_start, $ref_end, $query_start, $query_end];

    return [@starts_ends];

} # end _get_ref_query_positions

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
            if (!$comp && ($query_pos < $query_start || $query_pos > $query_end)) {
                next;
            }
            elsif ($comp && ($query_pos > $query_start || $query_pos < $query_end)) {
                next;
            }

            my $current_ref = $ref_start - 1;
            my $current_query = ($comp) ? $query_start + 1 : $query_start - 1;

            my $cigar_string = $rh_align->{cigar_string};
            #print "Cigar string is $cigar_string!\n";

            #my $first_op = 1;
            while ($cigar_string) {
                my ($bases, $op);
                if ($cigar_string =~ s/^(\d+)([MDINSHP])//) {
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
                            print STDERR "Setting ref_match within insertion to $current_ref\n" if ($self->{verbose});
                            last;
                        }
                    }
                    #elsif ($op eq 'S') {
                        #if ($first_op) {
                            #$current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                        #}
                    #}
                    #elsif ($op eq 'H') {
                        #if ($first_op) {
                            #$current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                        #}
#
                    #}
                    #$first_op = 0;
                }
                else {
                    die "Cigar string $cigar_string is of the wrong form!\n";
                }
            }
            if (!$align_match_found) { # check to be sure we got to the end
                print STDERR "Reached query $current_query (end is $query_end), ref $current_ref (end is $ref_end) without finding match\n" if ($self->{verbose});
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

    foreach my $rh_pairentry (@{$self->{entry_pairs}}) {
        next if ($ref && $rh_pairentry->{ref_entry} ne $ref);

        foreach my $rh_align (@{$rh_pairentry->{aligns}}) {
            my $query_start = $rh_align->{query_start};
            my $ref_start = $rh_align->{ref_start};
            my $query_end = $rh_align->{query_end};
            my $ref_end = $rh_align->{ref_end};
            my $comp = ($query_start < $query_end) ? 0 : 1;
            if (($ref_pos < $ref_start) || ($ref_pos > $ref_end)) {
                next;
            }

            my $current_ref = $ref_start - 1;
            my $current_query = ($comp) ? $query_start + 1 : $query_start - 1;

            my $cigar_string = $rh_align->{cigar_string};

            #my $first_op = 1;
            while ($cigar_string) {
                my ($bases, $op);
                if ($cigar_string =~ s/^(\d+)([MDINSHP])//) {
                    ($bases, $op) = ($1, $2);
                    if ($op eq 'M') {
                        if ($current_ref + $bases >= $ref_pos) {
                            $rh_align->{query_matches}->{$ref_pos} = ($comp) ? $current_query - ($ref_pos - $current_ref) : $current_query + ($ref_pos - $current_ref);
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
                    #elsif ($op eq 'S') {
                    #    if ($first_op) {
                    #        $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                    #    }
                    #}
                    #elsif ($op eq 'H') {
                    #    if ($first_op) {
                    #        $current_query = ($comp) ? $current_query - $bases : $current_query + $bases;
                    #    }
#
#                    }
                    #$first_op = 0;
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
  by the MiniMap2 show-diff program.

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
        my $query1 = $found_query_matches{$varcontig}->{'query1'};
        my $query2 = $found_query_matches{$varcontig}->{'query2'};
        my $comp1 = $found_query_matches{$varcontig}->{'comp1'};
        my $comp2 = $found_query_matches{$varcontig}->{'comp2'};
        if ($query1 && $query2 && $comp1 && $comp2 && ($comp1 == $comp2)) { # everything's here for this contig
            return ($varcontig, $query1, $query2, $comp1);
        }
    }

}

###########################################################

1;
__END__

=back
