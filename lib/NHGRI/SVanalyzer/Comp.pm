package NHGRI::SVanalyzer::Comp;
############################################################

=head1 NAME

Comp.pm - A Perl module to compare structural variant calls
          and report distance metrics between them.

=head1 DESCRIPTION

  Calls are passed to the module as hash references, along
  with a set of reference sequences. Methods are provided to
  align predicted alternate haplotypes to each other, and to
  report different distance metrics based on the resulting
  alignments and size differences between the alleles.

=head1 DATE

 July, 2017

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
use File::Temp qw/ tempdir /;

our $VERSION  = '0.01';

###########################################################

=over 4

=item B<new()>

  This method creates a Comp object from a set of reference
  sequences, passed as a valid reference fasta file path.

  Input:  -ref_fasta - path to a valid FASTA-formatted file
          containing the reference sequences that were used
          in calling structural variants
          -ref_db - path to a GTB::FASTA object (will be
          created from -ref_fasta file if not specified)
          -sv1_info - reference to a hash containing the
          first SV call (required)
          -sv2_info - reference to a hash containing the 
          second SV call (required)

  Output: New Comp object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $ref_fasta = $params{-ref_fasta};
    my $fasta_db = $params{-ref_db};
    my $rh_sv1 = $params{-sv1_info};
    my $rh_sv2 = $params{-sv2_info};

    my $self = { ref_fasta => $ref_fasta,
                 fasta_db => $fasta_db,
                 sv1_info => $rh_sv1,
                 sv2_info => $rh_sv2,
                  };

    if ((!$rh_sv1) || (!$rh_sv2)) {
        die "Two references to hashes of SV info must be passed as -sv1_info and -sv2_info to the NHGRI::SVanalyzer::Comp constructor!\n";
    }

    bless $self, $class;

    if (($ref_fasta) && (!$fasta_db)) {
        $self->{fasta_db} = GTB::FASTA->new($ref_fasta);
    }

    $self->_calculate_dependent_variables();

    return $self;
}

###########################################################

=item B<potential_match()>

  This method 
  
  Input:  -rel_size_diff - maximum allowable relative size difference
          -rel_shift - maximum allowable relative shift
          
  Output: 1 if there is a potential match between the 
          two SV's, 0 otherwise.

=cut

###########################################################
sub potential_match {
    my $self  = shift;
    my %params = @_;

    my $relshift = (defined($params{'-relshift'})) ? $params{'-relshift'} : 1.0;
    my $relsizediff = (defined($params{'-relsizediff'})) ? $params{'-relsizediff'} : 1.0;

    my $rh_sv1 = $self->{sv1_info};
    my $rh_sv2 = $self->{sv2_info};

    my $reflength1 = $rh_sv1->{reflength};
    my $reflength2 = $rh_sv2->{reflength};
    my $altlength1 = $rh_sv1->{altlength};
    my $altlength2 = $rh_sv2->{altlength};

    my $size1 = $altlength1 - $reflength1;
    my $size2 = $altlength2 - $reflength2;

    my $minsvsize = (abs($size1) < abs($size2)) ? abs($size1) : abs($size2);

    my $larger_allele1 = ($reflength1 > $altlength1) ? $reflength1 : $altlength1;
    my $larger_allele2 = ($reflength2 > $altlength2) ? $reflength2 : $altlength2;
    my $shared_denominator = ($larger_allele1 + $larger_allele2) / 2.0;

    if ((abs($size1 - $size2) > $shared_denominator * $relsizediff) ||
        (abs($size1 - $size2) > $shared_denominator * $relshift)) {
        return 0;
    }
    else {
        return 1;
    }

} ## end potential_match

###########################################################

=item B<calc_distance()>

  This method takes two hashes of structural variant info 
  and a Bio::DB::Sam::Fai object for the reference sequences,
  and returns a hash of different distance metrics 
  measuring how different the two variants are.

  Input:  If the second and third arguments are the key
          '-nocleanup' and a true value, will leave comparison
          fasta files and alignment output in place.
  Output: Reference to a hash containing distance metrics.

=cut

###########################################################
sub calc_distance {
    my $self  = shift;
    my %params = @_;

    my $rh_sv1 = $self->{sv1_info};
    my $rh_sv2 = $self->{sv2_info};

    my $ref1 = $rh_sv1->{ref};
    my $ref2 = $rh_sv2->{ref};
    my $alt1 = $rh_sv1->{alt};
    my $alt2 = $rh_sv2->{alt};
    my $chrom = $rh_sv1->{chrom};
    my $pos1 = $rh_sv1->{pos};
    my $pos2 = $rh_sv2->{pos};
    my $end1 = $rh_sv1->{end};
    my $end2 = $rh_sv2->{end};
    my $id1 = $rh_sv1->{id};
    my $id2 = $rh_sv2->{id};
    my $reflength1 = $rh_sv1->{reflength};
    my $reflength2 = $rh_sv2->{reflength};
    my $altlength1 = $rh_sv1->{altlength};
    my $altlength2 = $rh_sv2->{altlength};

    my $size1 = $altlength1 - $reflength1;
    my $size2 = $altlength2 - $reflength2;

    my $minsvsize = (abs($size1) < abs($size2)) ? abs($size1) : abs($size2);

    my $larger_allele1 = ($reflength1 > $altlength1) ? $reflength1 : $altlength1;
    my $larger_allele2 = ($reflength2 > $altlength2) ? $reflength2 : $altlength2;
    my $shared_denominator = ($larger_allele1 + $larger_allele2) / 2.0;

    # boundaries of constructed haplotype:

    my $left_bound = ($pos1 < $pos2) ? $pos1 : $pos2;
    my $right_bound = ($end1 < $end2) ? $end2 : $end1;
    my $rs_alt_hap1 = $self->construct_alt_hap(-left_bound => $left_bound, -right_bound => $right_bound, -svhashref => $rh_sv1);
    my $rs_alt_hap2 = $self->construct_alt_hap(-left_bound => $left_bound, -right_bound => $right_bound, -svhashref => $rh_sv2);
    my $minhaplength = (length(${$rs_alt_hap1}) < length(${$rs_alt_hap2})) ? length(${$rs_alt_hap1}) : length(${$rs_alt_hap2});
    my $althaplength_avg = (length(${$rs_alt_hap1}) + length(${$rs_alt_hap2}))/2.0;
    my $althaplength_diff = length(${$rs_alt_hap1}) - length(${$rs_alt_hap2});

    if (${$rs_alt_hap1} eq ${$rs_alt_hap2}) { # identical variants
        my $matchtype = 'EXACTMATCH';
        return {'edit_distance' => 0,
                'max_shift' => 0,
                'match_type' => $matchtype,
                'altlength_diff' => $althaplength_diff,
                'altlength_avg' => $althaplength_avg,
                'minhaplength' => $minhaplength,
                'size_diff' => $size2 - $size1,
                'size_avg' => ( $size1 + $size2 )/2.0,
                'shared_denominator' => $shared_denominator};
    }
    else {
        # align the alternative haplotypes to each other and evaluate
        my ($maxshift, $editdistance) = $self->compare_alt_haplotypes($rs_alt_hap1, $rs_alt_hap2, $id1, $id2, '-nocleanup',  $params{'-nocleanup'});
        #if (($editdistance/$minhaplength < 0.05) && (abs($maxshift) < $minsvsize)) {
        if (abs($maxshift) < $minsvsize) {
            my $matchtype = 'NWMATCH';
            return {'edit_distance' => $editdistance,
                    'max_shift' => $maxshift,
                    'match_type' => $matchtype,
                    'altlength_diff' => $althaplength_diff,
                    'altlength_avg' => $althaplength_avg,
                    'minhaplength' => $minhaplength,
                    'size_diff' => $size2 - $size1,
                    'size_avg' => ( $size1 + $size2 )/2.0,
                    'shared_denominator' => $shared_denominator};
        }
        else {
            my $matchtype = 'NWFAIL';
            return {'edit_distance' => $editdistance,
                    'max_shift' => $maxshift,
                    'match_type' => $matchtype,
                    'altlength_diff' => $althaplength_diff,
                    'altlength_avg' => $althaplength_avg,
                    'minhaplength' => $minhaplength,
                    'size_diff' => $size2 - $size1,
                    'size_avg' => ( $size1 + $size2 )/2.0,
                    'shared_denominator' => $shared_denominator};
        }
    }

} ## end calc_distance

###########################################################

=item B<construct_alt_hap()>

  This method constructs the alternate haplotype of a 
  structural variant between specified start and endpoints
  in the reference genome for the object.

  The following parameters should be passed to the method:
  
  Input:  -left_bound - coordinate of 5'-most reference base
          to be included in the haplotype sequence
          -right_bound - coordinate of the 3'-most reference
          base to be included in the haplotype sequence
          -svhashref - reference to a hash containing the SV
          call

  Output: reference to string containing the alternate haplotype.

=cut

###########################################################
sub construct_alt_hap {
    my $self  = shift;
    my %params = @_;

    my $rh_sv = $params{-svhashref};

    my $left_bound = $params{-left_bound};
    my $right_bound = $params{-right_bound};

    my $chrom = $rh_sv->{chrom};
    my $alt_allele = uc($self->{fasta_db}->seq("$chrom:$left_bound-$right_bound")); # widened ref
    my $svstart = $rh_sv->{'pos'};
    my $offset = $svstart - $left_bound; # number of positions to move from first included base to ref start
    my $svend = $rh_sv->{'end'};
    my $length = $svend - $svstart + 1;

    substr($alt_allele, $offset, $length) = $rh_sv->{'alt'};

    return \$alt_allele;

} ## end construct_alt_hap

###########################################################

=item B<compare_alt_haplotypes()>

  This method aligns alternate haplotype sequences from two 
  different variants, and calculates distance metrics between
  them.

  Input: Two references to alternate haplotype strings, 
         followed by two optional id strings (to be used in
         temporary file names). Finally, a hash of parameters
         can be passed with a "-nocleanup" option.
  Output: A hash of distance metrics

=cut

###########################################################
sub compare_alt_haplotypes {
    my $self  = shift;
    my $rs_alt1 = shift;
    my $rs_alt2 = shift;
    my $id1 = shift || 'seq1';
    my $id2 = shift || 'seq2';
    my %params = @_;

    if (!$params{'-nocleanup'}) {
        $params{'-nocleanup'} = 0;
    }

    my $workingdir = tempdir( CLEANUP => !($params{'-nocleanup'}), DIR => '.');
    my $tmpfasta1 = "$workingdir/$id1.$id2.fa";
    $self->write_fasta_file($tmpfasta1, "$id1\_alt", $rs_alt1);
    my $tmpfasta2 = "$workingdir/$id2.$id1.fa";
    $self->write_fasta_file($tmpfasta2, "$id2\_alt", $rs_alt2);

    my $edlib_aligner = 'edlib-aligner';
    my $nw_output = `$edlib_aligner -p -f CIG_STD $tmpfasta1 $tmpfasta2`;   
    unlink $tmpfasta1 unless ($params{'-nocleanup'});
    unlink $tmpfasta2 unless ($params{'-nocleanup'});
    if ($nw_output =~ /Cigar:\n(.*)/m) {
        my $cigar_string = $1;
        my $score = ($nw_output =~ /score = (\d+)/)  ? $1 : 'NA';
        my $maxshift = calc_max_shift($cigar_string);
        return ($maxshift, $score);
    }
    else {
        return ('NA', 'NA');
    }

} ## end compare_alt_haplotypes

###########################################################

=item B<_calculate_dependent_variables()>

  This method calculates the values of reflength, altlength,
  and size for each SV of the Comp object.

  Input: Comp object
  Output: 1 if successful

=cut

###########################################################
sub _calculate_dependent_variables {
    my $self  = shift;
   
    $self->{sv1_info}->{reflength} = length($self->{sv1_info}->{ref}); 
    $self->{sv2_info}->{reflength} = length($self->{sv2_info}->{ref}); 
    $self->{sv1_info}->{altlength} = length($self->{sv1_info}->{alt}); 
    $self->{sv2_info}->{altlength} = length($self->{sv2_info}->{alt}); 
    $self->{sv1_info}->{size} = $self->{sv1_info}->{reflength} - $self->{sv1_info}->{altlength};
    $self->{sv2_info}->{size} = $self->{sv2_info}->{reflength} - $self->{sv2_info}->{altlength};

    return 1;

} ## end _calculate_dependent_variables

###########################################################

=item B<write_fasta_file()>

  This method writes the specified FASTA-formatted sequence 
  to the specified file with the specified sequence id
  (yes, I know I should have called someone else's code to 
  do this).

  Input: Comp object, then FASTA file path (must be writeable),
         sequence id (string), and reference to sequence.
  Output: 1 if successful, dies otherwise.

=cut

###########################################################
sub write_fasta_file {
    my $self  = shift;
    my $fasta_file = shift;
    my $seq_id = shift;
    my $r_sequence = shift;
    
    my $seq_fh = Open("$fasta_file", "w");
    print $seq_fh ">$seq_id\n";

    my $r_alt_50 = format_50($r_sequence);
    print $seq_fh ${$r_alt_50};

    if (${$r_alt_50} !~ /\n$/) {
        print $seq_fh "\n";
    }

    close $seq_fh;

    return 1;

} ## end write_fasta_file

###########################################################

=item B<format_50()>

  This method takes a sequence string and splits it into
  lines 50 bases long.

  Input: Comp object, then reference to sequence string.
  Output: reference to formatted sequence string

=cut

###########################################################
sub format_50 {
    my $r_seq = shift;
    
    my $bases = 0;
    my $revseq = reverse(${$r_seq});
    my $formattedseq = '';
    while (my $nextbase = chop $revseq) {
        $formattedseq .= $nextbase;
        $bases++;
        if ($bases == 50) {
            $formattedseq .= "\n";
            $bases = 0;
        }
    }
    if ($bases) { # need an extra enter
        $formattedseq .= "\n";
    }

    return \$formattedseq;

} ## end format_50

###########################################################

=item B<calc_max_shift()>

  This method takes a cigar string as an argument and
  returns the maximum number of bases the alignment shifts
  from the diagonal.

  Input: Cigar string
  Output: Integer

=cut

###########################################################
sub calc_max_shift {
    my $cigar = shift;

    my ($max_shift, $current_shift) = (0, 0);
    #print "Cigar: $cigar\n";
    while ($cigar) {
        my ($bases, $op);
        if ($cigar =~ s/^(\d+)([MDI])//) {
            ($bases, $op) = ($1, $2);
            if ($op eq 'D') { # deleted from reference, decrease shift
                $current_shift -= $bases;
            }
            elsif ($op eq 'I') { # deleted from reference--advance ref coord
                $current_shift += $bases;
            }
            if (abs($current_shift) > abs($max_shift)) {
                $max_shift = $current_shift;
            }
        }
        else {
            die "Cigar string $cigar is of the wrong form!\n";
        }
    }
    #print "Max shift: $max_shift\n";

    return $max_shift;

} ## end calc_max_shift

###########################################################

1;

__END__

=back
