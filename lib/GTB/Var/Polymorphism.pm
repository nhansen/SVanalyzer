package GTB::Var::Polymorphism;
############################################################

=head1 NAME

Polymorphism.pm - A Perl module to contain information
          about genomic polymorphisms (snps and dips).

=head1 DESCRIPTION

  A Perl module for polymorphism characterization.

=head1 DATE

 April 11, 2007

=head1 AUTHORS

 Nancy Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=cut

############################################################
use strict;
use warnings;

use Carp;

our $VERSION  = '0.01';

use vars qw( );

###########################################################

=over 4

=item new()

  This method creates a Polymorphism object.

  Input:  -left_flank_end - coordinate of last reference
             base in the left flanking sequence.
          -right_flank_start - coordinate of the first
             reference base in the right flanking seq.
          -left_flank_seq - sequence flanking on the left
             side of the polymorphic region
          -right_flank_seq - sequence flanking on the right
             side of the polymorphic region
          -allele_seqs - reference to an array of sequences
             one for each allele, which represent the
             sequence in that allele between the left and
             right flanking sequences.  The first entry of
             the array is the reference sequence.
          -type - the type of polymorphism
          -score - a score predicting how likely this 
             polymorphism is to be real.
          -reported_position - position of this polymorphism
             as reported by calling program. 
          -reported_size - size of this polymorphism
             as reported by calling program (1 for SNPs)
          -amplimer_id - the id of the amplimer in which this
             polymorphism was called (if there is one).
          -assembly_id - the id value for the assembly this
             polymorphism was predicted in.
          -analysis_id - the id value for the analysis this
             polymorphism was predicted in.
          -genotypes - reference to a list of genotypes for
             this polymorphism (will be set to [] if not
             provided.)
          -strict_boundaries - if true, the object will
             maintain strict boundaries for a polymorphism
             even if the repetitive nature of the variation
             allows generalization (default false).

  Output: New Polymorphism object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $polymorphism_id = $params{-polymorphism_id};
    my $left_flank_end = $params{-left_flank_end};
    my $right_flank_start = $params{-right_flank_start};
    my $right_flank_seq = $params{-right_flank_seq};
    my $left_flank_seq = $params{-left_flank_seq};
    my $ra_allele_seqs = $params{-allele_seqs} || [];
    my $reported_position = $params{-reported_position};
    my $reported_size = $params{-reported_size};
    my $score = $params{-score};
    my $type = $params{-type};
    my $amplimer_id = $params{-amplimer_id};
    my $assembly_id = $params{-assembly_id};
    my $analysis_id = $params{-analysis_id};
    my $reference = $params{-reference};
    my $ra_genotypes = $params{-genotypes} || [];

    my $dummy_var = map {tr/a-z/A-Z/} @{$ra_allele_seqs};

    foreach my $seq (@{$ra_allele_seqs})
    {
        if ($seq =~ /[a-z]/) # something went wrong!
        {
             croak "Your funky translation didn\'t work!\n";
        }
    }

    my $rh_genotypes = {};
    foreach my $genotype (@{$ra_genotypes})
    {
        my $read_name = $genotype->read_name();
        $rh_genotypes->{$read_name} = $genotype;
    }
    
    my $strict_boundaries = $params{-strict_boundaries} || 0;
 
    my $self = {polymorphism_id => $polymorphism_id,
                left_flank_end => $left_flank_end,
                right_flank_start => $right_flank_start,
                left_flank_seq => $left_flank_seq,
                right_flank_seq => $right_flank_seq,
                allele_seqs => $ra_allele_seqs,
                score => $score,
                type => $type,
                reported_position => $reported_position,
                reported_size => $reported_size,
                amplimer_id => $amplimer_id,
                assembly_id => $assembly_id,
                analysis_id => $analysis_id,
                reference => $reference,
                genotype_hash => $rh_genotypes};

    bless $self, $class;

    unless ($strict_boundaries)
    {
        $self->_adjust_boundaries();
    }

    return $self;
}

###########################################################

=item polymorphism_id()

  This method gets or sets the value of polymorphism_id, which
  is the id value for this Polymorphism in the tgs_polymorphism
  table.

  Input: Optional argument sets value.
  Output: Value of "polymorphism_id" (scalar number).

=cut

###########################################################
sub polymorphism_id {
    my $self  = shift;

    if (defined (my $new_polymorphism_id = shift))
    {
        $self->{polymorphism_id} = $new_polymorphism_id;
    }

    return $self->{polymorphism_id};

} ## end polymorphism_id

###########################################################

=item reference()

  This method gets or sets the value of reference, which
  is a string defining the reference sequence.

  Input: Optional argument sets value.
  Output: Value of "reference" (scalar string).

=cut

###########################################################
sub reference {
    my $self  = shift;

    if (defined (my $new_reference = shift))
    {
        $self->{reference} = $new_reference;
    }

    return $self->{reference};

} ## end reference

###########################################################

=item amplimer_id()

  This method gets or sets the value of amplimer_id, which
  is a string defining the amplimer_id sequence.

  Input: Optional argument sets value.
  Output: Value of "amplimer_id" (scalar string).

=cut

###########################################################
sub amplimer_id {
    my $self  = shift;

    if (defined (my $new_amplimer_id = shift))
    {
        $self->{amplimer_id} = $new_amplimer_id;
    }

    return $self->{amplimer_id};

} ## end amplimer_id

###########################################################

=item assembly_id()

  This method gets or sets the value of assembly_id, which
  is the primary key value of the assembly record in tgs_assembly.

  Input: Optional argument sets value.
  Output: Value of "assembly_id" (scalar number).

=cut

###########################################################
sub assembly_id {
    my $self  = shift;

    if (defined (my $new_assembly_id = shift))
    {
        $self->{assembly_id} = $new_assembly_id;
    }

    return $self->{assembly_id};

} ## end assembly_id

###########################################################

=item analysis_id()

  This method gets or sets the value of analysis_id, which
  is the primary key value of the analysis record in tgs_analysis.

  Input: Optional argument sets value.
  Output: Value of "analysis_id" (scalar number).

=cut

###########################################################
sub analysis_id {
    my $self  = shift;

    if (defined (my $new_analysis_id = shift))
    {
        $self->{analysis_id} = $new_analysis_id;
    }

    return $self->{analysis_id};

} ## end analysis_id

###########################################################

=item type()

  This method gets or sets the value of type, which
  is the type of polymorphism (e.g., substitution, 
  insertion, or deletion).

  Input: Optional argument sets value.
  Output: Value of "type" (scalar string).

=cut

###########################################################
sub type {
    my $self  = shift;

    if (defined (my $new_type = shift))
    {
        $self->{type} = $new_type;
    }

    return $self->{type};

} ## end type

###########################################################

=item score()

  This method gets or sets the value of score, which
  is a number indicating how likely this polymorphism
  is to be real.

  Input: Optional argument sets value.
  Output: Value of "score", rounded to nearest integer 
        (scalar number).

=cut

###########################################################
sub score {
    my $self  = shift;

    if (defined (my $new_score = shift))
    {
        $self->{score} = $new_score;
    }

    return int($self->{score} + 0.5);

} ## end score

###########################################################

=item recalculate_score()

  This method recalculates and resets the polymorphism
  score based on the genotype scores.  Presently, it 
  assigns the polymorphism the highest genotype score for
  a non-reference genotype. 

  Input: None.
  Output: New value of "score", rounded to nearest integer 
        (scalar number).

=cut

###########################################################
sub recalculate_score {
    my $self  = shift;
    my %params = @_;
    my $verbose = $params{-verbose};

    my $ref_allele = $self->allele_seqs()->[0];
    print "Ref allele is $ref_allele\n" if ($verbose);
    my $high_score;
    foreach my $genotype (@{$self->genotypes()})
    {
        next if (($genotype->allele1() eq $ref_allele) &&
                 ($genotype->allele2() eq $ref_allele));

        my $gen_score = $genotype->score();
        my $read_name = $genotype->read_name();
        print "Non ref $read_name has score $gen_score\n" if ($verbose);
        if (!$high_score || $gen_score > $high_score)
        {
            $high_score = $gen_score;
        }
    }

    $self->score($high_score);

    return $high_score;

} ## end recalculate_score

###########################################################

=item reported_position()

  This method gets or sets the value of reported_position, which
  is the position of this polymorphism as reported by the 
  poly-calling program.

  Input: Optional argument sets value.
  Output: Value of "reported_position" (scalar string).

=cut

###########################################################
sub reported_position {
    my $self  = shift;

    if (defined (my $new_reported_position = shift))
    {
        $self->{reported_position} = $new_reported_position;
    }

    return $self->{reported_position};

} ## end reported_position

###########################################################

=item reported_size()

  This method gets or sets the value of reported_size, which
  is the size of this polymorphism as reported by the 
  poly-calling program.

  Input: Optional argument sets value.
  Output: Value of "reported_size" (scalar string).

=cut

###########################################################
sub reported_size {
    my $self  = shift;

    if (defined (my $new_reported_size = shift))
    {
        $self->{reported_size} = $new_reported_size;
    }

    return $self->{reported_size};

} ## end reported_size

###########################################################

=item left_flank_end()

  This method gets or sets the value of left_flank_end, which
  is the location of the last reference base in the left
  flanking sequence.

  Input: Optional argument sets value.
  Output: Value of "left_flank_end" (scalar number).

=cut

###########################################################
sub left_flank_end {
    my $self  = shift;

    if (defined (my $new_left_flank_end = shift))
    {
        $self->{left_flank_end} = $new_left_flank_end;
    }

    return $self->{left_flank_end};

} ## end left_flank_end

###########################################################

=item right_flank_start()

  This method gets or sets the value of right_flank_start, which
  is the location of the last reference base in the left
  flanking sequence.

  Input: Optional argument sets value.
  Output: Value of "right_flank_start" (scalar number).

=cut

###########################################################
sub right_flank_start {
    my $self  = shift;

    if (defined (my $new_right_flank_start = shift))
    {
        $self->{right_flank_start} = $new_right_flank_start;
    }

    return $self->{right_flank_start};

} ## end right_flank_start

###########################################################

=item left_flank_seq()

  This method gets or sets the value of left_flank_seq, which
  is the sequence to the left of the polymorphic region.

  Input: Optional argument sets value.
  Output: Value of "left_flank_seq" (scalar string).

=cut

###########################################################
sub left_flank_seq {
    my $self  = shift;

    if (defined (my $new_left_flank_seq = shift))
    {
        $self->{left_flank_seq} = $new_left_flank_seq;
    }

    return $self->{left_flank_seq};

} ## seq left_flank_seq

###########################################################

=item right_flank_seq()

  This method gets or sets the value of right_flank_seq, which
  is the sequence to the right of the polymorphic region.

  Input: Optional argument sets value.
  Output: Value of "right_flank_seq" (scalar string).

=cut

###########################################################
sub right_flank_seq {
    my $self  = shift;

    if (defined (my $new_right_flank_seq = shift))
    {
        $self->{right_flank_seq} = $new_right_flank_seq;
    }

    return $self->{right_flank_seq};

} ## seq right_flank_seq

###########################################################

=item allele_seqs()

  This method gets or sets the value of allele_seqs, which
  is a reference to an array containing the various 
  sequences that could potentially be between the reference
  bases left_flank_end and right_flank_start.

  Input: Optional argument sets value.
  Output: Value of "allele_seqs" (reference to an array of
          scalar strings).

=cut

###########################################################
sub allele_seqs {
    my $self  = shift;

    if (defined (my $new_allele_seqs = shift))
    {
        $self->{allele_seqs} = $new_allele_seqs;
    }

    return $self->{allele_seqs};

} ## end allele_seqs

###########################################################

=item allele_aliases()

  This method returns a reference to a hash containing
  "aliases" for each of the alleles for this polymorphism.
  For DIP polymorphisms, these aliases are currently a
  representation of each allele in the form n(X)Y, where 
  n is a number of repetitions, X is a repeated sequence,
  and Y is a sequence that is a subsequence of X. 

  Input: None.
  Output: Value of "allele_aliases" (reference to a hash of
          scalar strings).

=cut

###########################################################
sub allele_aliases {
    my $self  = shift;

    if (!defined($self->{allele_aliases}))
    {
        my $rh_allele_seqs = {};
        for (my $i=0; $i<=$#{$self->{allele_seqs}}; $i++)
        {
            my $allele = $self->{allele_seqs}->[$i];
            # default is to keep allele the same:
            $rh_allele_seqs->{$allele} = $allele;
            next if ($allele eq '*');

            for (my $i=1; $i<=length $allele; $i++)
            {
                my $tmp_allele = $allele;
                my $test_string = substr($tmp_allele, 0, $i);

                my $rep = 0;
                while ($tmp_allele =~ s/^$test_string//)
                {
                    $rep++;
                }
                if ($test_string =~ /^$tmp_allele/)
                {
                    #print "$allele is $rep reps of $test_string + $tmp_allele\n";
                    $rh_allele_seqs->{$allele} = ($rep==1) ? "$test_string$tmp_allele" : "$rep\($test_string\)$tmp_allele";
                    last;
                }
            }
        }

        $self->{allele_aliases} = $rh_allele_seqs;
    }

    return $self->{allele_aliases};

} ## allele_aliases

###########################################################

=item jamies_allele_aliases()

  This method returns a reference to a hash containing
  "aliases" for each of the alleles for this polymorphism,
  according to the code in Jamie's find_smallest_repeat.pl

  Input: None.
  Output: Value of "allele_aliases" (reference to a hash of
          scalar strings).

=cut

###########################################################
sub jamies_allele_aliases {
    my $self  = shift;

    my $rh_allele_seqs = {};
    for (my $i=0; $i<=$#{$self->{allele_seqs}}; $i++)
    {
        my $allele = $self->{allele_seqs}->[$i];
        next if ($rh_allele_seqs->{$allele});

        my $allele_len = length $allele;
        my @uniq_units;
        my %units_like;
        my $min_unit_size = $allele_len;
        
        for my $i (1..$allele_len) {
        
            # At a given chunk size, find all units - load to hash
            for (my $j = 1; ($j - ($allele_len % $i)) <= $allele_len; $j+=$i) {
            
                my $chunk;
                {
                    no warnings;
                    $chunk = substr $allele, $j-1, $i;
                }
                
                if ($chunk) {
                    $units_like{$i}{$chunk}++;
                }
            }
            
        
            # Look at hash of chunks, figure out which has the smallest number of chunks (not counting "spill over" / modulo)
            my $total_chunk = 0;
            my $max_chunk = 0;
            my $chunk_num = 0;
        
            for my $chunk (keys %{$units_like{$i}}) {
        
                # ignore modulo "spill-over"
                next if ( length($chunk) == ($allele_len % $i) );
        
                $chunk_num++;
                $total_chunk += $units_like{$i}{$chunk};
                if ($units_like{$i}{$chunk} > $max_chunk) {
                    $max_chunk = $units_like{$i}{$chunk};
                }
            }
        
            if ($chunk_num == 1 && $i < $min_unit_size && $max_chunk > 1) {
                $min_unit_size = $i;
            }
        }            
        for my $chnk (reverse sort {length $a <=> length $b} keys %{$units_like{$min_unit_size}}) {

            my $string = ($units_like{$min_unit_size}{$chnk} == 1) ? 
                           $chnk : "$units_like{$min_unit_size}{$chnk}($chnk)";
            $rh_allele_seqs->{$allele} .= $string;
        }
    }

    return $rh_allele_seqs;

} ## end jamies_allele_aliases

###########################################################

=item reference_allele()

  This method returns the first element of the "allele_seqs"
  array for this object.

  Input: None.
  Output: Reference allele (scalar string).

=cut

###########################################################
sub reference_allele {
    my $self  = shift;

    my $ra_allele_seqs = $self->allele_seqs();

    return $ra_allele_seqs->[0] || 'None';

} ## end reference_allele

###########################################################

=item non_ref_genotypes()

  This method returns a reference to a list of the 
  non-reference genotypes for this polymorphism object.

  Input: None.
  Output: Reference to an array of GTB::Var::Genotype 
       objects.

=cut

###########################################################
sub non_ref_genotypes {
    my $self  = shift;

    my $reference_allele = $self->reference_allele();
    my $ra_genotypes = $self->genotypes();

    my @non_ref_genotypes = ();

    foreach my $genotype (@{$ra_genotypes})
    {
        if (($genotype->allele1() !~ /^$reference_allele$/i) || 
            ($genotype->allele2() !~ /^$reference_allele$/i))
        {
            push @non_ref_genotypes, $genotype;
        }
    }

    return \@non_ref_genotypes;

} ## end non_ref_genotypes

###########################################################

=item genotypes()

  This method gets or sets the value of genotypes, which
  is a reference to an array containing Genotype objects
  for any known genotypes for this polymorphism.

  Input: None.
  Output: Value of "genotypes" (reference to an array of
          GTB::Var::Genotype objects).

=cut

###########################################################
sub genotypes {
    my $self  = shift;

    my @genotypes = values %{$self->genotype_hash()};
    return \@genotypes;

} ## end genotypes

###########################################################

=item genotype_hash()

  This method gets or sets the value of genotype_hash, which
  is a reference to a hash containing read_name, Genotype
  object pairs.

  Input: Optional argument sets value.
  Output: Value of "genotype_hash" (reference to a hash
          whose keys are read names and values are 
          GTB::Var::Genotype objects).

=cut

###########################################################
sub genotype_hash {
    my $self  = shift;

    if (defined (my $new_genotype_hash = shift))
    {
        $self->{genotype_hash} = $new_genotype_hash;
    }

    return $self->{genotype_hash} || {};

} ## end genotype_hash

###########################################################

=item add_genotype()

  This method adds a genotype to the genotype hash.

  Input: GTB::Var::Genotype object.
  Output: 1 if successful.

=cut

###########################################################
sub add_genotype {
    my $self  = shift;
    my $genotype = shift;

    my $readname = $genotype->read_name();
    $self->{genotype_hash}->{$readname} = $genotype;

    return 1;
}

###########################################################

=item poly_file_strings()

  This method returns a reference to a hash with two
  strings, a "polyfile_string", which is the string for
  this polymorphism that should be written to Pedro\'s
  "poly_info" file, and a "polycompfile_string", which 
  is the string to be written to the file with the
  polymorphism components for this polymorphism.

  Input: -project_id passes the project id to be written
         to the strings
         -analysis_id passes the analysis id
         -trace_ids passes a reference to a hash of 
         trace_id values, where the keys are the trace
         names from the assembly.
  Output: Array containing the poly file string and the
         poly comp file string.

=cut

###########################################################
sub poly_file_strings {
    my $self  = shift;
    my %params = @_;

    my $amplimer_id = $params{'-amplimer_id'} || $self->amplimer_id() || 'NA';
    my $project_id = $params{'-project_id'} || 'NA';
    my $analysis_id = $params{'-analysis_id'} || 'NA';
    my $rh_trace_ids = $params{'-trace_ids'};

    my $reported_position = $self->reported_position;
    my $left_flank_end = $self->left_flank_end;
    my $right_flank_start = $self->right_flank_start;

    my $allele1 = $self->allele_seqs()->[0] || '-';
    my $allele2 = $self->allele_seqs()->[1] || '-';

    my $type = $self->type();
    my $score = $self->score();

    my $type_char = ($type eq 'SNP') ? '' :
                    ($type eq 'Deletion') ? 'D' :
                    ($type eq 'Insertion') ? 'I' : 'ID'; 

    my $poly_string = "$amplimer_id:$reported_position$type_char\t$project_id\t$analysis_id\tAMPLIMER\t$amplimer_id\t$left_flank_end\t$right_flank_start\t$type\t$allele1\t$allele2\t$score\n";

    my $polycomp_string = '';
    foreach my $genotype (@{$self->genotypes()})
    {
        my $readname = $genotype->read_name();
        next if ($readname =~ /\.REF\./);
        my $trace_id = $rh_trace_ids ? $rh_trace_ids->{$readname} 
                                     : $readname; # will return undef if undefined in hash, but readname if no hash
        if (!defined ($trace_id))
        {
            croak "No trace id defined for $readname in hash passed as -trace_ids to poly_file_strings!\n";
        }
        my $sample_allele1 = $genotype->allele1() || '-';
        my $sample_allele2 = $genotype->allele2() || '-';

        next if ($sample_allele1 eq 'W' || $sample_allele2 eq 'W'); # don't want to print NN genotypes

        my $gen_score = $genotype->score();
        my $trace_allele1_left_end = $genotype->trace_allele1_left_end() || 0;
        my $trace_allele1_right_start = $genotype->trace_allele1_right_start() || 0;
        my $trace_allele2_left_end = $genotype->trace_allele2_left_end() || 0;
        my $trace_allele2_right_start = $genotype->trace_allele2_right_start() || 0;
        $polycomp_string .= "$amplimer_id:$reported_position$type_char\t$trace_id\t$analysis_id\t$sample_allele1\t$sample_allele2\t$gen_score\t$trace_allele1_left_end\t$trace_allele1_right_start\t$trace_allele2_left_end\t$trace_allele2_right_start\n";
    }

    return ($poly_string, $polycomp_string);

} ## file_strings

###########################################################

=item genotype_alleles()

  This method returns a list of the alleles seen in
  genotypes for this object.

  Input: -include_reference => true value will include 
         the reference allele, whether or not there is 
         a genotype that contains it.
  Output: List of allele seqs (Array of scalar strings).

=cut

###########################################################
sub genotype_alleles {
    my $self  = shift;
    my %params = @_;
    my $include_ref = defined ($params{-include_reference}) ? 
                      $params{-include_reference} : 0;

    my %alleles_here = ();
    if ($include_ref)
    {
        $alleles_here{$self->allele_seqs()->[0]} = 1;
    }

    foreach my $genotype (@{$self->genotypes()})
    {
        $alleles_here{$genotype->allele1()} = 1;
        $alleles_here{$genotype->allele2()} = 1;
    }

    my @all_alleles = values %alleles_here;
    return @all_alleles;

} ## end genotype_alleles

###########################################################

=item _adjust_boundaries()

  This method looks to see if this variation is non-unique
  with respect to translation of alleles, and enlarges the
  boundaries to make a unique specification of the variation.
  If empty alleles are passed as '*' character, they will be
  handled correctly, and returned as '*' if not expanded.

  Input: None.
  Output: 1 if object is altered, 0 otherwise.

=cut

###########################################################
sub _adjust_boundaries {
    my $self  = shift;
    my $left_flank_seq = $self->left_flank_seq();
    my $right_flank_seq = $self->right_flank_seq();

    my $return_value = 0;
  
    my $type = $self->type();

    #print "Type $type\n";
    return $return_value if ($type !~ /^(insertion|deletion|indel)$/i); # only expanding indels
    #print "$left_flank_seq\n$right_flank_seq\n";
    return $return_value if (!$left_flank_seq && !$right_flank_seq); # can't expand without sequence

    # assuming only two alleles for now--easiest way!
    my ($first_allele, $second_allele) = @{$self->allele_seqs()};
    $first_allele =~ s/\*//;
    $second_allele =~ s/\*//;
    #print "Original alleles: $first_allele, $second_allele\n";

    # first flush left and attempt to expand right:
    my $right_flank_start = $self->right_flank_start();
    my $rev_right_flank = reverse $right_flank_seq;

    # store original first and second alleles in case we need to alter genotypes:
    my ($orig_first_allele, $orig_second_allele) = ($first_allele, $second_allele);
    my $okay = 1;
    while (($okay) && (my $next_base = chop $rev_right_flank))
    {
        my $first_free_base = _get_first_free_base($first_allele, $second_allele, 'left');
        if (!$first_free_base) # alleles are not shiftable
        {
            $okay = 0;
            next;
        }
        if ($next_base =~ /^$first_free_base$/i)
        {
            $return_value = 1;
            $first_allele .= $next_base;
            $second_allele .= $next_base;
            $right_flank_start += 1;
            $right_flank_seq =~ s:^.::;
        }
        else
        {
            $okay = 0;
        }
    }

    # now flush right and attempt to expand left:
    my $left_flank_end = $self->left_flank_end();
    my $disp_left_flank = $left_flank_seq;

    $okay = 1;
    while (($okay) && (my $next_base = chop $disp_left_flank))
    {
        my $first_free_base = _get_first_free_base($first_allele, $second_allele, 'right');
        if (!$first_free_base) # alleles are not shiftable
        {
            $okay = 0;
            next;
        }
        if ($next_base =~ /^$first_free_base$/i)
        {
            $return_value = 1;
            $first_allele = "$next_base$first_allele";
            $second_allele = "$next_base$second_allele";
            $left_flank_end -= 1;
            chop $left_flank_seq;
        }
        else
        {
            $okay = 0;
        }
    }

    #print "New alleles: $first_allele, $second_allele\n";
    if ($return_value) # need to alter the object
    {
        $self->left_flank_end($left_flank_end);
        $self->right_flank_start($right_flank_start);
        $self->left_flank_seq($left_flank_seq);
        $self->right_flank_seq($right_flank_seq);
        $self->allele_seqs([$first_allele, $second_allele]);
        my $new_reported_position = ($type =~ /ins/i) ? $left_flank_end : $left_flank_end + 1;
        $self->reported_position($new_reported_position);

        # alter genotypes:
        my @all_genotypes = @{$self->genotypes()};
        $self->genotype_hash({}); # empty out genotypes--will replace with altered ones.
        foreach my $genotype (@all_genotypes)
        {
            my $new_allele1 = ($genotype->allele1() eq $orig_first_allele) ? $first_allele :
                              ($genotype->allele1() eq $orig_second_allele) ? $second_allele : 
                              $genotype->allele1(); # leave unchanged if not an original allele
            my $new_allele2 = ($genotype->allele2() eq $orig_first_allele) ? $first_allele : 
                              ($genotype->allele2() eq $orig_second_allele) ? $second_allele : 
                              $genotype->allele2(); # leave unchanged if not an original allele
            $genotype->allele1($new_allele1);
            $genotype->allele2($new_allele2);
            $self->add_genotype($genotype);
        }
    }

    return $return_value;

} ## end _adjust_boundaries

###########################################################

=item _get_first_free_base()

  This method takes two allele strings as arguments, and
  returns the first unaligned base.

  Input: Two allele strings, then 'left' or 'right'
      to specify how to flush the strings.
  Output: First unaligned base (scalar character)

=cut

###########################################################
sub _get_first_free_base {
    my $first_allele = shift;
    my $second_allele = shift;
    my $side = shift;

    my $first_seq = ($side eq 'left') ? reverse $first_allele : $first_allele;
    my $second_seq = ($side eq 'left') ? reverse $second_allele : $second_allele;

    while (my $first_base = chop $first_seq)
    {
        my $second_base = chop $second_seq;

        return $first_base if (!$second_base);
        return $second_base if (!$first_base);

        return '' if ($first_base !~ /^$second_base$/i);
    }

    # if there is no sequence in first_seq, need to return second seq if it's there:
    my $return_base = chop $second_seq || '';
    return $return_base;

} ## end _get_first_free_base

###########################################################

=item reconcile_genotypes()

  This method examines each genotype, grouping them 
  according to their source, then alters genotypes as
  necessary so all genotypes from a single source agree
  (the genotype with the highest score is used).
  The method also adjusts scores to be consistent for
  all genotypes from a source, using an algorithm that
  penalizes for high-scoring, disagreeing genotypes.

  Input: -source => $source_char can specify a delimiter
     for which the read name string up to, but not 
     including that character, is considered to be the
     sample source (default '-').
         -source_pos => $position can specify which 
     position the source will be in after splitting on 
     the source character (default 1).
  Output: 1 if object is altered, 0 otherwise.

=cut

###########################################################
sub reconcile_genotypes {
    my $self  = shift;
    my %params = @_;
    my $source_char = $params{-source} || '-';
    my $source_pos = defined($params{-source_pos}) ? $params{-source_pos} : 1;
    my $source_position = $source_pos - 1;
    my $verbose = $params{-verbose};

    my $rh_source_genotypes = {};
    foreach my $genotype (@{$self->genotypes()})
    {
        my @source_parts = split /[$source_char]/, $genotype->read_name();
        my $source = $source_parts[$source_position] || 'None';
        push @{$rh_source_genotypes->{$source}}, $genotype;
    }

    my $rh_new_genotypes = {};
    foreach my $source (keys %{$rh_source_genotypes})
    {
        if ($source eq 'None')
        {
            foreach my $ns_genotype (@{$rh_source_genotypes->{'None'}})
            {
                my $read_name = $ns_genotype->read_name();
                $rh_new_genotypes->{$read_name} = $ns_genotype;
            }
            next;
        }

        # figure out high score genotype:
        my ($high_score, $high_score_read, $best_allele1, $best_allele2);
        foreach my $genotype (@{$rh_source_genotypes->{$source}})
        {
            my $allele1 = $genotype->allele1();
            my $allele2 = $genotype->allele2();
            my $readname = $genotype->read_name();
            my $score = $genotype->score();

            print "Source $source: $allele1/$allele2 score $score\n" if ($verbose);

            if (!defined ($high_score) || ($score > $high_score))
            {
                 $best_allele1 = $allele1;
                 $best_allele2 = $allele2;
                 $high_score = $score;
                 $high_score_read = $readname;
            }
        }

        # adjust the h.s.according to whether other genotypes agree or disagree:

        my $score = $high_score;
        foreach my $genotype (@{$rh_source_genotypes->{$source}})
        {
            next if ($genotype->read_name() eq $high_score_read);
            $score += (($genotype->allele1() eq $best_allele1) &&
                       ($genotype->allele2() eq $best_allele2)) ? 
                      0.25*($genotype->score()) : -1*($genotype->score()); 
        } 

        # now assign all genotypes
        
        foreach my $genotype (@{$rh_source_genotypes->{$source}})
        {
            my $read_name = $genotype->read_name();
            $genotype->allele1($best_allele1);
            $genotype->allele2($best_allele2);
            $genotype->score($score);
            $rh_new_genotypes->{$read_name} = $genotype;
        }
    }
    $self->genotype_hash($rh_new_genotypes);

} ## end reconcile_genotypes

###########################################################

=item db_select()

  This method selects the specified Polymorphism objects from
  the database table tgs.tgs_polymorphism.

  Input: Class name, 
         -dbh => database handle
         -project_id => desired project_id
         -assembly_id => desired assembly id
         -analysis_id => desired analysis id
         -qc_poly_id => QCPoly record Polymorphisms should be associated with
  Output: Reference to an array of Polymorphism objects.

=cut

###########################################################

sub db_select {

    my $class = shift;
    my %params = @_;
    my $project_id = $params{'-project_id'};
    my $assembly_id = $params{'-assembly_id'};
    my $analysis_id = $params{'-analysis_id'};
    my $qc_poly_id = $params{'-qc_poly_id'};
    my $dbh = $params{'-dbh'};

    my $select = qq! SELECT p.polymorphism_id, p.analysis_id, a.assembly_id, 
                            p.left_ref_end, p.right_ref_start,
                            p.score, p.reference_id, pd.label,
                            s.project_id, ad.value
                     FROM   tgs.tgs_poly_dict pd, tgs.tgs_polymorphism p, 
                            tgs.tgs_analysis a, tgs.tgs_assembly s,
                            tgs.tgs_poly_allele pa, tgs.tgs_allele_dict ad
                     WHERE  p.analysis_id = a.analysis_id
                     AND    p.poly_dict_id = pd.poly_dict_id
                     AND    p.polymorphism_id = pa.polymorphism_id
                     AND    pa.allele_id = ad.allele_dict_id
                     AND    a.assembly_id = s.assembly_id !;

    if (defined ($project_id))
    {
        $select .= qq! AND s.project_id = $project_id !;
    }

    if (defined ($assembly_id))
    {
        $select .= qq! AND a.assembly_id = $assembly_id !;
    }

    if (defined ($analysis_id))
    {
        $select .= qq! AND p.analysis_id = $analysis_id !;
    }

    if (defined ($qc_poly_id))
    {
        $select .= qq! AND p.polymorphism_id IN (SELECT
                           poly_id FROM tgs.tgs_poly_qc_poly
                           WHERE poly_qc_id = $qc_poly_id )  !;
    }

    $select .= qq! ORDER BY pa.reference DESC !; # to get reference allele first

    my $sth = $dbh->prepare($select)
         or croak "Couldn\'t prepare statement $select!\n";
    my $rv = $sth->execute()
         or croak "Couldn\'t execute statement $select!\n";

    my %return_polys = (); # will contain Polymorphism's by id
    my $rh_allele_seqs = {}; # to hold refs to lists of aleles for each poly
    while (my $ra_poly = $sth->fetchrow_arrayref())
    {
        my ($polymorphism_id, $analysis_id, $assembly_id,
             $left_ref_end, $right_ref_start, $score, $amplimer_id, 
             $type, $project_id, $allele_seq) = @{$ra_poly};

        if (!$return_polys{$polymorphism_id})
        {
            $return_polys{$polymorphism_id} = $class->new(
                           -polymorphism_id => $polymorphism_id,
                           -left_flank_end => $left_ref_end,
                           -right_flank_start => $right_ref_start,
                           -type => $type,
                           -score => $score,
                           -amplimer_id => $amplimer_id,
                           -assembly_id => $assembly_id,
                           -analysis_id => $analysis_id,
                           -strict_boundaries => 1,
                           -project_id => $project_id, # not currently used
                           -dbh => $dbh);
        }
        push @{$rh_allele_seqs->{$polymorphism_id}}, $allele_seq;
    }

    my @return_pols = ();
    foreach my $id (keys %return_polys)
    {
         $return_polys{$id}->allele_seqs($rh_allele_seqs->{$id});
         push @return_pols, $return_polys{$id};
    }

    return \@return_pols;

} ## end db_select


1;
__END__

=back
