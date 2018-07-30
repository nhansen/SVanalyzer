package NHGRI::SVanalyzer::Align2SV;
############################################################

=head1 NAME

Align2SV.pm - A Perl module to derive SV coordinates and
          tags from broken alignments

=head1 DESCRIPTION

    Get SVs ready to print in VCF format by looking at
    overlap of alignments, inferring repetitive structure.

=head1 DATE

 October, 2017

=head1 AUTHORS

 Nancy Fisher Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=cut

############################################################
use strict;
use warnings;

use Carp;
use GTB::File qw(Open);
use File::Temp qw/ tempdir /;

our $VERSION  = '0.01';

###########################################################

=over 4

=item B<new()>

  This method creates a Align2SV object from a set of
  alignment coordinates and other variables.

  Input:  -delta_obj - NHGRI::MUMmer::AlignSet object 
          containing alignments of alternate to reference
          sequences.
          -left_align - alignment hash from delta file
          -right_align - alignment hash from delta file

  Output: New Align2SV object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $delta_obj = $params{-delta_obj};
    my $rh_left_align = $params{-left_align};
    my $rh_right_align = $params{-right_align};

    my $self = { delta_obj => $delta_obj,
                 left_align => $rh_left_align,
                 right_align => $rh_right_align,
                  };

    if ((!$rh_left_align) || (!$rh_right_align)) {
        die "Two references to hashes of alignment info must be passed as -left_align and -right_align to the NHGRI::SVanalyzer::Align2SV constructor!\n";
    }

    if ((!$delta_obj) || (ref($delta_obj) ne 'NHGRI::MUMmer::AlignSet')) {
        die "A NHGRI::MUMmer::AlignSet object must be passed to NHGRI::SVanalyzer::Align2SV constructor as -delta_obj!\n";
    }

    bless $self, $class;

    $self->_calculate_type(); # populates "type" field

    return $self;
}

###########################################################

=item B<widen_insertion()>

  This method calculates widened endpoints for SVs of type
  "SIMPLEINS" or "DUP"
  
  Input:  Object.
  Output: 0 if successful, dies otherwise

=cut

###########################################################
sub widen_insertion {
    my $self  = shift;

    my $left_align = $self->{left_align};
    my $right_align = $self->{right_align};
    my $comp = $left_align->{comp};
    my $ref1 = $self->{ref1};
    my $ref2 = $self->{ref2};
    my $query1 = $self->{query1};
    my $query2 = $self->{query2};

    my $chrom = $left_align->{ref_entry};
    $self->{delta_obj}->find_query_coords_from_ref_coord($ref1, $chrom);
    $self->{delta_obj}->find_query_coords_from_ref_coord($ref2, $chrom);

    if (!$left_align->{query_matches}->{$ref1} || $query1 != $left_align->{query_matches}->{$ref1}) {
        die "QUERY1 $query1 doesn\'t match $left_align->{query_matches}->{$ref1}\n";
    }

    if (!$right_align->{query_matches}->{$ref2} || $query2 != $right_align->{query_matches}->{$ref2}) {
        die "QUERY2 $query2 doesn\'t match $right_align->{query_matches}->{$ref2}\n";
    }

    $self->{query1p} = $right_align->{query_matches}->{$ref1};
    $self->{query2p} = $left_align->{query_matches}->{$ref2};

    #print STDERR "query1p $query1p corresponds to ref $ref1 in second align, query2p $query2p corresponds to ref2 $ref2 in first align\n" if ($Opt{verbose});

    if (!defined($self->{query1p})) {
        #print STDERR "NONREPETITIVE INSERTION--need to check\n";
        $self->{query1p} = $query2 - 1;
    }
    if (!defined($self->{query2p})) {
        #print STDERR "NONREPETITIVE INSERTION--need to check\n";
        $self->{query2p} = $query1 - 1;
    }

    return 0;

} ## end widen_insertion

###########################################################

=item B<widen_deletion()>

  This method calculates widened endpoints for SVs of type
  "SIMPLEDEL" or "CONTRAC"
  
  Input:  Object.
  Output: 0 if successful, dies otherwise

=cut

###########################################################
sub widen_deletion {
    my $self  = shift;

    my $left_align = $self->{left_align};
    my $right_align = $self->{right_align};
    my $comp = $left_align->{comp};
    my $ref1 = $self->{ref1};
    my $ref2 = $self->{ref2};
    my $query1 = $self->{query1};
    my $query2 = $self->{query2};

    my $contig = $left_align->{query_entry};
    $self->{delta_obj}->find_ref_coords_from_query_coord($query1, $contig);
    $self->{delta_obj}->find_ref_coords_from_query_coord($query2, $contig);

    if (!$left_align->{ref_matches}->{$query1} || $ref1 != $left_align->{ref_matches}->{$query1}) {
        die "REF1 $ref1 doesn\'t match $left_align->{ref_matches}->{$query1}\n";
    }

    if (!$right_align->{ref_matches}->{$query2} || $ref2 != $right_align->{ref_matches}->{$query2}) {
        die "REF2 $ref2 doesn\'t match $right_align->{ref_matches}->{$query2}\n";
    }

    $self->{ref1p} = $right_align->{ref_matches}->{$query1};
    $self->{ref2p} = $left_align->{ref_matches}->{$query2};

    if (!defined($self->{ref1p})) {
        #print STDERR "NONREPETITIVE DELETION--need to check\n";
        $self->{ref1p} = $ref2 - 1;
    }
    if (!defined($self->{ref2p})) {
        #print STDERR "NONREPETITIVE DELETION--need to check\n";
        $self->{ref2p} = $ref1 + 1;
    }

    return 0;

} ## end widen_deletion

###########################################################

=item B<_calculate_type()>

  This method uses the coordinates of the left and right
  alignment to determine the type of the SV
  
  Input:  Object.
  Output: Value of type, which will be set as a property
          of the object.

=cut

###########################################################
sub _calculate_type {
    my $self  = shift;

    my $left_align = $self->{left_align};
    my $right_align = $self->{right_align};
    my $comp = $left_align->{comp};

    my $ref1 = $left_align->{ref_end};
    my $ref2 = $right_align->{ref_start};
    my $refjump = $ref2 - $ref1 - 1;
               
    my $query1 = $left_align->{query_end}; 
    my $query2 = $right_align->{query_start}; 
    my $queryjump = ($comp) ? $query1 - $query2 - 1 : $query2 - $query1 - 1; # from show-diff code--number of unaligned bases

    my $svsize = $refjump - $queryjump; # negative for insertions/positive for deletions
 
    my $type = ($refjump < 0 && $queryjump < 0) ? (($svsize < 0) ? 'DUP' : 'CONTRAC') :
               (($svsize < 0 ) ? 'SIMPLEINS' : 'SIMPLEDEL'); # we reverse this later on so that deletions have negative SVLEN

    if (($svsize < 0) && ($refjump > 0)) {
        $type = 'SUBSINS';
    }
    elsif (($svsize > 0) && ($queryjump > 0)) {
        $type = 'SUBSDEL';
    }
    elsif ($svsize == 0) { # is this even possible?--yes, but these look like alignment artifacts
        $type = 'ARTIFACT';
        print STDERR "ref1 $ref1 ref2 $ref2 query1 $query1 query2 $query2\n";
        #exit;
    }
    $svsize = abs($svsize);
    my $repeat_bases = ($type eq 'SIMPLEINS' || $type eq 'DUP') ? -1*$refjump : 
                        (($type eq 'SIMPLEDEL' || $type eq 'CONTRAC') ? -1*$queryjump : 0);

    $self->{type} = $type;
    $self->{ref1} = $ref1;
    $self->{ref2} = $ref2;
    $self->{refjump} = $refjump;
    $self->{query1} = $query1;
    $self->{query2} = $query2;
    $self->{queryjump} = $queryjump;
    $self->{svsize} = $svsize;
    $self->{repeat_bases} = $repeat_bases;

    return $type;

} ## end _calculate_type

###########################################################

1;
__END__

=back
