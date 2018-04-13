#!/usr/bin/perl -w
# $Id:$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use GTB::FASTA;
use NHGRI::MUMmer::AlignSet;
our %Opt;

our $VERSION=0.1;

=head1 NAME

delta2sam.pl - Read MUMmer output and write a SAM file with edgelets from a delta formatted alignment file

=head1 SYNOPSIS

  delta2sam.pl --delta <path to delta file of alignments> --ref_fasta <path to reference multi-FASTA file> --query_fasta <path to query multi-FASTA file> --out <path to output SAM file>

For complete documentation, run C<delta2sam.pl -man>

=head1 DESCRIPTION

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $query_fasta = $Opt{query_fasta}; # will use values from delta file if not supplied as arguments
my $ref_fasta = $Opt{ref_fasta}; # will use values from delta file if not supplied as arguments

my $delta_file = $Opt{delta};

my $outsam = $Opt{out};
my $sam_fh = Open($outsam, "w");

my $delta_obj = NHGRI::MUMmer::AlignSet->new(-delta_file => $delta_file,
                                             -extend_exact => 1);

# set up assembly FASTA objects:
$ref_fasta = $delta_obj->{reference_file} if (!$ref_fasta);
$query_fasta = $delta_obj->{query_file} if (!$query_fasta);

my $ref_db = GTB::FASTA->new($ref_fasta);
write_header($sam_fh, $ref_db) if (!$Opt{noheader});

if ($Opt{index} && !(-e ("$query_fasta.fai"))) {
    system("samtools faidx $query_fasta") == 0
        or die "Unable to index fasta file $query_fasta: $!\n";
}
elsif (!(-e ("$query_fasta.fai"))) {
    die "Query fasta file $query_fasta must be indexed with samtools--use --index option to preindex with this program.\n";
}

my $query_db = GTB::FASTA->new($query_fasta);

write_edgelets($delta_obj, $ref_db, $query_db, $sam_fh);

close $sam_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( );
    GetOptions(\%Opt, qw( delta=s ref_fasta=s query_fasta=s out=s refname=s noheader index alignlengths manual help+ verbose version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "mummer2vcf_SVs.pl, ", q$Revision:$, "\n"; }

    if (!($Opt{delta})) {
        print STDERR "Must specify a delta file path with --delta option!\n"; 
        pod2usage(0);
    }

    if (!($Opt{out})) {
        print STDERR "Must specify a SAM file path to write to with the --out option!\n"; 
        pod2usage(0);
    }
}

sub write_header {
    my $fh = shift;
    my $ref_db = shift;

    print $fh "\@HD\tVN:1.5\tSO:coordinate\n";

    foreach my $ref_entry (sort $ref_db->ids()) {
        my $chrlen = $ref_db->len($ref_entry);
        print $fh "\@SQ\tSN:$ref_entry\tLN:$chrlen\n";
    }

    print $fh "\@PG\tID:delta2sam.pl\tVN:$VERSION\n";
}

sub write_edgelets {
    my $delta_obj = shift;
    my $ref_db = shift;
    my $query_db = shift;
    my $sam_fh = shift;

    my %ref_entry_hash = map {$_->{ref_entry} => 1} @{$delta_obj->{entry_pairs}};
    my @unique_ref_entries = keys %ref_entry_hash;
    @unique_ref_entries = sort @unique_ref_entries;

    foreach my $ref_entry (@unique_ref_entries) {
        my @ref_entry_pairs = grep {$_->{ref_entry} eq $ref_entry} @{$delta_obj->{entry_pairs}};
        my @ref_aligns = ();
        foreach my $rh_entry_pair (@ref_entry_pairs) {
            my $query_entry = $rh_entry_pair->{query_entry};
            my $query_length = $rh_entry_pair->{query_length};
            foreach my $rh_align (@{$rh_entry_pair->{aligns}}) {
                my $query_start = $rh_align->{query_start};
                my $query_end = $rh_align->{query_end};
                my $comp = ($query_start > $query_end) ? 1 : 0;
                my $query_remaining = ($comp) ? $query_length - $query_start : $query_length - $query_end;
                push @ref_aligns, [$query_entry, $rh_align->{ref_start}, $rh_align->{cigar_string}, $query_start, $query_end, $query_remaining];

                if ($Opt{'alignlengths'}) {
                    my $cigar = $rh_align->{'cigar_string'};
                    print "CIGAR $cigar\n";
                    my $longest_match = 0;
                    my $longest_del = 0;
                    my $longest_ins = 0;
                    while ($cigar && $cigar =~ s/^(\d+[MDI])//) {
                        my $op = $1;
                        #print "$op\n"; 
                        if ($op =~ /(\d+)M/) {
                            my $matchlength = $1;
                            if ($matchlength > $longest_match) {
                                $longest_match = $matchlength;
                            }
                        }
                        elsif ($op =~ /(\d+)I/) {
                            my $inslength = $1;
                            if ($inslength > $longest_ins) {
                                $longest_ins = $inslength;
                            }
                        }
                        elsif ($op =~ /(\d+)D/) {
                            my $dellength = $1;
                            if ($dellength > $longest_del) {
                                $longest_del = $dellength;
                            }
                        }
                    }
                    print "LONGESTMATCH $longest_match\n";
                    print "LONGESTINS $longest_ins\n";
                    print "LONGESTDEL $longest_del\n";
                }
            }
        }

        my @sorted_ref_aligns = sort {$a->[1] <=> $b->[1]} @ref_aligns;
        foreach my $ra_ref_align (@sorted_ref_aligns) {
            my ($query_entry, $ref_start, $cigar, $query_start, $query_end, $query_remaining) = @{$ra_ref_align};
            my $flag = ($query_start > $query_end) ? 16 : 0;
            my $comp = ($query_start > $query_end) ? 1 : 0;
            if ((!$comp && $query_start > 1) || ($comp && $query_end > 1)) { # hard clip start
                my $query_hc = ($comp) ? $query_end - 1 : $query_start - 1;
                $cigar = $query_hc."H".$cigar;
            }
            if ($query_remaining > 0) { # hard clip end
                $cigar = $cigar.$query_remaining."H";
            }
            my $seq = $query_db->seq($query_entry, $query_start, $query_end);
            print $sam_fh "$query_entry\t$flag\t$ref_entry\t$ref_start\t0\t$cigar\t*\t*\t*\t$seq\t*\n";
        }
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--delta <path to delta file>>

Specify a delta file produced by MUMmer with the alignments to be used for
retrieving SV sequence information.  Generally, one would use the same
filtered delta file that was used to create the "diff" file (see below).
(Required).

=item B<--diffs <path to diff file>>

Specify a diff file produced by MUMmer's "show-diff" command run with the
-r option on the delta file passed with the --delta option to this program.
(Required).

=item B<--ref_fasta <path to reference multi-fasta file>>

Specify the path to a multi-fasta file containing the sequences used as 
reference in the MUMmer alignment.  If not specified on the command line,
the script uses the reference path obtained by parsing the delta file's
first line.

=item B<--query_fasta <path to query multi-fasta file>>

Specify the path to a multi-fasta file containing the sequences used as 
the query in the MUMmer alignment.  If not specified on the command line,
the script uses the query path obtained by parsing the delta file's
first line.

=item B<--out <path to which to write a new SAM-formatted file>>

Specify the path to which to write a new SAM file containing the edgelet
alignments from this delta file.  BEWARE: if this file already exists, it 
will be overwritten!

=item B<--refname <string to include as the reference name in the VCF header>>

Specify a string to be written as the reference name in the ##reference line 
of the VCF header.

=item B<--noheader>

Flag option to suppress printout of the SAM header.

=item B<--index>

Flag option to specify that unindexed query fasta files should be indexed with samtools by this program.

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
