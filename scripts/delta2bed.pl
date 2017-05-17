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

delta2bed.pl - Read MUMmer output and write a BED file with coverage

=head1 SYNOPSIS

  delta2bed.pl --delta <path to delta file of alignments> --out <path to output BED file>

For complete documentation, run C<delta2bed.pl -man>

=head1 DESCRIPTION

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $delta_file = $Opt{delta};

my $outbed = $Opt{out};
my $bed_fh = Open($outbed, "w");

my $delta_obj = NHGRI::MUMmer::AlignSet->new(-delta_file => $delta_file);
write_header($bed_fh) if ($Opt{header});
write_edgelets($delta_obj, $bed_fh);

close $bed_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( );
    GetOptions(\%Opt, qw( delta=s out=s header manual help+ verbose version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "delta2bed.pl, ", q$Revision:$, "\n"; }

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

    print $fh "track name=\"delta2bed\"\n";
}

sub write_edgelets {
    my $delta_obj = shift;
    my $bed_fh = shift;

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
                my $ref_start = $rh_align->{ref_start};
                my $ref_end = $rh_align->{ref_end};
                my $mismatches = $rh_align->{mismatches};
                my $score = sprintf("%3.1f", 100*(1 - $mismatches/($ref_end - $ref_start + 1)));
                push @ref_aligns, [$ref_entry, $ref_start, $ref_end, $query_entry, $query_start, $query_end, $score];

            }
        }

        my @sorted_ref_aligns = sort {$a->[1] <=> $b->[1]} @ref_aligns;
        foreach my $ra_ref_align (@sorted_ref_aligns) {
            my ($ref_entry, $ref_start, $ref_end, $query_entry, $query_start, $query_end, $score) = @{$ra_ref_align};
            my $comp = ($query_start > $query_end) ? 1 : 0;
            my $rsmo = $ref_start - 1;
            my $queryname = "$query_entry:$query_start-$query_end";
            my $strand = ($comp) ? '+' : '-';
            print $bed_fh "$ref_entry\t$rsmo\t$ref_end\t$queryname\t$score\t$strand\n";
        }
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--delta <path to delta file>>

Specify a delta file produced by MUMmer with the alignments to be displayed
in BED format.  (Required).

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

=item B<--out <path to which to write a new BED-formatted file>>

Specify the path to which to write a new BED file containing the edgelet
alignments from this delta file.  BEWARE: if this file already exists, it 
will be overwritten!

=item B<--header>

Flag option to write a BED file header (currently just identifying program name).

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
