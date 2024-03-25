#!/usr/bin/perl -w
# $Id: $

use strict;
use Getopt::Long;
use Pod::Usage;
use FileHandle;

use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr(q$Rev: 0$, 4)/10000;

=head1 NAME

vcf2tdf.pl - script to convert a VCFv4.0 file to a tab-delimited format

=head1 SYNOPSIS

This script reads in a VCF file (minimum v4.0) and outputs tab-delimited format with fields for info tags

=head1 USAGE

vcf2tdf.pl --vcf <vcf file> --reg <optional region> --seqs --noheader

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $vcf_file = $Opt{'vcf'};
my $reg = $Opt{'reg'};
my $seqs = $Opt{'seqs'};

my $tabix_exe = 'tabix';
my $bgzip_exe = 'bgzip';

my $vcf_open = ($Opt{'reg'}) ? "$tabix_exe $vcf_file $Opt{'reg'} | " : (($vcf_file =~ /\.gz$/) ? "$bgzip_exe -d -c $vcf_file | " : $vcf_file);

my $vcf_fh = FileHandle->new("$vcf_open")
    or die "Couldn\'t open $vcf_open for reading: $!\n";

my $ra_info_fields = [];
my @all_prefixes = ();
my @all_rows = ();
my %all_detected_fields = ();
while (<$vcf_fh>) {
    if (/^#/) { # must be first line of header
        $_ = parse_header($vcf_fh, $_, $ra_info_fields);

        if (!$Opt{'noheader'}) { # write a header line
            my $vcf_field_names = ($seqs) ? "#Chrom\tPos\tID\tRef\tAlt\tQual" : "#Chrom\tPos\tID\tRefLength\tAltLength\tQual";
            print "$vcf_field_names";
            my $info_field_names = join "\t", @{$ra_info_fields};
            print "\t$info_field_names" if ($info_field_names);
            print "\n";
        }
    }
    if (/^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        my ($chr, $pos, $var_id, $ref_allele, $alt_allele_string, $score, 
                  $filter, $infotag) =  
                           ($1, $2, $3, $4, $5, $6, $7, $8);

        $ref_allele = ($seqs) ? $ref_allele : length($ref_allele);
        my @alt_alleles = split /,/, $alt_allele_string;
        my @alt_allele_lengths = map { length($_) } @alt_alleles;

        $alt_allele_string = ($seqs) ? $alt_allele_string : join ',', @alt_allele_lengths;

        my @info_pairs = split /;/, $infotag;
        my %info_values = ();
        foreach my $info_field (@info_pairs) {
            if ($info_field =~ /^([^=]+)=([^=]+)$/) {
                $info_values{$1} = $2;
            }
        }

        if (@{$ra_info_fields}) {
            my @all_info_fields = map {defined($info_values{$_}) ? $info_values{$_} : '.'} @{$ra_info_fields};
            my $info_string = join "\t", @all_info_fields;
            print join("\t", $chr, $pos, $var_id, $ref_allele,
                    $alt_allele_string, $score, $info_string), "\n";
        }
        else {
            push @all_prefixes, "$chr\t$pos\t$var_id\t$ref_allele\t$alt_allele_string\t$score";
            push @all_rows, {%info_values};
            foreach my $info_field (keys %info_values) {
                $all_detected_fields{$info_field} = 1;
            }
        }
    }
}

$vcf_fh->close();

if (!(@{$ra_info_fields})) {
    my @info_fields = keys %all_detected_fields;
    foreach my $rh_row (@all_rows) {
        my @all_info_values = map {$rh_row->{$_} || '.'} @info_fields;
        my $info_string = join "\t", @all_info_values;
        my $prefix = shift @all_prefixes;
        print join("\t", $prefix, $info_string), "\n";
    }
}

#------------
# End MAIN
#------------

sub process_commandline {
    
    # Set defaults here
    GetOptions(\%Opt, qw(
                manual help+ version 
                verbose vcf=s
                reg=s seqs noheader
                )) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "vcf2tdf.pl, ", q$Revision: $, "\n"; }

    if (!exists ($Opt{vcf})) { die "Specify input vcf filename with option --vcf.\n"; }

}

sub parse_header {
    my $fh = shift;
    my $first_line = shift;
    my $ra_info_fields = shift;

    # recklessly ignoring the header for now
    while ($first_line) {
        if ($first_line !~ /^#/) { # end of header
             return $first_line;
        }
        elsif ($first_line =~ /INFO=\<ID=([^,]+)/) {
            push @{$ra_info_fields}, $1;
        }
        $first_line = <$fh>;
    }

    return undef;
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=over 4

=item B<--vcf>

Specifies the VCF file to read.

=back

=over 4

=item B<--reg>

Specify an optional region to extract from the VCF file.

=back

=over 4

=item B<--seqs>

Write allele sequences, rather than lengths of sequences, to tab-delimited output.

=back

=over 4

=item B<--noheader>

Suppress header line.

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
