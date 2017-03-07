#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl -w
# $Id$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use NHGRI::MUMmer::AlignSet;

our %Opt;

=head1 NAME

SVbackgenotype.pl - Create a single VCF with genotypes after using SVrefine.pl to report structural variants in multiple assemblies for the same candidate regions.

=head1 SYNOPSIS

  SVbackgenotype.pl --target target_vcffile.vcf --info file_of_svrefine_files.txt --ref reference.fasta

=head1 DESCRIPTION

The program creates a single VCF file for variants contained in the target VCF file specified in the first argument, using the SVrefine.pl output detailed in the "file_of_svrefine_files.txt" file.  For each SV contained in the target VCF, SVbackgenotype.pl examines output in the reference coverage and SV VCF files specified in the SVrefine info file.

=cut

#------------
# Begin MAIN 
#------------

$|=1;

process_commandline();

my $target_vcffile = $Opt{target};
my $svrefine_infofile = $Opt{info};
my $ref_fasta = $Opt{ref};
my $fai_obj = Bio::DB::Sam::Fai->load($ref_fasta);

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( );
    GetOptions(\%Opt, qw( ref=s target=s info=s manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVbackgenotype.pl, ", q$Revision:$, "\n"; }
    # If non-option arguments are required, uncomment next line
    pod2usage("SVbackgenotype.pl --ref ref.fasta --target target_vcf.vcf --samples svrefine_info.txt") if !@ARGV;

}

sub parse_vcf_line {
    my $vcf_line = shift;

    if ($vcf_line =~ /^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        my ($chr, $start, $id, $ref, $alt, $info) = ($1, $2, $3, $4, $5, $8);
        $ref = uc($ref);
        $alt = uc($alt);

        my $end; # if we have sequence, end will be determined as last base of REF, otherwise, from END=
        if (($ref =~ /^([ATGC]+)$/) && ($alt =~ /^([ATGC]+)$/)) {
            $end = $start + length($ref) - 1;
        }
        else { # check for END= INFO tag
            if ($info =~ /END=\s*(\d+)/) {
                $end = $1;
            }
            else {
                die "Variants without END in INFO field must have sequence alleles in REF and ALT fields!\n";
            }
        }

        if (!$Opt{'ignore_length'} && ($ref =~ /^([ATGC]+)$/) && ($end - $start + 1 != length($ref))) {
            die "Length of reference allele does not match provided POS, END!  Use --ignore_length option to ignore this discrepancy.\n";
        }

        return ($chr, $start, $end, $id, $ref, $alt);
    }
    else {
        die "Unexpected VCF line:\n$vcf_line";
    }
}

# format a sequence string for writing to FASTA file
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
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=head1 AUTHOR

 Nancy Hansen - nhansen@mail.nih.gov

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
