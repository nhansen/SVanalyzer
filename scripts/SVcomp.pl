#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl -w
# $Id: SVcomp.pl 7772 2017-03-06 20:25:29Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use NHGRI::MUMmer::AlignSet;
use NHGRI::SVanalyzer::Comp;

our %Opt;

=head1 NAME

SVcomp.pl - calculate "distances" between structural variants in VCF format by constructing their alternate haplotypes and comparing them.

=head1 SYNOPSIS

  SVcomp.pl --ref reference.fasta first_vcf.vcf second_vcf.vcf

=head1 DESCRIPTION

The program steps through two VCF files, processing pairs of variants that are on the same numbered, non-comment line in each file.  For each pair of lines, the program reports whether the variants result in haplotypes that are alignable with discrepancy rates within specified ranges.

Note that the VCF files don't need to be sorted, and in fact, depending on the exact variants being compared, it's possible that they shouldn't be sorted.

=cut

#------------
# Begin MAIN 
#------------

$|=1;

my $edlib = `which edlib-aligner`;
if (!$edlib) {
    die "Running SVcomp requires that the edlib aligner (http://martinsosic.com/edlib/) be in your Linux path.\n";
}
else {
    chomp $edlib;
}

process_commandline();

my $vcf_file1 = $ARGV[0];
my $vcf_file2 = $ARGV[1];

my $workingdir = $Opt{workdir};
my $ref_fasta = $Opt{ref};

my $vcf1_fh = Open($vcf_file1);
my $vcf2_fh = Open($vcf_file2);

my $varpairindex = 1;
while (<$vcf1_fh>) {
    next if (/^#/);

    my $vcf1_line = $_;
    my ($chr1, $pos1, $end1, $id1, $ref1, $alt1) = parse_vcf_line($vcf1_line);

    my $vcf2_line = <$vcf2_fh>;
    while ($vcf2_line =~ /^#/) {
        $vcf2_line = <$vcf2_fh>;
    }

    my ($chr2, $pos2, $end2, $id2, $ref2, $alt2) = parse_vcf_line($vcf2_line);

    my $reflength1 = length($ref1);
    my $reflength2 = length($ref2);
    my $altlength1 = length($alt1);
    my $altlength2 = length($alt2);

    if ($chr1 ne $chr2) {
        print "$varpairindex\t$id1\t$id2\tDIFFCHROM\t$chr1\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\n";
        $varpairindex++;
        next;
    }

    my $size1 = $reflength1 - $altlength1;
    my $size2 = $reflength2 - $altlength2;
 
    my $rh_sv1 = { chrom => $chr1, pos => $pos1, end => $end1, id => $id1, ref => $ref1, alt => $alt1};
    my $rh_sv2 = { chrom => $chr2, pos => $pos2, end => $end2, id => $id2, ref => $ref2, alt => $alt2};

    my $comp_obj = NHGRI::SVanalyzer::Comp->new(-sv1_info => $rh_sv1, -sv2_info => $rh_sv2, -ref_fasta => $ref_fasta);
    my $minsvsize = (abs($size1) < abs($size2)) ? abs($size1) : abs($size2);

    if ($comp_obj->potential_match()) {
        my $rh_distance_metrics = $comp_obj->calc_distance(); 
        my $edit_distance = $rh_distance_metrics->{edit_distance};
        my $max_shift = $rh_distance_metrics->{max_shift};
        my $minhaplength = $rh_distance_metrics->{minhaplength};
        my $matchtype = $rh_distance_metrics->{matchtype};
    
        print "$varpairindex\t$id1\t$id2\t$matchtype\t$chr1\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\n";
    }
    else {
        print "$varpairindex\t$id1\t$id2\tDIFFLENGTHS\t$chr1\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\n";
    }
    $varpairindex++;
}

close $vcf1_fh;
close $vcf2_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( workdir => '.', match_percent => 90, shift_percent => 50 );
    GetOptions(\%Opt, qw( ref=s ignore_length check_ref workdir=s cleanup manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVcomp.pl, ", q$Revision:$, "\n"; }
    # If non-option arguments are required, uncomment next line
    pod2usage("SVcomp.pl --ref first_vcf.vcf second_vcf.vcf") if !@ARGV;

    if ($Opt{workdir} ne '.') {
        mkdir $Opt{workdir}; # don't stress if it was already there, so not checking return value
    }
}

sub parse_vcf_line {
    my $vcf_line = shift;

    if ($vcf_line =~ /^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        my ($chr, $start, $id, $ref, $alt, $info) = ($1, $2, $3, $4, $5, $8);
        $ref = uc($ref);
        $alt = uc($alt);
        #print "REF: $ref, alt $alt\n";

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

        #print "Parsed $chr:$start-$end ($id, $ref, $alt)\n";
        return ($chr, $start, $end, $id, $ref, $alt);
    }
    else {
        die "Unexpected VCF line:\n$vcf_line";
    }
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
