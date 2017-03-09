#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl -w
# $Id: SVcomp.pl 7772 2017-03-06 20:25:29Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use Bio::DB::Sam;
use NHGRI::MUMmer::AlignSet;

our %Opt;

=head1 NAME

SVcomp.pl - compare structural variants in VCF format by constructing alternate haplotypes and comparing them.

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

process_commandline();

my $vcf_file1 = $ARGV[0];
my $vcf_file2 = $ARGV[1];

my $workingdir = $Opt{workdir};
my $ref_fasta = $Opt{ref};
my $fai_obj = Bio::DB::Sam::Fai->load($ref_fasta);

my $vcf1_fh = Open($vcf_file1);
my $vcf2_fh = Open($vcf_file2);

my $varpairindex = 1;
while (<$vcf1_fh>) {
    next if (/^#/);

    my $vcf1_line = $_;
    my ($chr1, $start1, $end1, $id1, $ref1, $alt1) = parse_vcf_line($vcf1_line);

    my $vcf2_line = <$vcf2_fh>;
    while ($vcf2_line =~ /^#/) {
        $vcf2_line = <$vcf2_fh>;
    }

    my ($chr2, $start2, $end2, $id2, $ref2, $alt2) = parse_vcf_line($vcf2_line);

    # boundaries of constructed haplotype:

    my $left_bound = ($start1 < $start2) ? $start1 - $Opt{buffer} : $start2 - $Opt{buffer};
    my $right_bound = ($end1 < $end2) ? $end2 + $Opt{buffer} : $end1 + $Opt{buffer};
    my $alt_hap1 = construct_alt_hap($fai_obj, $chr1, $left_bound, $right_bound, $start1, $end1, \$ref1, \$alt1);
    my $alt_hap2 = construct_alt_hap($fai_obj, $chr2, $left_bound, $right_bound, $start2, $end2, \$ref2, \$alt2);

    # align the alternative haplotypes to each other and evaluate

    compare_haplotypes(\$alt_hap1, \$alt_hap2, $varpairindex, $id1, $id2);
    $varpairindex++;

}

close $vcf1_fh;
close $vcf2_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( buffer => 100, workdir => '.' );
    GetOptions(\%Opt, qw( ref=s ignore_length check_ref buffer=s workdir=s cleanup manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVcomp.pl, ", q$Revision: 7772 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    pod2usage("SVcomp.pl --ref first_vcf.vcf second_vcf.vcf") if !@ARGV;

    if ($Opt{workdir} ne '.') {
        mkdir $Opt{workdir}; # don't stress if it was already there, so not checking return value
    }
}

sub parse_vcf_line {
    my $vcf_line = shift;

    #print "$vcf_line\n";
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

sub construct_alt_hap {
    my $fai_obj = shift;
    my $chr = shift;
    my $left = shift;
    my $right = shift;
    my $start = shift;
    my $end = shift;
    my $r_ref = shift;
    my $r_alt = shift;

    #my $alt_allele = uc($gtb_ref->seq("$chr:$left-$right"));
    #print "Extracting sequence $chr:$left-$right\n";
    my $alt_allele = uc($fai_obj->fetch("$chr:$left-$right"));
    my $offset = $start - $left;
    my $length = $end - $start + 1;
    #print "Subseting offset $offset, length $length\n";
    my $extracted_ref = substr($alt_allele, $offset, $length); 
    if (($Opt{check_ref}) && ($extracted_ref ne ${$r_ref})) {
        die "Extracted reference for $chr:$start-$end does not match provided REF!\n";
    }

    #print "Replacing ref at offset $offset, length $length with ${$r_alt}\n";
    substr($alt_allele, $offset, $length) = ${$r_alt}; 

    return $alt_allele;
}

sub compare_haplotypes {
    my $r_alt1 = shift;
    my $r_alt2 = shift;
    my $pair_id = shift;
    my $id1 = shift;
    my $id2 = shift;

    my $length1 = length(${$r_alt1});
    my $length2 = length(${$r_alt2});

    my $seq1_fh = Open("$workingdir/alt1.$pair_id.fa", "w"); 
    print $seq1_fh ">Pair$pair_id.alt1\n";
    my $r_alt1_50 = format_50($r_alt1);
    print $seq1_fh ${$r_alt1_50};
    if (${$r_alt1_50} !~ /\n$/) {
        print $seq1_fh "\n";
    }
    close $seq1_fh;

    my $seq2_fh = Open("$workingdir/alt2.$pair_id.fa", "w"); 
    print $seq2_fh ">Pair$pair_id.alt2\n";
    my $r_alt2_50 = format_50($r_alt2);
    print $seq2_fh ${$r_alt2_50};
    if (${$r_alt2_50} !~ /\n$/) {
        print $seq2_fh "\n";
    }
    close $seq2_fh;

    my $delta_file = run_mummer("$workingdir/alt1.$pair_id.fa", "$workingdir/alt2.$pair_id.fa", "pair$pair_id");
    my $mummer_obj = NHGRI::MUMmer::AlignSet->new(-delta_file => $delta_file);
    my $ra_entrypairs = $mummer_obj->{entry_pairs} || [];
    my $no_pairs = @{$ra_entrypairs};
    if ($no_pairs != 1) {
        print "NOMATCH: pair $pair_id $no_pairs entry pairs $id1 $id2 $length1 $length2\n";
    }
    else {
        my $rh_pairentry = $ra_entrypairs->[0];
        my $match_found = 0;
        foreach my $rh_align (@{$rh_pairentry->{aligns}}) {
            my $start1 = $rh_align->{ref_start};
            my $end1 = $rh_align->{ref_end};
            my $start2 = $rh_align->{query_start};
            my $end2 = $rh_align->{query_end};
            my $cigar_string = $rh_align->{cigar_string};
            my $mismatches = $rh_align->{mismatches};
            if ($start1 < 5 && $start2 < 5 && $length1 - $end1 < 5 and $length2 - $end2 < 5) {
                print "Pair$pair_id $id1 $id2 $start1-$end1 (length $length1) matches $start2-$end2 (length $length2) $mismatches mismatches $cigar_string\n";
                $match_found = 1;
            }
            else {
                #print "pair $pair_id NONGLOBAL\n";
            }
        }
        if (!$match_found) {
            print "NOMATCH: pair $pair_id $no_pairs entry pairs $id1 $id2 $length1 $length2\n";
        }
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

sub run_mummer {
    my $fasta1 = shift;
    my $fasta2 = shift;
    my $nucmer_name = shift;

    my $mummer_cmd_file = "$workingdir/run_mummer.$nucmer_name.sh";
    my $mummer_fh = Open("$mummer_cmd_file", "w");
    print $mummer_fh "#!/bin/bash\n/home/nhansen/projects/MUMmer/MUMmer3.23/nucmer -o -p $workingdir/nucmer.$nucmer_name -maxmatch $fasta1 $fasta2\n";
    close $mummer_fh;
    chmod 0755, $mummer_cmd_file;
    
    my $cmd = "$mummer_cmd_file > $workingdir/nucmer.$nucmer_name.out 2> $workingdir/nucmer.$nucmer_name.err";
    system($cmd) == 0
        or print STDERR "Something went wrong running $mummer_cmd_file!\n";

    if ($Opt{'cleanup'}) {
        unlink $fasta1;
        unlink $fasta2;
        unlink "$workingdir/nucmer.$nucmer_name.err";
        unlink "$workingdir/nucmer.$nucmer_name.coords";
        unlink "$workingdir/nucmer.$nucmer_name.out";
        unlink "$workingdir/run_mummer.$nucmer_name.sh";
    }

    return "$workingdir/nucmer.$nucmer_name.delta";
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
