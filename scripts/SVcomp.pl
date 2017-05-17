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

    my $minsvsize = (abs($size1) < abs($size2)) ? abs($size1) : abs($size2);

    if (potential_aligning_lengths($size1, $size2)) { # essentially, one SV's size is not more than twice the size of the other
        # boundaries of constructed haplotype:
    
        my $left_bound = ($pos1 < $pos2) ? $pos1 : $pos2;
        my $right_bound = ($end1 < $end2) ? $end2 : $end1;
        my $alt_hap1 = construct_alt_hap($fai_obj, $chr1, $left_bound, $right_bound, $pos1, $end1, \$ref1, \$alt1);
        my $alt_hap2 = construct_alt_hap($fai_obj, $chr2, $left_bound, $right_bound, $pos2, $end2, \$ref2, \$alt2);
        my $minhaplength = (length($alt_hap1) < length($alt_hap2)) ? length($alt_hap1) : length($alt_hap2);
    
        if ($alt_hap1 eq $alt_hap2) { # identical variants
            print "$varpairindex\t$id1\t$id2\tEXACTMATCH\t$chr1\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\n";
        }
        else {
            # align the alternative haplotypes to each other and evaluate
            my ($maxshift, $editdistance) = compare_haplotypes(\$alt_hap1, \$alt_hap2, $varpairindex);
            if (($editdistance/$minhaplength < 0.05) && (abs($maxshift) < $minsvsize)) {
                print "$varpairindex\t$id1\t$id2\tNWMATCH\t$chr1\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\t$minhaplength\t$maxshift\t$editdistance\n";
            }
            else {
                print "$varpairindex\t$id1\t$id2\tNWFAIL\t$chr1\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\t$minhaplength\t$maxshift\t$editdistance\n";
            }
        }
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
    my $extracted_ref = substr($alt_allele, $offset, $length); 
    if (($Opt{check_ref}) && ($extracted_ref ne ${$r_ref})) {
        die "Extracted reference for $chr:$start-$end does not match provided REF!\n";
    }

    substr($alt_allele, $offset, $length) = ${$r_alt}; 

    return $alt_allele;
}

sub compare_haplotypes {
    my $r_alt1 = shift;
    my $r_alt2 = shift;
    my $pair_id = shift;

    my $tmpfasta1 = "$workingdir/alt1.$pair_id.fa";
    write_fasta_file($tmpfasta1, "Pair$pair_id.alt1", $r_alt1);
    my $tmpfasta2 = "$workingdir/alt2.$pair_id.fa";
    write_fasta_file($tmpfasta2, "Pair$pair_id.alt2", $r_alt2);

    my $edlib_aligner = '/home/nhansen/projects/SVanalyzer/edlib/edlib-1.1.2/build/bin/edlib-aligner';
    my $nw_output = `$edlib_aligner $tmpfasta1 $tmpfasta2 -p -f CIG_STD`;   
    if ($nw_output =~ /Cigar:\n(.*)/m) {
        my $cigar_string = $1;
        my $score = ($nw_output =~ /score = (\d+)/)  ? $1 : 'NA';
        my $maxshift = calc_max_shift($cigar_string);
        return ($maxshift, $score);
    }
    else {
        return ('NA', 'NA');
    }
}

sub potential_aligning_lengths {
    my $size1 = shift;
    my $size2 = shift;

    my $minsize = (abs($size1) < abs($size2)) ? abs($size1) : abs($size2);
    if (abs($size1 - $size2) >= $minsize) {
        return 0;
    }
    else {
        return 1;
    }
}

sub write_fasta_file {
    my $fasta_file = shift;
    my $seq_id = shift;
    my $r_sequence = shift;

    my $seq_fh = Open("$fasta_file", "w"); 
    print $seq_fh ">seq_id\n";
    my $r_alt_50 = format_50($r_sequence);
    print $seq_fh ${$r_alt_50};
    if (${$r_alt_50} !~ /\n$/) {
        print $seq_fh "\n";
    }
    close $seq_fh;
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
