#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl -w
# $Id:$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use Bio::DB::Sam;

our %Opt;

=head1 NAME

SVmerge.pl - group structural variants from a VCF file by calculating a distance matrix, then finding connected components of a graph.

=head1 SYNOPSIS

  SVmerge.pl --ref <reference FASTA file> --vcf <variant file>

=head1 DESCRIPTION

The program steps through a VCF file, calculating distances to other variants in the file that are nearby by comparing alternate haplotypes. It then reports cluters of variants, and prints a VCF file of unique variants.

The VCF file must be sorted by position within each chromosome.

=cut

#------------
# Begin MAIN 
#------------

$|=1;

process_commandline();

my $workingdir = $Opt{workdir}; # good to allow use of a temporary file

my $ref_fasta = $Opt{ref};
my $vcf_file = $Opt{vcf};
my $fai_obj = Bio::DB::Sam::Fai->load($ref_fasta);
my $max_distance = $Opt{max_dist};

my $vcf_fh = Open($vcf_file);

my ($current_chrom, $current_lastpos); # to check for sorting
my $total_svs = 0;
my @current_neighborhood_svs = ();
print "DIST\tID1\tID2\tAVGALTLENGTH\tALTLENGTHDIFF\tAVGSIZE\tSIZEDIFF\tEDITDIST\tMAXSHIFT\tRELSHIFT\tRELSIZEDIFF\tRELDIST\n";
while (<$vcf_fh>) {
    next if (/^#/);

    my $vcf_line = $_;
    my ($chrom, $pos, $end, $id, $ref, $alt, $reflength, $altlength) = parse_vcf_line($vcf_line);

    next if (!$chrom);

    # check sorted:
    check_sort(\$current_chrom, \$current_lastpos, $chrom, $pos);

    $total_svs++;

    my %current_sv = ( chrom => $chrom, pos => $pos, end => $end, id => $id, ref => $ref, alt => $alt, 
                       reflength => $reflength, altlength => $altlength, size => $altlength - $reflength );
    if (!@current_neighborhood_svs) {
        push @current_neighborhood_svs, { %current_sv };
        next;
    }
    else { # compare current SV to all SV's in the neighborhood

        # remove different chromosome or distant SVs:
        while (@current_neighborhood_svs && ($current_neighborhood_svs[0]->{chrom} ne $chrom || $current_neighborhood_svs[0]->{pos} < $pos - $max_distance)) {
            shift @current_neighborhood_svs;
        }

        # if chance of match, compare alternate alleles with edlib:
        for (my $i=0; $i<=$#current_neighborhood_svs; $i++) {
            my $rh_neighborhood_sv = $current_neighborhood_svs[$i];
            if (potential_match($rh_neighborhood_sv, { %current_sv })) {
                my ($edit_dist, $max_shift, $altlength_diff, $altlength_avg, $size_diff, $size_avg) = calc_distance($rh_neighborhood_sv, { %current_sv }, $fai_obj);
                my $pos_diff = abs($rh_neighborhood_sv->{pos} - $current_sv{pos});
                my $d1 = abs($max_shift)/(minimum(abs((2*$size_avg - $size_diff)/2.0), abs((2*$size_avg + $size_diff)/2.0)) + 1);
                my $d2 = abs($size_diff)/(abs($size_avg) + 1);
                my $d3 = abs($edit_dist)/(minimum((2*$altlength_avg - $altlength_diff)/2.0, (2*$altlength_avg + $altlength_diff)/2.0) + 1);
                print "DIST\t$rh_neighborhood_sv->{id}\t$id\t$altlength_avg\t$altlength_diff\t$size_avg\t$size_diff\t$edit_dist\t$max_shift\t$pos_diff\t$d1\t$d2\t$d3\n";
            }
        }

        push @current_neighborhood_svs, { %current_sv };
    }
}

close $vcf_fh;

print "Read $total_svs SVs from file $vcf_file\n";

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( workdir => '.', match_percent => 90, shift_percent => 50, max_dist => 100000 );
    GetOptions(\%Opt, qw( ref=s vcf=s check_length check_ref workdir=s max_dist=i cleanup manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVmerge.pl, ", q$Revision:$, "\n"; }

    if ($Opt{workdir} ne '.') {
        mkdir $Opt{workdir}; # don't stress if it was already there, so not checking return value
    }

    if (!$Opt{ref} || !$Opt{vcf}) {
        pod2usage();
    }
}

sub minimum {
    my $a = shift;
    my $b = shift;

    return ($a < $b) ? $a : $b;
}

sub parse_vcf_line {
    my $vcf_line = shift;

    if ($vcf_line =~ /^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        my ($chr, $start, $id, $ref, $alt, $info) = ($1, $2, $3, $4, $5, $8);
        $ref = uc($ref);
        $alt = uc($alt);
 
        my $reflength = length($ref);
        my $altlength = length($alt);
 
        my $end; # if we have sequence, end will be determined as last base of REF, otherwise, from END=
        if (($ref =~ /^([ATGCN]+)$/) && ($alt =~ /^([ATGCN]+)$/)) {
            $end = $start + $reflength - 1;
        }
        elsif ($info =~ /SVTYPE=BND/ || $info =~ /SVTYPE=UNK/) {
            return;
        }
        else { # check for END= INFO tag
            warn "Skipping non-ATGC ref or alt:$vcf_line";
            return;
        }

        return ($chr, $start, $end, $id, $ref, $alt, $reflength, $altlength);
    }
    else {
        die "Unexpected VCF line:\n$vcf_line";
    }
}

sub check_sort {
    my $rs_lastchrom = shift;
    my $rs_lastpos = shift;
    my $chrom = shift;
    my $pos = shift;

    if (($$rs_lastchrom) && ($chrom eq $$rs_lastchrom) && ($pos < $$rs_lastpos)) {
        die "VCF file is not sorted ($chrom, $$rs_lastpos, $pos)\n";
    }

    $$rs_lastchrom = $chrom;
    $$rs_lastpos = $pos;
}

sub calc_distance {
    my $rh_sv1 = shift;
    my $rh_sv2 = shift;
    my $fai_obj = shift;

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

    # boundaries of constructed haplotype:

    my $left_bound = ($pos1 < $pos2) ? $pos1 : $pos2;
    my $right_bound = ($end1 < $end2) ? $end2 : $end1;
    my $alt_hap1 = construct_alt_hap($fai_obj, $chrom, $left_bound, $right_bound, $pos1, $end1, \$ref1, \$alt1);
    my $alt_hap2 = construct_alt_hap($fai_obj, $chrom, $left_bound, $right_bound, $pos2, $end2, \$ref2, \$alt2);
    my $minhaplength = (length($alt_hap1) < length($alt_hap2)) ? length($alt_hap1) : length($alt_hap2);
    my $althaplength_avg = (length($alt_hap1) + length($alt_hap2))/2.0;
    my $althaplength_diff = length($alt_hap1) - length($alt_hap2);

    if ($alt_hap1 eq $alt_hap2) { # identical variants
        print "$id1\t$id2\tEXACTMATCH\t$chrom\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\n";
        return (0, 0, $althaplength_diff, $althaplength_avg, $size2 - $size1, ($size1 + $size2)/2.0);
    }
    else {
        # align the alternative haplotypes to each other and evaluate
        my ($maxshift, $editdistance) = compare_haplotypes(\$alt_hap1, \$alt_hap2, $id1, $id2);
        if (($editdistance/$minhaplength < 0.05) && (abs($maxshift) < $minsvsize)) {
            print "$id1\t$id2\tNWMATCH\t$chrom\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\t$minhaplength\t$maxshift\t$editdistance\n";
            return ($editdistance, $maxshift, $althaplength_diff, $althaplength_avg, $size2 - $size1, ($size1 + $size2)/2.0);
        }
        else {
            print "$id1\t$id2\tNWFAIL\t$chrom\t$pos1\t$pos2\t$reflength1\t$reflength2\t$altlength1\t$altlength2\t$minhaplength\t$maxshift\t$editdistance\n";
            return ($editdistance, $maxshift, $althaplength_diff, $althaplength_avg, $size2 - $size1, ($size1 + $size2)/2.0);
        }
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
    my $id1 = shift;
    my $id2 = shift;

    my $tmpfasta1 = "$workingdir/$id1.$id2.fa";
    write_fasta_file($tmpfasta1, $id1, $r_alt1);
    my $tmpfasta2 = "$workingdir/$id2.$id1.fa";
    write_fasta_file($tmpfasta2, "id2", $r_alt2);

    my $edlib_aligner = '/home/nhansen/projects/SVanalyzer/edlib/edlib-1.1.2/build/bin/edlib-aligner';
    my $nw_output = `$edlib_aligner $tmpfasta1 $tmpfasta2 -p -f CIG_STD`;   
    unlink $tmpfasta1;
    unlink $tmpfasta2;
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

sub potential_match {
    my $rh_sv1 = shift;
    my $rh_sv2 = shift;

    # First score is the max shift of the alignment divided by the min absolute value SV size
    # Max shift is greater than the difference in alternate alleles, which is the difference in SV sizes

    my $altlength1 = $rh_sv1->{altlength};
    my $altlength2 = $rh_sv2->{altlength};
    my $minshift = abs($altlength1 - $altlength2);

    my $size1 = $rh_sv1->{size};
    my $size2 = $rh_sv2->{size};

    my $minsize = (abs($size1) < abs($size2)) ? abs($size1) : abs($size2);

    if ($minshift >= $minsize) {
        return 0;
    }
    
    my $pos1 = $rh_sv1->{pos};
    my $pos2 = $rh_sv2->{pos};

    if (abs($pos2 - $pos1) > 6.0*$minsize) {
        return 0;
    }

    return 1;
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
