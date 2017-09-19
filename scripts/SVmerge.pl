#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl -w
# $Id:$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use Bio::DB::Sam;
use NHGRI::SVanalyzer::Comp;

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
print "DIST\tID1\tID2\tAVGALTLENGTH\tALTLENGTHDIFF\tAVGSIZE\tSIZEDIFF\tEDITDIST\tMAXSHIFT\tPOSDIFF\tRELSHIFT\tRELSIZEDIFF\tRELDIST\n";

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
            my $comp_obj = NHGRI::SVanalyzer::Comp->new(-ref_fasta => $ref_fasta, 
                                                        -sv1_info => $rh_neighborhood_sv,
                                                        -sv2_info => {%current_sv},
                                                        -workdir => $workingdir);
            if ($comp_obj->potential_match()) {
                my $rh_distance_metrics = $comp_obj->calc_distance();
                my $edit_dist = $rh_distance_metrics->{'edit_distance'};
                my $max_shift = $rh_distance_metrics->{'max_shift'};
                my $altlength_diff = $rh_distance_metrics->{'altlength_diff'};
                my $altlength_avg = $rh_distance_metrics->{'altlength_avg'};
                my $size_diff = $rh_distance_metrics->{'size_diff'};
                my $size_avg = $rh_distance_metrics->{'size_avg'};
                my $shared_denominator = $rh_distance_metrics->{'shared_denominator'};

                my $pos_diff = abs($rh_neighborhood_sv->{pos} - $current_sv{pos});
                # divide maximum shift by the minimum absolute size of the two variants:
                my $d1 = ($Opt{olddist}) ? abs($max_shift)/(minimum(abs((2*$size_avg - $size_diff)/2.0), abs((2*$size_avg + $size_diff)/2.0)) + 1) :
                           abs($max_shift)/$shared_denominator;
                # divide the size difference of the two indels by the average absolute size of the difference
                my $d2 = ($Opt{olddist}) ? abs($size_diff)/(abs($size_avg) + 1) : abs($size_diff)/$shared_denominator;
                # divide edit distance by the minimum alternate haplotype length:
                my $d3 = ($Opt{olddist}) ? abs($edit_dist)/(minimum((2*$altlength_avg - $altlength_diff)/2.0, (2*$altlength_avg + $altlength_diff)/2.0) + 1) :
                           abs($edit_dist)/$shared_denominator;
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
    GetOptions(\%Opt, qw( ref=s vcf=s check_length check_ref workdir=s max_dist=i olddist nocleanup manual help+ version)) || pod2usage(0);
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
