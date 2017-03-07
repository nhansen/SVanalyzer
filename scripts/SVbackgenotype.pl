#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl -w
# $Id$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use GTB::FASTA;
use NHGRI::MUMmer::AlignSet;

our %Opt;

=head1 NAME

SVbackgenotype.pl - Create a single VCF with genotypes after using SVrefine.pl to report structural variants in multiple assemblies for the same candidate regions.

=head1 SYNOPSIS

  SVbackgenotype.pl --target target_samplename --info file_of_svrefine_files.txt --ref reference.fasta

=head1 DESCRIPTION

The program creates a single VCF file for variants contained in the VCF file for the target sample specified with the --target option, using the SVrefine.pl output detailed in the "file_of_svrefine_files.txt" file specified with the --info option.  For each SV contained in the target VCF, SVbackgenotype.pl examines output in the SV VCF and reference region files specified in the SVrefine info file.

=cut

#------------
# Begin MAIN 
#------------

$|=1;

process_commandline();

my $target_sample = $Opt{target};
my $svrefine_infofile = $Opt{info};
my $ref_fasta = $Opt{ref};
my $ref_obj = GTB::FASTA->new($ref_fasta);

# Read in info file:
my ($ra_sample_order, $rh_sample_info) = read_info_file($svrefine_infofile);

if (!$rh_sample_info->{$target_sample}) {
    die "Info file passed with --info option must contain a line for the target sample!\n";
}

my $out_vcf_fh = Open("out.multi.vcf", "w");
write_multisample_vcf($target_sample, $ra_sample_order, $rh_sample_info, $out_vcf_fh);

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( maxdist => 1000 );
    GetOptions(\%Opt, qw( ref=s target=s info=s maxdist=i manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVbackgenotype.pl, ", q$Revision:$, "\n"; }
    # If non-option arguments are required, uncomment next line
    pod2usage("SVbackgenotype.pl --ref ref.fasta --target target_vcf.vcf --samples svrefine_info.txt") if ($#ARGV >= 0);

}

sub read_info_file {
    my $infofile = shift;

    my $info_fh = Open($infofile);

    my %sample_info = ();
    my @sample_order = ();
    while (<$info_fh>) {
        chomp;
        my ($sample, $sample_vcf, $sample_cov_file) = split /\s/, $_;
        $sample_info{$sample} = {vcffile => $sample_vcf, covfile => $sample_cov_file};
        $sample_info{$sample}->{svs} = read_vcf_variants($sample_vcf);
        my $no_vars = @{$sample_info{$sample}->{svs}};
        print "Read $no_vars variants for $sample\n";
        $sample_info{$sample}->{hr_regions} = read_ref_coverage($sample_cov_file);
        my $no_hr_regions = @{$sample_info{$sample}->{hr_regions}};
        print "Read $no_hr_regions HOMREF regions for $sample\n";
        push @sample_order, $sample;
    }

    return ([@sample_order], {%sample_info});
}

sub read_vcf_variants {
    my $vcffile = shift;

    my $vcf_fh = Open($vcffile);

    my @variants = ();
    while (<$vcf_fh>) {
        next if (/^#/);
        push @variants, parse_vcf_line($_);
    }
    close $vcf_fh;

    return [@variants];
}

sub read_ref_coverage {
    my $covfile = shift;

    my $cov_fh = Open($covfile);

    my @ref_regions = ();
    while (<$cov_fh>) {
        chomp;
        if (/^HOMREF/) {
            my ($homref, $chrom, $start, $end, $hr_chrom, $hr_start, $hr_end, $hr_contig, $rest) = split /\t/, $_;
            if (!grep {$_->{chrom} eq $hr_chrom && $_->{start}==$hr_start && $_->{end}==$hr_end && $_->{contig} eq $hr_contig} @ref_regions) {
                push @ref_regions, {chrom => $hr_chrom, start => $hr_start, end => $hr_end, contig => $hr_contig};
            }
        }
    }
    close $cov_fh;

    return [@ref_regions];
}

sub parse_vcf_line {
    my $vcf_line = shift;

    if ($vcf_line =~ /^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        my ($chrom, $start, $id, $ref, $alt, $info) = ($1, $2, $3, $4, $5, $8);
        $ref = uc($ref);
        $alt = uc($alt);

        my %variant = ();
        $variant{chrom} = $chrom;
        $variant{start} = $start;
        $variant{refseq} = $ref;
        $variant{altseq} = $alt;
        $variant{info} = $info;

        if (($ref =~ /^([ATGC]+)$/) && ($alt =~ /^([ATGC]+)$/)) {
            $variant{end} = $variant{start} + length($variant{refseq}) - 1;
        }
        else { # check for END= INFO tag
            if ($info =~ /END=\s*(\d+)/) {
                $variant{end} = $1;
            }
            else {
                die "Variants without END in INFO field must have sequence alleles in REF and ALT fields!\n";
            }
        }

        if ($info =~ /REFWIDENED=([^:]+):(\d+)\-(\d+)/) {
            $variant{widechr} = $1;
            $variant{widestart} = $2;
            $variant{wideend} = $3;
            $variant{widemean} = ($variant{widestart} + $variant{wideend})/2.0;
        }
        else {
            print STDERR "Skipping variant at $chrom:$start with no REFWIDENED tag in INFO field!\n";
            return;
        }

        if ($info =~ /SVTYPE=([^;]+)/) {
            $variant{svtype} = $1;
        }
        else {
            print STDERR "Skipping variant at $variant{chrom}:$variant{start} with no SVTYPE tag in INFO field!\n";
            return;
        }

        if ($info =~ /SVLEN=([^;]+)/) {
            $variant{svlen} = $1;
        }
        else {
            print STDERR "Skipping variant at $variant{chrom}:$variant{start} with no SVLEN tag in INFO field!\n";
            return;
        }
        if (!$Opt{'ignore_length'} && ($ref =~ /^([ATGC]+)$/) && ($variant{end} - $variant{start} + 1 != length($variant{refseq}))) {
            die "Length of reference allele does not match provided POS, END!  Use --ignore_length option to ignore this discrepancy.\n";
        }

        $variant{vcfline} = $vcf_line;
        return { %variant };
    }
    else {
        die "Unexpected VCF line:\n$vcf_line";
    }
}

sub write_multisample_vcf {
    my $target_sample = shift;
    my $ra_sample_order = shift;
    my $rh_sample_info = shift;
    my $out_vcf_fh = shift;

    foreach my $rh_variant (@{$rh_sample_info->{$target_sample}->{svs}}) {
        my $chrom = $rh_variant->{chrom};
        my $widestart = $rh_variant->{widestart};
        my $wideend = $rh_variant->{wideend};
        my $widemean = $rh_variant->{widemean};
        print "Processing $chrom:$widestart-$wideend\n";
        my $svtype = $rh_variant->{svtype};
        my $svlen = $rh_variant->{svlen};
        my $vcfline = $rh_variant->{vcfline};

        my %genotype = ();
        foreach my $sample (@{$ra_sample_order}) {
            my @alleles = ();

            # any HomRef coverage?
            my @hr_regions = grep {$_->{chrom} eq $chrom && $_->{start} < $widestart && $_->{end} > $wideend} 
                                 @{$rh_sample_info->{$sample}->{hr_regions}};
            my $no_contigs = @hr_regions;
            print "$sample has $no_contigs HR contigs\n";
            foreach my $rh_region (@hr_regions) {
                print "$rh_region->{chrom}\t$rh_region->{start}\t$rh_region->{end}\t$rh_region->{contig}\n";
                push @alleles, '0';
            }

            # what about evidence of the variant itself?
            if ($sample eq $target_sample) {
                push @alleles, '1';
            }
            else {
                my @exact_matches = grep {$_->{chrom} eq $chrom && $_->{widestart} == $widestart && $_->{wideend} == $wideend &&
                                            $_->{svtype} eq $svtype && $_->{svlen} == $svlen } @{$rh_sample_info->{$sample}->{svs}};
                if (@exact_matches) {
                    foreach my $rh_exact_match (@exact_matches) {
                        print "EXACT MATCH $chrom:$widestart-$wideend($svtype size $svlen) vs. $rh_exact_match->{chrom}:$rh_exact_match->{widestart}-$rh_exact_match->{wideend}($rh_exact_match->{svtype} size $rh_exact_match->{svlen})\n";
                        push @alleles, '1';
                    }
                }
                else {
                    my @potential_variants = grep {$_->{chrom} eq $chrom && abs($_->{widemean} - $widemean) < $Opt{maxdist}} @{$rh_sample_info->{$sample}->{svs}};
                    my $no_potential_matches = @potential_variants;
                    print "Sample $sample has $no_potential_matches potential matching variants\n";
                    foreach my $rh_pot_variant (@potential_variants) {
                        print "POTENTIAL MATCH $chrom:$widestart-$wideend($svtype size $svlen) vs. $rh_pot_variant->{chrom}:$rh_pot_variant->{widestart}-$rh_pot_variant->{wideend}($rh_pot_variant->{svtype} size $rh_pot_variant->{svlen})\n";
                        push @alleles, '2';
                    }
                }
            }
            $genotype{$sample} = join '/', @alleles;
            $genotype{$sample} = '.' if ($genotype{$sample} eq '');
        }
        print $out_vcf_fh "$vcfline";
        foreach my $sample (@{$ra_sample_order}) {
            print $out_vcf_fh "\t$genotype{$sample}";
        }
        print $out_vcf_fh "\n";
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
