#!/usr/bin/perl -w
# $Id:$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);

use NHGRI::SVanalyzer::Comp;

our %Opt;

=head1 NAME

SVfindmatches - construct two VCF files of side-by-side variants for comparison, based on windowing of a specified size, and then compare the corresponding nearby variants with SVcomp.

=head1 SYNOPSIS

  SVfindmatches --ref reference.fasta first_vcf.vcf second_vcf.vcf

=head1 DESCRIPTION

The program steps through two VCF files, writing two new VCF files with corresponding variants from the two original files that are on the same chromosome and within a defined distance of each other (as measured by their position).  Then, run compSV.pl on the resulting pair of files and report matches.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $vcf_file1 = $ARGV[0];
my $vcf_file2 = $ARGV[1];

my $crossvcf_file1 = $vcf_file1;
$crossvcf_file1 =~ s/\.vcf/.cross.vcf/;
$crossvcf_file1 =~ s/.*\///;
my $crossvcf_file2 = $vcf_file2;
$crossvcf_file2 =~ s/\.vcf/.cross.vcf/;
$crossvcf_file2 =~ s/.*\///;

if ($crossvcf_file1 eq $vcf_file1 || $crossvcf_file2 eq $vcf_file2) {
    die "VCF file names must contain the string \"vcf\" for this program to work!\n";
}

my $ref_fasta = $Opt{ref};

my $maxwindow = $Opt{'maxwindow'};

my $output_file = $Opt{'output'} || 'compsv.out';
my $workdir = $Opt{'workdir'};

my $rh_variants1 = read_variants($vcf_file1);
my $rh_variants2 = read_variants($vcf_file2);

my $output_fh = Open($output_file, "w");
print $output_fh "DIST\tID1\tID2\tAVGALTLENGTH\tALTLENGTHDIFF\tAVGSIZE\tSIZEDIFF\tEDITDIST\tMAXSHIFT\tPOSDIFF\tRELSHIFT\tRELSIZEDIFF\tRELDIST\n";

foreach my $chr (sort keys %{$rh_variants1}) {
    my $rh_chrvars1 = $rh_variants1->{$chr};
    my $rh_chrvars2 = $rh_variants2->{$chr} || {};

    foreach my $var1_id (sort {$a cmp $b} keys %{$rh_chrvars1}) {
        my $position1 = $rh_chrvars1->{$var1_id}->{'pos'}; 
        my @nearby_ids2 = grep {abs($rh_chrvars2->{$_}->{'pos'} - $position1) < $maxwindow} keys %{$rh_chrvars2};
        my $vcf1_line = $rh_chrvars1->{$var1_id}->{'vcf_line'};

        foreach my $nearby_id2 (@nearby_ids2) {
            my $vcf2_line = $rh_chrvars2->{$nearby_id2}->{'vcf_line'};
            print $crossvcf1_fh "$vcf1_line";
            print $crossvcf2_fh "$vcf2_line";
        }
    }
}

close $crossvcf1_fh;
close $crossvcf2_fh;

####TODO#####
# Use perl module rather than system call #
#system("/home/nhansen/projects/SVanalyzer/github/SVanalyzer/scripts/SVcomp.pl --ref $ref_fasta --cleanup --ignore_length $crossvcf_file1 $crossvcf_file2 --workdir $workdir > $output_file");
#
#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( maxwindow => 100000, workdir => '.' );
    GetOptions(\%Opt, qw( ref=s maxwindow=i output=s workdir=s manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVfindmatches.pl, ", q$Revision:$, "\n"; }
    # If non-option arguments are required, uncomment next line
    pod2usage("SVfindmatches.pl --ref <ref fasta> first_vcf.vcf second_vcf.vcf") if !@ARGV;
}

sub read_variants {
    my $vcf_file = shift;

    my %return_variants = (); # keys chromosome, then id, then hash with fields 'pos' and 'vcf_line'

    my $vcf_fh = Open($vcf_file);
    my $varcount = 1;
    while (<$vcf_fh>) {
        my $vcf_line = $_;
        next if ($vcf_line =~ /^#/);

        if ($vcf_line =~ /^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
            my ($chr, $start, $id, $ref, $alt, $info) = ($1, $2, $3, $4, $5, $8);
            if ($id eq '.') {
                $id = $varcount;
                $varcount++;
            }
            $return_variants{$chr}->{$id} = {'pos' => $start, 'vcf_line' => $vcf_line};
        }
        else {
            die "Unexpected VCF line:\n$vcf_line";
        }
    }
    close $vcf_fh;

    return {%return_variants};
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
