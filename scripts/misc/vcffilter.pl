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

vcffilter.pl - script to filter a VCFv4.0 file based on positions and INFO field values

=head1 SYNOPSIS

This script reads in a VCF file (minimum v4.0) and outputs a filtered VCF file

=head1 USAGE

vcffilter.pl --vcf <vcf file> --reg <optional region> --noheader --minsvlen <minimum SV length> --maxsvlen <maximum SV length> --svtype <desired SV type>

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

my $vcf_open = ($Opt{'reg'}) ? "$tabix_exe $vcf_file $Opt{'reg'} | " : "$bgzip_exe -d -c $vcf_file | ";

my $vcf_fh = FileHandle->new("$vcf_open")
    or die "Couldn\'t open $vcf_open for reading: $!\n";

my $ra_info_fields = [];
while (<$vcf_fh>) {
    if (/^#/) { # must be first line of header
        $_ = parse_header($vcf_fh, $_, $ra_info_fields);
    }
    my $vcf_line = $_;
    if ($vcf_line =~ /^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        my ($chr, $pos, $var_id, $ref_allele, $alt_allele_string, $score, 
                  $filter, $infotag) =  
                           ($1, $2, $3, $4, $5, $6, $7, $8);

        my @info_pairs = split /;/, $infotag;
        my %info_values = ();
        foreach my $info_field (@info_pairs) {
            if ($info_field =~ /^([^=]+)=([^=]+)$/) {
                $info_values{$1} = $2;
            }
        }
        next if (defined($Opt{minsvlen}) && defined($info_values{'SVLEN'}) &&
                  ($info_values{'SVLEN'} < $Opt{minsvlen}));

        next if (defined($Opt{maxsvlen}) && defined($info_values{'SVLEN'}) &&
                  ($info_values{'SVLEN'} > $Opt{maxsvlen}));

        next if (defined($Opt{svtype}) && defined($info_values{'SVTYPE'}) &&
                  ($info_values{'SVTYPE'} ne $Opt{svtype}));

        print $vcf_line;
    }
}

$vcf_fh->close();


#------------
# End MAIN
#------------

sub process_commandline {
    
    # Set defaults here
    GetOptions(\%Opt, qw(
                manual help+ version 
                verbose vcf=s
                reg=s noheader minsvlen=i 
                maxsvlen=i svtype=s
                )) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "vcffilter.pl, ", q$Revision: $, "\n"; }

    if (!exists ($Opt{vcf})) { die "Specify input vcf filename with option --vcf.\n"; }

}

sub parse_header {
    my $fh = shift;
    my $first_line = shift;
    my $ra_info_fields = shift;

    while ($first_line) {
        if ($first_line !~ /^#/) { # end of header
             return $first_line;
        }
        elsif ($first_line =~ /INFO=\<ID=([^,]+)/) {
            push @{$ra_info_fields}, $1;
            if (!$Opt{'noheader'}) {
                print $first_line;
            }
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

=item B<--noheader>

Suppress VCF header.

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
