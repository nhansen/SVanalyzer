.. _svcomp:

SVcomp
===============

NAME
    SVcomp.pl - calculate "distances" between structural variants in VCF
    format by constructing their alternate haplotypes and comparing them.

SYNOPSIS
      SVcomp.pl --ref reference.fasta first_vcf.vcf second_vcf.vcf

DESCRIPTION
    The program steps through two VCF files, processing pairs of variants that
    are on the same numbered, non-comment line in each file. For each pair of
    lines, the program reports whether the variants result in haplotypes that
    are alignable with discrepancy rates within specified ranges.

    Note that the VCF files don't need to be sorted, and in fact, depending on
    the exact variants being compared, it's possible that they shouldn't be
    sorted.

OPTIONS
    --help|--manual
        Display documentation. One "--help" gives a brief synopsis, "-h -h"
        shows all options, "--manual" provides complete documentation.

