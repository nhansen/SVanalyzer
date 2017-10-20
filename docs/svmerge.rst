.. _svmerge:

NAME
    SVmerge.pl - group structural variants from a VCF file by calculating a
    distance matrix, then finding connected components of a graph.

SYNOPSIS
      SVmerge.pl --ref <reference FASTA file> --vcf <variant file>

DESCRIPTION
    The program steps through a VCF file, calculating distances to other
    variants in the file that are nearby by comparing alternate haplotypes. It
    then reports cluters of variants, and prints a VCF file of unique
    variants.

    The VCF file must be sorted by position within each chromosome.

OPTIONS
    --help|--manual
        Display documentation. One "--help" gives a brief synopsis, "-h -h"
        shows all options, "--manual" provides complete documentation.


