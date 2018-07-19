.. _svmerge:

NAME
===============

svanalyzer merge

    group structural variants from a VCF file by calculating a
    distance matrix, then finding connected components of a graph.

SYNOPSIS
===============

      svanalyzer merge --ref <reference FASTA file> --variants <VCF-formatted variant file>
      svanalyzer merge --ref <reference FASTA file> --fof <file of paths to VCF-formatted variant files>

DESCRIPTION
===============

    The program steps through a set of structural variants, calculating distances to other
    nearby variants by comparing their alternate haplotypes. The program
    then reports cluters of variants, and prints a VCF file of "unique"
    variants.

OPTIONS
===============

    --help|--manual
        Display documentation. One "--help" gives a brief synopsis, "-h -h"
        shows all options, "--manual" provides complete documentation.


