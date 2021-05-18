.. _svbenchmark:

===============
**SVbenchmark**
===============

SVbenchmark compares a set of "test" structural variants in VCF format to a known
truth set (also in VCF format) and outputs estimates of sensitivity and specificity.

Usage
------------
::

   svanalyzer benchmark --ref <reference FASTA file> --test <VCF-formatted file of variants to test> --truth <VCF-formatted file of true variants>

Options
------------

==========================     =======================================================================================================
 Option                          Description
==========================     =======================================================================================================
**--help|--manual**               Display documentation.
**--ref**                         The reference FASTA file for the supplied VCF file or files (required).
**--test**                        A VCF-formatted file of structural variants to test (required).
**--truth**                       A VCF-formatted file of variants to compare against (required).
**--maxdist**                     Disallow matches if positions of two variants are more than maxdist bases from each other (default 100,000).
**--normshift**                   Disallow matches if alignments between alternate alleles have normalized shift greater than normshift (default 0.2)
**--normsizediff**                Disallow matches if alternate alleles have normalized size difference greater than normsizediff (default 0.2)
**--normdist**                    Disallow matches if alternate alleles have normalized edit distance greater than normdist (default 0.2)
**--minsize**                     Only include true variants of size >= minsize for recall calculation and test variants >= minsize for precision calculation (default 0)
**--prefix**                      Prefix for output file names (default: "benchmark")
==========================     =======================================================================================================

Description
------------

For sequence-specified test and truth structural variants in VCF files (i.e., files with ATGC sequences in the REF and
ALT fields), SVbenchmark aligns constructed alternate haplotypes of each test/truth variant pair separated by no more
than the distance specified by the --maxdist option to determine if the pair
represent two equivalent variants.

In the false positive output VCF file, the program reports all test variants that are not equivalent to any true
variant. In the false negative output VCF file, the program reports all true variants that are not equivalent to
any test variant. The recall rate is reported in the report file as the percentage of true variants that are not
false negatives, and the precision is reported as the percentage of test variants that are not false positives.

As of SVanalyzer v0.33, SVbenchmark will include non-sequence-specified deletions in its comparisons so long as
the ALT field values of the VCF deletion records are "<DEL>" and an END value is include in the INFO field (e.g.,
END=5289355).

