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
**--ref**                         The reference FASTA file for the supplied VCF file or files.
**--test**                        A VCF-formatted file of structural variants to test.
**--truth**                       A VCF-formatted file of variants to compare against.
==========================     =======================================================================================================

