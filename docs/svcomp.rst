.. _svcomp:

===============
**SVcomp**
===============

SVcomp calculates "distances" between pairs of structural variants in VCF
format by constructing their alternate haplotypes and aligning them to each other.

Usage
------------
::

   svanalyzer comp --ref reference.fasta --first <first VCF-formatted file> --second <second VCF-formatted file>

Options
------------

==========================     =======================================================================================================
 Option                          Description
==========================     =======================================================================================================
**--help|--manual**               Display documentation.
**--ref**                         The reference FASTA file for the supplied VCF file or files.
**--first**                       A VCF-formatted file of variants to compare
**--second**                      Second VCF-formatted file of variants to compare--must have the same number of variants as the first file
==========================     =======================================================================================================

