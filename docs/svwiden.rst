.. _svwiden:

===============
**SVwiden**
===============

SVwiden reads a VCF file and uses MUMmer to determine widened
coordinates for structural variants, adding custom tags to the VCF record.

Usage
------------
::

   svanalyzer widen --ref <reference FASTA file> --variants <VCF-formatted variant file> --prefix <prefix for output files>

Options
------------

==========================     =======================================================================================================
 Option                          Description
==========================     =======================================================================================================
**--help|--manual**               Display documentation.
**--ref**                         The reference FASTA file for the supplied VCF file or files.
**--variants**                    A VCF-formatted file of (possibly equivalent) variants to merge.
**--fof**                         A file of paths to VCF-formatted files to merge.
**--prefix**                      Prefix for output file names (default: "widened")
==========================     =======================================================================================================

