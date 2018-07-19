.. _svmerge:

===============
**SVmerge**
===============

SVmerge groups structural variants from a VCF file by calculating a
distance matrix, then finding connected components of a graph in 
which the nodes are the variants, and edges exist when the distances
are below the specified maximum values.

The program steps through a set of structural variants, calculating distances to other
nearby variants by comparing their alternate haplotypes. The program
then reports cluters of variants, and prints a VCF file of "unique"
variants.

===============
Usage
===============
::

   svanalyzer merge --ref <reference FASTA file> --variants <VCF-formatted variant file>
   svanalyzer merge --ref <reference FASTA file> --fof <file of paths to VCF-formatted variant files>

===============
Options
===============

==========================     =======================================================================================================
 Option                          Description
==========================     =======================================================================================================
**--help|--manual**               Display documentation.
**--ref**                         The reference FASTA file for the supplied VCF file or files.
**--variants**                    A VCF-formatted file of (possibly equivalent) variants to merge.
**--fof**                         A file of paths to VCF-formatted files to merge.
==========================     =======================================================================================================

