.. _svmerge:

===============
**SVmerge**
===============

SVmerge groups structural variants from a VCF file by calculating a
distance matrix, then finding connected components of a graph in 
which the nodes are the variants and edges exist when the distances
are below the specified maximum values.

The program steps through a set of structural variants, calculating distances to other
nearby variants by comparing their alternate haplotypes. The program
then reports clusters of variants, and prints a VCF file of "unique"
variants, where the variant reported in the VCF record is a randomly-chosen
representative from the largest cluster (or a randomly selected largest
cluster, in the case of a tie among cluster sizes) of exactly matching variants.

Alternatively, a file of previously-calculated distances can be provided
with the --distance_file option, and the clustering can be skipped with the option
--skip_clusters.

Usage
------------
::

   svanalyzer merge --ref <reference FASTA file> --variants <VCF-formatted variant file>
   svanalyzer merge --ref <reference FASTA file> --fof <file of paths to VCF-formatted variant files>

Options
------------

==========================     =======================================================================================================
 Option                          Description
==========================     =======================================================================================================
**--help|--manual**               Display documentation.
**--ref**                         The reference FASTA file for the supplied VCF file or files.
**--variants**                    A VCF-formatted file of (possibly equivalent) variants to merge.
**--fof**                         A file of paths to VCF-formatted files to merge.
==========================     =======================================================================================================

