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

NOTE: SVmerge only clusters and merges sequence-specific variants, i.e., structural
variants with ATGCN sequences for their REF and ALT alleles, or deletions with a
valid "END" INFO tag. These variants will be printed as singletons unless the 
--seqspecific option is specified (see below).

Usage
------------
::

   svanalyzer merge --ref <reference FASTA file> --variants <VCF-formatted variant file> --prefix <prefix for output files>
   svanalyzer merge --ref <reference FASTA file> --fof <file of paths to VCF-formatted variant files> --prefix <prefix for output files>

Options
------------

==========================     =======================================================================================================
 Option                          Description
==========================     =======================================================================================================
**--help|--manual**               Display documentation.
**--ref**                         The reference FASTA file for the supplied VCF file or files.
**--variants**                    A VCF-formatted file of (possibly equivalent) variants to merge.
**--fof**                         A file of paths to VCF-formatted files to merge.
**--prefix**                      Prefix for output file names (default "merged")
**--max_dist**                    Maximum distance between pairs of variants to perform comparison for potential merging (default: 2000)
**--reldist**                     Maximum allowable edit distance, normalized by the mean length of larger allele for the two variants, in an alignment used to merge two variants
**--relsizediff**                 Maximum allowable alt allele size difference, normalized by the mean length of larger allele for the two variants, to merge two variants
**--relshift**                    Maximum allowable shift, normalized by the mean length of the larger allele for the two variants, in an alignment used to merge two variants.
**--seqspecific**                 With this option, SVmerge will fail to print out any SV that does not have an ATGCN sequence for REF and ALT in the input VCF files.
==========================     =======================================================================================================

