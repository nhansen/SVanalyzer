.. _svrefine:

===============
**SVrefine**
===============

SVrefine reads a delta-formatted file of MUMmer alignments of an assembly
to the reference to call structural variants (or refine variants in chosen
genomic regions) and print them out in VCF format.

Usage
------------
::

      SVrefine.pl --delta <path to delta file of alignments> --regions <path to BED-formatted file of regions> --ref_fasta <path to reference multi-FASTA file> --query_fasta <path to query multi-FASTA file> --outvcf <path to output VCF file> --svregions <path to output BED file of SV regions> --outref <path to bed file of homozygous reference regions> --nocov <path to bed file of regions with no coverage>

Options
------------

==========================     =======================================================================================================
 Option                          Description
==========================     =======================================================================================================
**--help|--manual**               Display documentation.
**--delta**                       Path to a delta file produced by MUMmer with alignments to be used for retrieving SVs.
**--regions**                     Path to a BED file of regions to be investigated for structural variants in the assembly (Optional).
**--ref_fasta**                   Path to a multi-fasta file containing the sequences used as a reference in the MUMmer alignment (Optional).
**--query_fasta**                 Path to a multi-fasta file containing the sequences used as a query in the MUMmer alignment (Optional).
**--outvcf**                      Path to which to write a new VCF-formatted file of structural variants.
**--refname**                     String to include as the reference name in the VCF header.
**--samplename**                  String to include as the sample name in the output VCF file.
**--maxsize**                     Specify an integer for the maximum size of SV to report.
**--noheader**                    Flag option to suppress printout of the VCF header.
**--nocov**                       Path to write a BED file with "no coverage" regions (only used when --regions option is specified).
==========================     =======================================================================================================

