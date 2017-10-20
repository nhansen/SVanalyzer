.. _svrefine:

SVrefine
===============

NAME
    SVrefine.pl - Read regions from a BED file and use MUMmer alignments of an
    assembly to the reference to refine structural variants in those regions
    and print them out in VCF format.

SYNOPSIS
      SVrefine.pl --delta <path to delta file of alignments> --regions <path to BED-formatted file of regions> --ref_fasta <path to reference multi-FASTA file> --query_fasta <path to query multi-FASTA file> --outvcf <path to output VCF file> --outref <path to bed file of homozygous reference regions> --nocov <path to bed file of regions with no coverage>

    For complete documentation, run "SVrefine.pl -man"

OPTIONS
    --delta <path to delta file>
        Specify a delta file produced by MUMmer with the alignments to be used
        for retrieving SV sequence information. Generally, one would use the
        same filtered delta file that was used to create the "diff" file (see
        below). (Required).

    --regions <path to a BED file of regions>
        Specify a BED file of regions to be investigated for structural
        variants in the assembly (i.e., the query fasta file). (Required).

    --ref_fasta <path to reference multi-fasta file>
        Specify the path to a multi-fasta file containing the sequences used
        as reference in the MUMmer alignment. If not specified on the command
        line, the script uses the reference path obtained by parsing the delta
        file's first line.

    --query_fasta <path to query multi-fasta file>
        Specify the path to a multi-fasta file containing the sequences used
        as the query in the MUMmer alignment. If not specified on the command
        line, the script uses the query path obtained by parsing the delta

    --outvcf <path to which to write a new VCF-formatted file>
        Specify the path to which to write a new VCF file containing the
        structural variants discovered in this comparison. BEWARE: if this
        file already exists, it will be overwritten!

    --refname <string to include as the reference name in the VCF header>
        Specify a string to be written as the reference name in the
        ##reference line of the VCF header.

    --samplename <string to include as the sample name in the "CHROM" line>
        Specify a string to be written as the sample name in the header
        specifying a genotype column in the VCF line beginning with "CHROM".

    --maxsize <maximum size of SV to report>
        Specify an integer for the maximum size of SV to report.

    --noheader
        Flag option to suppress printout of the VCF header.

    --help|--manual
        Display documentation. One "--help" gives a brief synopsis, "-h -h"
        shows all options, "--manual" provides complete documentation.


