.. _svwiden:

NAME
    SVwiden.pl - Read a VCF file and use MUMmer to determine widened
    coordinates for structural variants, adding custom tags to the VCF record.

SYNOPSIS
      SVwiden.pl --invcf <path to input VCF file> --ref <path to reference multi-FASTA file> --outvcf <path to output VCF file>

    For complete documentation, run "SVwiden.pl -man"

OPTIONS
    --ref <path to reference multi-fasta file>
        Specify the path to the multi-fasta file that serves as a reference
        for the structural variants in the VCF file.

    --outvcf <path to which to write a new VCF-formatted file>
        Specify the path to which to write a new VCF file containing the
        structural variants from the input VCF file, but now with tags
        specifying widened coordinates

    --refname <string to include as the reference name in the VCF
    header>
        Specify a string to be written as the reference name in the
        ##reference line of the VCF header.

    --noheader
        Flag option to suppress printout of the VCF header.

    --help|--manual

