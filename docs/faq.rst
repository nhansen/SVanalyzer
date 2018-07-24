
.. _faq:

SVanalyzer FAQ
================


.. contents::
  :local:


What does my VCF file need to contain?
----------------------------------------------------------------------------------------------
    To run SVmerge, SVcomp, SVwiden, or SVbenchmark on a VCF file or files, the file needs
    to contain meaningful values in each record's "REF", "ALT", and "ID" columns. The REF
    and ALT sequences are used to determine the break symmetry associated with a variant,
    so are necessary to determine the similarity of events that are not exactly equivalent.
    The ID values are necessary to distinguish variants from each other.

