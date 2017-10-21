.. SVanalyzer documentation master file, created by
   sphinx-quickstart on Fri Aug  4 13:33:33 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SVanalyzer
==========

.. toctree::
   :hidden:
   :maxdepth: 2

   svrefine
   svcomp
   svmerge
   svwiden

`SVanalyzer <http://github.com/nhansen/SVanalyzer>`_ is a software package for the 
analysis of large insertions, deletions, and inversions in DNA. SVanalyzer tools
use repeat-aware methods to refine, compare, and cluster different structural 
variant calls.

Install
========

SVanalyzer can be installed by downloading a `release tarball <https://github.com/nhansen/SVanalyzer/releases>`_
or cloning the `<github repository https://github.com/nhansen/SVanalyzer>`_::

git clone git://github.com/nhansen/SVanalyzer.git

After unzipping the tarball or cloning the directory, build SVanalyzer:

::

  cd SVanalyzer
  perl Build.PL
  ./Build
  ./Build test
  ./Build install

To install SVanalyzer to an alternate location (e.g., if you do not have root
permissions), call "perl Build.PL --install_base $HOME".

Command documentation
=======

*  :ref:`SVrefine       <svrefine>` - Call sequence-resolved structural variants (SVs) from assembly consensus
*  :ref:`SVcomp         <svcomp>` - Compare sequence-resolved SVs to each other
*  :ref:`SVmerge        <svmerge>` - Merge similar sequence-resolved SVs in VCF format
*  :ref:`SVwiden        <svwiden>` - Add tags to a VCF file of sequence-resolved SVs detailing surrounding repetitive genomic context


