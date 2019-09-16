.. SVanalyzer documentation master file, created by
   sphinx-quickstart on Fri Aug  4 13:33:33 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SVanalyzer
==========

.. toctree::
   :hidden:
   :maxdepth: 2
   :name: mastertoc

   SVbenchmark <svbenchmark>
   SVmerge <svmerge>
   SVcomp <svcomp>
   SVwiden <svwiden>
   SVrefine <svrefine>

`SVanalyzer <http://github.com/nhansen/SVanalyzer>`_ is a software package for the 
analysis of large insertions, deletions, and inversions in DNA. SVanalyzer tools
use repeat-aware methods to refine, compare, and cluster different structural 
variant calls.

Install
========

Conda Installation
-------------------

SVanalyzer can be installed using the conda package manager with the bioconda channel.
For details on setting up conda/bioconda, see the `Bioconda user docs <https://bioconda.github.io/user/install.html>`_.

::

  conda create -n svanalyzer
  conda activate svanalyzer
  conda install svanalyzer

Tarball/Github Clone Installation
-------------------

SVanalyzer can also be installed by downloading a `release tarball <https://github.com/nhansen/SVanalyzer/releases>`_
or cloning the `github repository <https://github.com/nhansen/SVanalyzer>`_:

::

  git clone https://github.com/nhansen/SVanalyzer.git

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
======================

*  :ref:`SVbenchmark    <svbenchmark>` - Compare a set of "test" structural variants in VCF format to a known truth set and report sensitivity and specificity
*  :ref:`SVmerge        <svmerge>` - Merge similar sequence-resolved SVs in VCF format
*  :ref:`SVcomp         <svcomp>` - Compare sequence-resolved SVs to each other
*  :ref:`SVwiden        <svwiden>` - Add tags to a VCF file of sequence-resolved SVs detailing surrounding repetitive genomic context
*  :ref:`SVrefine       <svrefine>` - Call sequence-resolved structural variants (SVs) from assembly consensus


