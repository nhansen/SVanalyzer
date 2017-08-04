.. SVanalyzer documentation master file, created by
   sphinx-quickstart on Fri Aug  4 13:33:33 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SVanalyzer
==========

.. toctree::
   :hidden:
   :maxdepth: 2

`SVanalyzer <http://github.com/nhansen/SVanalyzer>`_ is a software package for the 
analysis of large insertions, deletions, and inversions in DNA. SVanalyzer tools
use repeat-aware methods to refine, compare, and cluster different structural 
variant calls.

Install
========

SVanalyzer can be installed by downloading a `release tarball <https://github.com/nhansen/SVanalyzer/releases>`_
or cloning the github repository::

git clone git://github.com/nhansen/SVanalyzer.git

After unzipping the tarball or cloning the directory, build SVanalyzer::

cd SVanalyzer
perl Makefile.PL
make
make test
make install

