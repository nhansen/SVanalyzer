#!/bin/bash
nucmer -o -p t/test/nucmer.prefix -maxmatch t/test/ref.fasta t/test/alt.fasta
delta-filter -q t/test/nucmer.prefix.delta > t/test/nucmer.prefix.qdelta
