# atlas2-helpers
[![GitHub release](https://img.shields.io/badge/Version-0.2.0-blue.svg)](https://github.com/jmtsuji/atlas2-helpers)

Scripts for post-processing of output from the ATLAS pipeline, approx. version 2.14 and higher

Copyright Jackson M. Tsuji, 2023

## Overview

[ATLAS](https://github.com/metagenome-atlas/atlas) is an open-source and extendable pipeline for 
metagenomic analysis. This repo contains a collection of helper scripts to perform post-processing 
on ATLAS output.

## Scripts

Scripts are available in the `scripts` directory and can be run from the Linux command 
line. Try running the script with the `-h` flag to get usage and install instructions.

Scripts contained in this repo include:
- `generate_MAG_table.py`: generates something like an OTU-table (in amplicon sequencing) 
for genome bins output by ATLAS. Requires Python3 and the pandas package. 
Default settings are optimized for ATLAS 2.14 (approx.) and higher.
NOTE: revert to version 0.1.0 for use with older ATLAS2 versions (e.g., v2.2.0) 
