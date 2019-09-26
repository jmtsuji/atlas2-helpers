# atlas2-helpers
Scripts for post-processing of output from the ATLAS pipeline, version 2 and higher

Copyright Jackson M. Tsuji, 2019

**NOTE: These scripts are still in progress - early development only.

## Overview

[ATLAS](https://github.com/metagenome-atlas/atlas) is an open-source and extendable pipeline for 
metagenomic analysis. This repo contains a collection of helper scripts to perform post-processing 
on ATLAS output

## Scripts

Scripts are available in the `scripts` directory and can be run from the Linux command 
line. Try running the script with the `-h` flag to get usage and install instructions.

Scripts contained in this repo include:
- `generate_MAG_table.py`: generates something like an OTU-table (in amplicon sequencing) 
for genome bins output by ATLAS. Requires Python3 and the pandas package.
