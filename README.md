# atlas2-extensions
Extensions to the ATLAS pipeline, version 2 and higher

Copyright Jackson M. Tsuji, 2019

**NOTE: These scripts are still in progress - early development only.

## Overview

[ATLAS](https://github.com/pnnl/atlas) is an open-source and extendable pipeline for 
metagenomic analysis. This repo contains a collection of scripts to apply custom 
extensions on top of the existing ATLAS framework for personal analyses.

## Scripts

Scripts are available in the `scripts` directory and can be run from the Linux command 
line. Try running the script with the `-h` flag to get usage and install instructions.

Scripts contained in this repo include:
- `generate_MAG_table.py`: generates something like an OTU-table (in amplicon sequencing) 
for genome bins output by ATLAS
