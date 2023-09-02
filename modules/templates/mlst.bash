#!/usr/bin/env bash

set -euo pipefail

mlst --scheme ${scheme} --novel Novel_alleles --nopath  ${fasta_files} > MLST.txt