#!/usr/bin/env bash

set -euo pipefail

abricate --db ${database} --quiet ${fasta} > AMR.txt
abricate --summary AMR.txt > abricate_summary.txt