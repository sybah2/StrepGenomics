#!/usr/bin/env bash

set -euo pipefail

abricate --db vfdb --quiet ${fasta} > VirulenceFactorDB.txt
abricate --summary VirulenceFactorDB.txt > VirulenceFactorDB_summary.txt


abricate --db resfinder --quiet ${fasta} > resfinder.txt
abricate --summary resfinder.txt > resfinder_summary.txt

abricate --db plasmidfinder --quiet ${fasta} > plasmidfinder.txt
abricate --summary plasmidfinder.txt > plasmidfinder_summary.txt

abricate --db ncbi --quiet ${fasta} > ncbi.txt
abricate --summary ncbi.txt > ncbi_summary.txt

abricate --db megares --quiet ${fasta} > megares.txt
abricate --summary megares.txt > megares_summary.txt

abricate --db ecoli_vf --quiet ${fasta} > ecoli_vf.txt
abricate --summary ecoli_vf.txt > ecoli_vf_summary.txt

abricate --db ecoh --quiet ${fasta} > ecoh.txt
abricate --summary ecoh.txt > ecoh_summary.txt

abricate --db card --quiet ${fasta} > card.txt
abricate --summary card.txt > card_summary.txt

abricate --db argannot --quiet ${fasta} > argannot.txt
abricate --summary argannot.txt > argannot_summary.txt
