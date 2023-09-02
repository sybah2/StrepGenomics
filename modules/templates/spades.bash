#!/usr/bin/env bash

set -euo pipefail

spades.py -k 21,33,55,77 --pe1-1 ${reads[0]}  --pe1-2 ${reads[1]} --careful -o ${sample_id}
cp ${sample_id}/contigs.fasta ${sample_id}.fasta

#awk -v name='$sample_id' '/^>/{print ">name." ++i; next}{print}' < ${sample_id}/contigs.fasta > ${sample_id}_contigs.fasta

#awk -v f="$sample_id" '/^>/{print ">"f"." ++i; next}{print}' < ${sample_id}/contigs.fasta > ${sample_id}.fasta