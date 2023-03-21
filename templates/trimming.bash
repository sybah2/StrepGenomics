#!/usr/bin/env bash

set -euo pipefail

java -jar $projectDir/Trimmomatic-0.39/trimmomatic-0.39.jar  PE ${reads[0]} ${reads[1]} ${sample_id}_paired_1.fq.gz  ${sample_id}_unpaired_1.fq.gz ${sample_id}_paired_2.fq.gz ${sample_id}_unpaired_2.fq.gz \
ILLUMINACLIP:$projectDir/Trimmomatic-0.39/adapters/All_adapter-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20