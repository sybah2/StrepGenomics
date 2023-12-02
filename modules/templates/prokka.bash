#!/usr/bin/env bash

set -euo pipefail

prokka --outdir ${sample_id} --addgenes  --genus ${genus} --species ${spps}  --usegenus --force  --prefix ${sample_id} ${assembly}


cp ${sample_id}/${sample_id}.gff ${sample_id}.gff
