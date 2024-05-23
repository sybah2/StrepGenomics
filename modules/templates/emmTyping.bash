#!/usr/bin/env bash

set -euo pipefail
TMP=${params.emmdb}
cp -r \$TMP GSQ_DB
DB="\$PWD/GSQ_DB"

As=${assembly}
cp  \$As assembly.fa
df="\$PWD/assembly.fa"

perl ${projectDir}/data/GAS_Scripts_Reference/my_version_emm_Typer4.pl -z "\${df}" -r "\${DB}" -o ${sample_id} -n ${sample_id}
