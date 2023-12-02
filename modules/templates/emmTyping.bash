#!/usr/bin/env bash

set -euo pipefail

perl ${projectDir}/bin/my_version_emm_Typer4.pl -z ${assembly} -r ${emmDb} -o ${sample_id} -n ${sample_id}
