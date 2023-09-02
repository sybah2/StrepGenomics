#!/usr/bin/env bash

set -euo pipefail

prokka --outdir ${sample_id} --addgenes  --genus ${genus} --species ${spps}  --usegenus --force  --prefix ${sample_id} ${assembly}