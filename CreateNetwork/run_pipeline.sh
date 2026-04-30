#!/bin/bash

set -euo pipefail

SAMPLE_SHEET="config/E10_sample.tsv"

declare -A clean_jobs
declare -A group_stage

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r sample_id group stage tissue genotype rep ref; do
    group_stage["$group"]="$stage"

    jid_align=$(sbatch --parsable scripts/01_align.sh "$sample_id" "$stage")
    jid_clean=$(sbatch --parsable --dependency=afterok:${jid_align} scripts/02_clean.sh "$sample_id" "$stage")

    if [[ -z "${clean_jobs[$group]+x}" ]]; then
        clean_jobs["$group"]="$jid_clean"
    else
        clean_jobs["$group"]="${clean_jobs[$group]}:$jid_clean"
    fi
done

for group in "${!clean_jobs[@]}"; do
    stage="${group_stage[$group]}"
    sbatch --dependency=afterok:${clean_jobs[$group]} scripts/03_callpeak.sh "$group" "$stage"
done
