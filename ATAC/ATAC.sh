#!/bin/bash

set -euo pipefail

SAMPLE_SHEET="config/samples.tsv"

declare -A seen_groups

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r sample_id group stage tissue genotype rep ref; do
    echo "Submitting align for ${sample_id} (${stage}, ${ref})"
    jid_align=$(sbatch --parsable scripts/01_align.sh "$sample_id" "$stage" "$ref")

    echo "Submitting clean for ${sample_id} after align job ${jid_align}"
    jid_clean=$(sbatch --parsable --dependency=afterok:${jid_align} scripts/02_clean.sh "$sample_id" "$stage")

    if [[ -z "${seen_groups[$group]+x}" ]]; then
        seen_groups[$group]="$jid_clean"
    else
        seen_groups[$group]="${seen_groups[$group]}:$jid_clean"
    fi
done

for group in "${!seen_groups[@]}"; do
    stage=$(echo "$group" | grep -o 'E[0-9]\+')

    echo "Submitting peak calling for ${group} (${stage}) after clean jobs ${seen_groups[$group]}"
    sbatch --dependency=afterok:${seen_groups[$group]} scripts/03_callpeak.sh "$group" "$stage"
done
