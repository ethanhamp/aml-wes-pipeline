#!/usr/bin/env bash
# prep_star_samplelist.sh
# Run this ONCE on OSC before submitting the array job.
# Expands the BATCH_4402 entry into per-sample rows, writes star_samples_expanded.tsv,
# and prints the correct --array range to use when submitting.
#
# Usage:
#   bash ~/prep_star_samplelist.sh

SAMPLE_LIST="$HOME/star_samples.tsv"
OUT="$HOME/star_samples_expanded.tsv"

echo "sample_id	run	fastq_path" > "$OUT"

while IFS=$'\t' read -r sample_id run fastq_path; do
    [[ "$sample_id" == "sample_id" ]] && continue  # skip header

    if [[ "$sample_id" == "BATCH_4402" ]]; then
        # Enumerate actual per-sample subdirectories
        while IFS= read -r -d '' subdir; do
            sid=$(basename "$subdir")
            echo -e "${sid}\t${run}\t${subdir}" >> "$OUT"
        done < <(find "$fastq_path" -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null)
    else
        echo -e "${sample_id}\t${run}\t${fastq_path}" >> "$OUT"
    fi
done < "$SAMPLE_LIST"

N=$(( $(wc -l < "$OUT") - 1 ))  # subtract header
echo ""
echo "=== Expanded sample list written to: $OUT ==="
echo "=== Total samples: $N ==="
echo ""
echo "Submit with:"
echo "  sbatch --array=1-${N}%20 ~/star_align_array.sh"
