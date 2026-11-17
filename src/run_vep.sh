#!/usr/bin/env bash
# Annotate every variant of both libraries with Ensembl VEP, run locally in Docker.
#
# Input   data/vep_input.vcf    written by src/derive_variants.py
# Output  data/vep_output.txt   one row per variant per transcript, tab separated

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VEP_DATA="$HOME/vep_data"
IMAGE="ensemblorg/ensembl-vep:latest"
FORKS="${FORKS:-4}"
BUFFER="${BUFFER:-500}"
PLUGIN_DATA="$VEP_DATA/plugin_data"
PLUGINS=()
if [ -f "$PLUGIN_DATA/AlphaMissense_hg38.tsv.gz" ]; then
    PLUGINS+=(--plugin "AlphaMissense,file=/opt/vep/.vep/plugin_data/AlphaMissense_hg38.tsv.gz")
fi
for dbnsfp in "$PLUGIN_DATA"/dbNSFP*.gz; do
    if [ -f "$dbnsfp" ]; then
        PLUGINS+=(--plugin "dbNSFP,/opt/vep/.vep/plugin_data/$(basename "$dbnsfp"),ALL,consequence=ALL")
        break
    fi
done
echo "plugins: ${PLUGINS[*]:-none}"

for required in "$VEP_DATA/homo_sapiens" "$VEP_DATA/GRCh38.fa" "$PROJECT_ROOT/data/vep_input.vcf"; do
    if [ ! -e "$required" ]; then
        echo "missing: $required" >&2
        exit 1
    fi
done

docker run --rm \
    -v "$VEP_DATA":/opt/vep/.vep \
    -v "$PROJECT_ROOT/data":/work \
    "$IMAGE" \
    vep \
    --input_file /work/vep_input.vcf --format vcf \
    --output_file /work/vep_output.txt --tab --force_overwrite \
    --cache --offline --dir_cache /opt/vep/.vep \
    --fasta /opt/vep/.vep/GRCh38.fa \
    --species homo_sapiens --assembly GRCh38 \
    --everything --show_ref_allele --no_stats \
    --fork "$FORKS" --buffer_size "$BUFFER" \
    "${PLUGINS[@]}"

echo "rows written: $(grep -vc '^#' "$PROJECT_ROOT/data/vep_output.txt")"
grep -m1 '^#Uploaded_variation' "$PROJECT_ROOT/data/vep_output.txt" | tr '\t' '\n' | wc -l |
    sed 's/^/columns: /'
