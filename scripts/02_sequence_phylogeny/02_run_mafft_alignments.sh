#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MANIFEST=""
LOG_FILE=""
THREADS="1"
MODE="auto"
OVERWRITE="false"

usage() {
  cat <<'EOF'
Usage:
  02_run_mafft_alignments.sh \
    --input-dir <fasta_dir> \
    --output-dir <alignment_dir> \
    --manifest <manifest_csv> \
    --log <log_file> \
    [--threads <n>] \
    [--mode auto|localpair|globalpair] \
    [--overwrite]

Description:
  Runs MAFFT on all FASTA files in the input directory and writes aligned FASTA files.

Outputs:
  - aligned FASTA files
  - manifest CSV
  - log file
EOF
}

log() {
  local msg="$1"
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  echo "[$ts] $msg"
  if [[ -n "${LOG_FILE}" ]]; then
    mkdir -p "$(dirname "$LOG_FILE")"
    echo "[$ts] $msg" >> "$LOG_FILE"
  fi
}

count_sequences() {
  local fasta="$1"
  grep -c '^>' "$fasta" || true
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --manifest)
      MANIFEST="$2"
      shift 2
      ;;
    --log)
      LOG_FILE="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --mode)
      MODE="$2"
      shift 2
      ;;
    --overwrite)
      OVERWRITE="true"
      shift 1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$MANIFEST" ]]; then
  echo "ERROR: --input-dir, --output-dir and --manifest are required." >&2
  usage
  exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: MAFFT not found in PATH." >&2
  echo "Install MAFFT or activate the environment where it is available." >&2
  exit 1
fi

case "$MODE" in
  auto|localpair|globalpair)
    ;;
  *)
    echo "ERROR: --mode must be one of: auto, localpair, globalpair" >&2
    exit 1
    ;;
esac

mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname "$MANIFEST")"

if [[ -n "$LOG_FILE" ]]; then
  mkdir -p "$(dirname "$LOG_FILE")"
  : > "$LOG_FILE"
fi

echo "input_fasta,output_alignment,n_sequences,mode,threads,status,error" > "$MANIFEST"

log "Starting MAFFT alignments"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"
log "Manifest: $MANIFEST"
log "Mode: $MODE"
log "Threads: $THREADS"

shopt -s nullglob

FASTA_FILES=("$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.faa)

if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  log "No FASTA files found in $INPUT_DIR"
  exit 1
fi

for fasta in "${FASTA_FILES[@]}"; do
  base="$(basename "$fasta")"
  stem="${base%.*}"
  out="$OUTPUT_DIR/${stem}.aligned.fasta"
  tmp="$out.tmp"
  nseq="$(count_sequences "$fasta")"

  status=""
  error=""

  if [[ "$nseq" -eq 0 ]]; then
    status="skipped_empty"
    error="no_sequences"
    log "$base: skipped because no sequences were found"
    echo "\"$fasta\",\"$out\",$nseq,$MODE,$THREADS,$status,\"$error\"" >> "$MANIFEST"
    continue
  fi

  if [[ -s "$out" && "$OVERWRITE" != "true" ]]; then
    status="skipped_existing"
    log "$base: skipped_existing"
    echo "\"$fasta\",\"$out\",$nseq,$MODE,$THREADS,$status,\"$error\"" >> "$MANIFEST"
    continue
  fi

  log "$base: running MAFFT on $nseq sequences"

  set +e

  if [[ "$MODE" == "auto" ]]; then
    mafft --auto --thread "$THREADS" "$fasta" > "$tmp" 2>> "$LOG_FILE"
  elif [[ "$MODE" == "localpair" ]]; then
    mafft --localpair --maxiterate 1000 --thread "$THREADS" "$fasta" > "$tmp" 2>> "$LOG_FILE"
  elif [[ "$MODE" == "globalpair" ]]; then
    mafft --globalpair --maxiterate 1000 --thread "$THREADS" "$fasta" > "$tmp" 2>> "$LOG_FILE"
  fi

  rc=$?
  set -e

  if [[ "$rc" -eq 0 && -s "$tmp" ]]; then
    mv "$tmp" "$out"
    status="aligned"
    log "$base: aligned"
  else
    rm -f "$tmp"
    status="failed"
    error="mafft_exit_code_${rc}"
    log "$base: failed with exit code $rc"
  fi

  echo "\"$fasta\",\"$out\",$nseq,$MODE,$THREADS,$status,\"$error\"" >> "$MANIFEST"
done

log "Finished MAFFT alignments"
log "Status counts:"
cut -d',' -f6 "$MANIFEST" | tail -n +2 | sort | uniq -c | while read -r count status; do
  log "$status: $count"
done
