#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MANIFEST=""
LOG_FILE=""
FASTTREE_BIN=""
MODEL="wag"
GAMMA="true"
MIN_SEQUENCES="3"
OVERWRITE="false"

usage() {
  cat <<'EOF'
Usage:
  03_run_fasttree.sh \
    --input-dir <alignment_dir> \
    --output-dir <tree_dir> \
    --manifest <manifest_csv> \
    --log <log_file> \
    [--fasttree-bin FastTree|FastTreeMP] \
    [--model wag|jtt|auto] \
    [--gamma true|false] \
    [--min-sequences <n>] \
    [--overwrite]

Description:
  Runs FastTree on aligned protein FASTA files and writes Newick trees.

Notes:
  - Default model: WAG
  - Default gamma correction: enabled
  - Alignments with fewer than --min-sequences are skipped.
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
    --fasttree-bin)
      FASTTREE_BIN="$2"
      shift 2
      ;;
    --model)
      MODEL="$2"
      shift 2
      ;;
    --gamma)
      GAMMA="$2"
      shift 2
      ;;
    --min-sequences)
      MIN_SEQUENCES="$2"
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

if [[ -z "$FASTTREE_BIN" ]]; then
  if command -v FastTreeMP >/dev/null 2>&1; then
    FASTTREE_BIN="FastTreeMP"
  elif command -v FastTree >/dev/null 2>&1; then
    FASTTREE_BIN="FastTree"
  elif command -v fasttree >/dev/null 2>&1; then
    FASTTREE_BIN="fasttree"
  else
    echo "ERROR: FastTree not found in PATH." >&2
    echo "Install FastTree or activate the environment where it is available." >&2
    exit 1
  fi
fi

case "$MODEL" in
  wag|jtt|auto)
    ;;
  *)
    echo "ERROR: --model must be one of: wag, jtt, auto" >&2
    exit 1
    ;;
esac

case "$GAMMA" in
  true|false)
    ;;
  *)
    echo "ERROR: --gamma must be true or false" >&2
    exit 1
    ;;
esac

mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname "$MANIFEST")"

if [[ -n "$LOG_FILE" ]]; then
  mkdir -p "$(dirname "$LOG_FILE")"
  : > "$LOG_FILE"
fi

echo "input_alignment,output_tree,n_sequences,fasttree_bin,model,gamma,min_sequences,status,error" > "$MANIFEST"

log "Starting FastTree tree inference"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"
log "Manifest: $MANIFEST"
log "FastTree binary: $FASTTREE_BIN"
log "Model: $MODEL"
log "Gamma: $GAMMA"
log "Minimum sequences: $MIN_SEQUENCES"

shopt -s nullglob

ALIGNMENTS=("$INPUT_DIR"/*.aligned.fasta "$INPUT_DIR"/*.aln "$INPUT_DIR"/*.afa)

if [[ ${#ALIGNMENTS[@]} -eq 0 ]]; then
  log "No alignment files found in $INPUT_DIR"
  exit 1
fi

for aln in "${ALIGNMENTS[@]}"; do
  base="$(basename "$aln")"
  stem="${base%.*}"
  stem="${stem%.aligned}"
  out="$OUTPUT_DIR/${stem}.nwk"
  tmp="$out.tmp"
  nseq="$(count_sequences "$aln")"

  status=""
  error=""

  if [[ "$nseq" -lt "$MIN_SEQUENCES" ]]; then
    status="skipped_too_few_sequences"
    error="n_sequences_less_than_minimum"
    log "$base: skipped because it has $nseq sequences"
    echo "\"$aln\",\"$out\",$nseq,$FASTTREE_BIN,$MODEL,$GAMMA,$MIN_SEQUENCES,$status,\"$error\"" >> "$MANIFEST"
    continue
  fi

  if [[ -s "$out" && "$OVERWRITE" != "true" ]]; then
    status="skipped_existing"
    log "$base: skipped_existing"
    echo "\"$aln\",\"$out\",$nseq,$FASTTREE_BIN,$MODEL,$GAMMA,$MIN_SEQUENCES,$status,\"$error\"" >> "$MANIFEST"
    continue
  fi

  log "$base: running FastTree on $nseq sequences"

  cmd=("$FASTTREE_BIN")

  if [[ "$MODEL" == "wag" ]]; then
    cmd+=("-wag")
  elif [[ "$MODEL" == "jtt" ]]; then
    cmd+=("-jtt")
  fi

  if [[ "$GAMMA" == "true" ]]; then
    cmd+=("-gamma")
  fi

  cmd+=("$aln")

  set +e
  "${cmd[@]}" > "$tmp" 2>> "$LOG_FILE"
  rc=$?
  set -e

  if [[ "$rc" -eq 0 && -s "$tmp" ]]; then
    mv "$tmp" "$out"
    status="tree_built"
    log "$base: tree_built"
  else
    rm -f "$tmp"
    status="failed"
    error="fasttree_exit_code_${rc}"
    log "$base: failed with exit code $rc"
  fi

  echo "\"$aln\",\"$out\",$nseq,$FASTTREE_BIN,$MODEL,$GAMMA,$MIN_SEQUENCES,$status,\"$error\"" >> "$MANIFEST"
done

log "Finished FastTree tree inference"
log "Status counts:"
cut -d',' -f8 "$MANIFEST" | tail -n +2 | sort | uniq -c | while read -r count status; do
  log "$status: $count"
done
