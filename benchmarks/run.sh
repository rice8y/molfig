#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
manifest="$repo_root/package/typst.toml"
benchmark_source="$repo_root/benchmarks/benchmark.typ"

if [[ "${MOLFIG_BENCH_ENV:-0}" != "1" ]]; then
  echo "Run this benchmark through its locked Nix environment:" >&2
  echo "  nix develop ./benchmarks --command benchmarks/run.sh $*" >&2
  exit 127
fi

export LC_ALL=C

typst_version="$(typst --version)"
hyperfine_version="$(hyperfine --version)"
expected_typst_version="${MOLFIG_BENCH_TYPST_VERSION:-}"

if [[ -z "$expected_typst_version" ]]; then
  echo "The locked benchmark environment did not declare its Typst version." >&2
  exit 1
fi

if [[ "$typst_version" != "typst $expected_typst_version "* ]]; then
  echo "The benchmark requires Typst $expected_typst_version, found: $typst_version" >&2
  exit 1
fi

if [[ "$hyperfine_version" != "hyperfine 1.20.0" ]]; then
  echo "The benchmark requires hyperfine 1.20.0, found: $hyperfine_version" >&2
  exit 1
fi

default_version="$(awk '
  /^\[package\]$/ { in_package = 1; next }
  /^\[/ { in_package = 0 }
  in_package && /^version[[:space:]]*=/ {
    value = $0
    sub(/^[^=]*=[[:space:]]*"/, "", value)
    sub(/"[[:space:]]*$/, "", value)
    print value
    exit
  }
' "$manifest")"

version="${MOLFIG_VERSION:-$default_version}"
namespace="${MOLFIG_NAMESPACE:-preview}"
mode="${MOLFIG_BENCH_MODE:-export}"
warmup="${MOLFIG_BENCH_WARMUP:-1}"
runs="${MOLFIG_BENCH_RUNS:-5}"

all_cases=(
  1crn-bcif-spacefill
  1fyy-cif-surface
  9r1o-pdb-cartoon
  9z4o-pdb-spacefill
  9m1u-pdb-cartoon-auto
  9q12-pdb-cartoon
)
selected_cases=()

usage() {
  cat <<'EOF'
Usage: benchmarks/run.sh [options] [case ...]

Options:
  --version VERSION       Molfig package version (default: package/typst.toml)
  --namespace NAMESPACE   Typst package namespace: preview or local (default: preview)
  --mode MODE             export or render (default: export)
  --warmup COUNT          Hyperfine warmup runs (default: 1)
  --runs COUNT            Hyperfine measured runs (default: 5)
  -h, --help              Show this help

The same settings can be supplied through MOLFIG_VERSION, MOLFIG_NAMESPACE,
MOLFIG_BENCH_MODE, MOLFIG_BENCH_WARMUP, and MOLFIG_BENCH_RUNS.
EOF
}

while (($# > 0)); do
  case "$1" in
    --version)
      version="$2"
      shift 2
      ;;
    --namespace)
      namespace="$2"
      shift 2
      ;;
    --mode)
      mode="$2"
      shift 2
      ;;
    --warmup)
      warmup="$2"
      shift 2
      ;;
    --runs)
      runs="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      selected_cases+=("$1")
      shift
      ;;
  esac
done

if [[ "$namespace" != "preview" && "$namespace" != "local" ]]; then
  echo "Namespace must be preview or local: $namespace" >&2
  exit 2
fi

if [[ "$mode" != "export" && "$mode" != "render" ]]; then
  echo "Mode must be export or render: $mode" >&2
  exit 2
fi

if ((${#selected_cases[@]} == 0)); then
  selected_cases=("${all_cases[@]}")
fi

for selected in "${selected_cases[@]}"; do
  found=false
  for known in "${all_cases[@]}"; do
    if [[ "$selected" == "$known" ]]; then
      found=true
      break
    fi
  done
  if [[ "$found" != true ]]; then
    echo "Unknown benchmark case: $selected" >&2
    exit 2
  fi
done

output_dir="${MOLFIG_BENCH_OUTPUT_DIR:-/tmp/molfig-benchmark}"
mkdir -p "$output_dir"

hyperfine_args=(
  --warmup "$warmup"
  --runs "$runs"
  --style basic
)

echo "Molfig @$namespace version: $version"
echo "Nixpkgs revision: ${MOLFIG_BENCH_NIXPKGS_REV:-unknown}"
echo "$(nix --version)"
echo "$typst_version"
echo "$hyperfine_version"
echo "System: $(uname -sm)"
echo "Benchmark mode: $mode"
echo "Warmup runs: $warmup; measured runs: $runs"

benchmark_index=0
for selected in "${selected_cases[@]}"; do
  benchmark_index=$((benchmark_index + 1))
  output="$output_dir/$selected-$mode.pdf"
  printf -v command \
    'typst compile --root %q --input %q --input %q --input %q --input %q %q %q' \
    "$repo_root" \
    "molfig-version=$version" \
    "molfig-namespace=$namespace" \
    "case=$selected" \
    "mode=$mode" \
    "$benchmark_source" \
    "$output"
  hyperfine "${hyperfine_args[@]}" --command-name "$selected" "$command" \
    | sed "s/^Benchmark 1: /Benchmark $benchmark_index: /"
done
