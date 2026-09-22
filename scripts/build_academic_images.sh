#!/usr/bin/env bash
# Build academic Docker images used by the PD-L1 campaign.
#
# Usage:
#   bash scripts/build_academic_images.sh afcyc
#   bash scripts/build_academic_images.sh pyrosetta
#   bash scripts/build_academic_images.sh all
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TARGET="${1:-}"
if [[ -z "${TARGET}" ]]; then
  echo "Usage: $0 {afcyc|pyrosetta|all}" >&2
  exit 2
fi

build_afcyc() {
  local tag="${AFCYC_IMAGE:-pdl1-afcyc:latest}"
  echo "=== build ${tag} ==="
  bash "${ROOT}/scripts/fetch_colabdesign.sh"
  if [[ ! -f "${ROOT}/third_party/ColabDesign/setup.py" \
     && ! -f "${ROOT}/third_party/ColabDesign/pyproject.toml" ]]; then
    echo "Missing third_party/ColabDesign after fetch." >&2
    exit 1
  fi
  docker build \
    -f "${ROOT}/docker/afcyc/Dockerfile" \
    -t "${tag}" \
    "${ROOT}"
  echo "OK → ${tag}"
  if [[ -f "${ROOT}/third_party/ColabDesign/.colabdesign_pin" ]]; then
    echo "ColabDesign pin:"
    cat "${ROOT}/third_party/ColabDesign/.colabdesign_pin"
  fi
}

build_pyrosetta() {
  local tag="${PYROSETTA_IMAGE:-pdl1-pyrosetta:latest}"
  # Optional host wheel cache (not required). Override: PYROSETTA_WHEELS_DIR=/path/to/dir
  local wheels="${PYROSETTA_WHEELS_DIR:-${ROOT}/models/academic}"
  local wheel_name="${PYROSETTA_WHEEL:-pyrosetta-0-cp310-cp310-linux_x86_64.whl}"
  local ctx="${wheels}"
  local tmp=""

  echo "=== build ${tag} ==="
  if [[ -f "${wheels}/${wheel_name}" ]]; then
    echo "Using host wheel: ${wheels}/${wheel_name}"
  else
    echo "No host wheel at ${wheels}/${wheel_name}"
    echo "  → Docker will download from official quarterly find-links (slow)."
    echo "  → To prefetch later: place the .whl under models/academic/ (see README there)."
    if [[ ! -d "${wheels}" ]]; then
      tmp="$(mktemp -d "${TMPDIR:-/tmp}/pdl1-pyrosetta-wheels.XXXXXX")"
      ctx="${tmp}"
      trap 'rm -rf "${tmp}"' RETURN
    fi
  fi

  docker build \
    -f "${ROOT}/docker/pyrosetta/Dockerfile" \
    --build-context "wheels=${ctx}" \
    -t "${tag}" \
    "${ROOT}"
  echo "OK → ${tag}"
}

case "${TARGET}" in
  afcyc) build_afcyc ;;
  pyrosetta) build_pyrosetta ;;
  all)
    build_afcyc
    build_pyrosetta
    ;;
  *)
    echo "Unknown target: ${TARGET} (want afcyc|pyrosetta|all)" >&2
    exit 2
    ;;
esac
