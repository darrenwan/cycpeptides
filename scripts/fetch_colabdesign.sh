#!/usr/bin/env bash
# Prefetch ColabDesign into third_party/ for docker/afcyc builds (CN github-friendly).
#
# Usage:
#   bash scripts/fetch_colabdesign.sh
#   COLABDESIGN_REF=v1.1.2 bash scripts/fetch_colabdesign.sh
#   FORCE=1 bash scripts/fetch_colabdesign.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DST="${ROOT}/third_party/ColabDesign"
REF="${COLABDESIGN_REF:-v1.1.3}"
REPO="${COLABDESIGN_REPO:-https://github.com/sokrypton/ColabDesign.git}"
FORCE="${FORCE:-0}"

if [[ -f "${DST}/setup.py" || -f "${DST}/pyproject.toml" ]]; then
  if [[ "${FORCE}" != "1" ]]; then
    echo "skip ${DST} (exists; FORCE=1 to re-clone)"
    if [[ -f "${DST}/.colabdesign_pin" ]]; then
      cat "${DST}/.colabdesign_pin"
    fi
    exit 0
  fi
  rm -rf "${DST}"
fi

mkdir -p "${ROOT}/third_party"
echo "→ clone ${REPO} @ ${REF} → third_party/ColabDesign"
git clone --depth 1 --branch "${REF}" "${REPO}" "${DST}"

commit="$(git -C "${DST}" rev-parse HEAD)"
{
  echo "repo=${REPO}"
  echo "ref=${REF}"
  echo "commit=${commit}"
  echo "fetched=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} | tee "${DST}/.colabdesign_pin"

echo "OK ${DST} @ ${REF} (${commit})"
