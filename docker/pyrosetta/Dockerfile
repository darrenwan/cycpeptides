# Academic PyRosetta (Non-Commercial) — independent of rosettacommons/rfdiffusion.
# Build: bash scripts/build_academic_images.sh pyrosetta
#
# Packages live in uv venv at /opt/venv (NOT system site-packages).
# Python pin rationale (official https://www.pyrosetta.org/downloads):
#   - "you will need to use Python-3.6 or later"
#   - conda packages: Python 3.6–3.14; wheels published per ABI (python39/310/311/…)
#   - Ubuntu wheels include PyRosetta4.Release.python310.ubuntu.wheel (first-class)
#   - Align with afcyc image on 3.10.15 for one academic Python line
#
# Preferred path: optional host-prefetched wheel under models/academic/
# (see models/academic/README.md). BuildKit bind-mounts it so the ~1.5 GiB
# file is not re-downloaded inside the build and is not stuffed into the context.
# Fallback: curl from official quarterly find-links (west/east).
# syntax=docker/dockerfile:1
ARG BASE_IMAGE=docker.m.daocloud.io/library/python:3.10-slim-bookworm
FROM ${BASE_IMAGE}

ARG PYTHON_VERSION=3.10.15
ARG UV_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
# Official recommended quarterly release mirror (West)
ARG PYROSETTA_FIND_LINKS=https://west.rosettacommons.org/pyrosetta/quarterly/release
# Stable "latest" alias for cp310 linux (~1.5 GiB). Override if ABI changes.
ARG PYROSETTA_WHEEL=pyrosetta-0-cp310-cp310-linux_x86_64.whl
# uv default HTTP timeout is 30s — far too low for a 1.5 GiB wheel.
ARG UV_HTTP_TIMEOUT=7200
# Official quarterly alias is named pyrosetta-0-*.whl while METADATA Version is
# e.g. 2026.29+releasequarterly.… — uv rejects that unless skipped.
ARG UV_SKIP_WHEEL_FILENAME_CHECK=1

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    UV_INDEX_URL=${UV_INDEX_URL} \
    UV_HTTP_TIMEOUT=${UV_HTTP_TIMEOUT} \
    UV_SKIP_WHEEL_FILENAME_CHECK=${UV_SKIP_WHEEL_FILENAME_CHECK} \
    VIRTUAL_ENV=/opt/venv \
    PATH="/opt/venv/bin:/root/.local/bin:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN curl -LsSf https://astral.sh/uv/install.sh | sh \
    && ln -sf /root/.local/bin/uv /usr/local/bin/uv

# Install exact CPython, then create venv (do not use bare python3)
RUN uv python install "${PYTHON_VERSION}" \
    && uv venv /opt/venv --python "${PYTHON_VERSION}" \
    && python -c "import sys; assert sys.version_info[:2]==(3,10), sys.version; print(sys.version)"

WORKDIR /opt/pdl1
RUN uv pip install --index-url "${UV_INDEX_URL}" \
      "numpy" "pyyaml"

# Official install target (pyrosetta.org/downloads): quarterly find-links.
# Dead URL (do NOT use): .../pyrosetta/release/release.html → 404
# Prefer host-prefetched wheel via additional build-context `wheels`
# (see build_academic_images.sh --build-context wheels=models/academic).
# Fallback: curl from official quarterly find-links (west/east).
RUN --mount=type=bind,from=wheels,source=.,target=/wheels,ro \
    set -eux; \
    echo "PYROSETTA_FIND_LINKS=${PYROSETTA_FIND_LINKS}"; \
    echo "PYROSETTA_WHEEL=${PYROSETTA_WHEEL}"; \
    echo "UV_HTTP_TIMEOUT=${UV_HTTP_TIMEOUT}"; \
    case "${PYROSETTA_FIND_LINKS}" in \
      *release/release.html*) \
        echo "ERROR: obsolete find-links URL (404). Use quarterly/release." >&2; exit 1 ;; \
    esac; \
    if [ -f "/wheels/${PYROSETTA_WHEEL}" ]; then \
      echo "Using host-prefetched wheel /wheels/${PYROSETTA_WHEEL}"; \
      ls -lh "/wheels/${PYROSETTA_WHEEL}"; \
      uv pip install "/wheels/${PYROSETTA_WHEEL}"; \
    else \
      BASE="${PYROSETTA_FIND_LINKS%/}"; \
      WHEEL_URL="${BASE}/${PYROSETTA_WHEEL}"; \
      echo "Host wheel missing; downloading ${WHEEL_URL}"; \
      curl -fL --retry 8 --retry-all-errors --retry-delay 15 \
        --connect-timeout 60 --max-time 0 \
        -C - \
        -o /tmp/pyrosetta.whl \
        "${WHEEL_URL}" \
      || (echo "ERROR: curl download failed: ${WHEEL_URL}" >&2; \
          echo "  Run scripts/download_pyrosetta_wheel.sh on the host first." >&2; \
          exit 1); \
      ls -lh /tmp/pyrosetta.whl; \
      uv pip install /tmp/pyrosetta.whl; \
      rm -f /tmp/pyrosetta.whl; \
    fi

RUN python -c "import sys, pyrosetta; \
print('python', sys.version); \
print('venv', sys.prefix); \
print('pyrosetta', getattr(pyrosetta,'__version__','?'))"

WORKDIR /work
ENTRYPOINT ["python"]
