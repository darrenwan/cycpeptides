# AfCycDesign runtime (ColabDesign + JAX CUDA). Independent of RFdiffusion image.
# Build: bash scripts/build_academic_images.sh afcyc
# Requires host: third_party/ColabDesign (bash scripts/fetch_colabdesign.sh).
# AF2 params NOT baked — mount models/afcyc → /models/afcyc at runtime.
#
# Pins (models/SOURCES.md): CPython 3.10.15, jax[cuda12]==0.4.30, dm-haiku==0.0.12,
# flax==0.8.5, numpy==1.26.4, biopython==1.84, ColabDesign v1.1.3.
# haiku≥0.0.17 needs jax.interpreters.xla.xe (gone in jax≥0.5) → stay on 0.4.30.
# jax[cuda12] pulls cudnn9 pip wheels (base image is cudnn8) — intentional.
# syntax=docker/dockerfile:1
ARG BASE_IMAGE=docker.m.daocloud.io/nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04
FROM ${BASE_IMAGE}

ARG PYTHON_VERSION=3.10.15
ARG UV_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
ARG JAX_VERSION=0.4.30
ARG DM_HAIKU_VERSION=0.0.12
ARG FLAX_VERSION=0.8.5
ARG NUMPY_VERSION=1.26.4
ARG BIOPYTHON_VERSION=1.84
ARG UV_HTTP_TIMEOUT=7200

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    UV_INDEX_URL=${UV_INDEX_URL} \
    UV_HTTP_TIMEOUT=${UV_HTTP_TIMEOUT} \
    VIRTUAL_ENV=/opt/venv \
    PATH="/opt/venv/bin:/root/.local/bin:${PATH}" \
    AFCYC_PARAMS_DIR=/models/afcyc

RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN curl -LsSf https://astral.sh/uv/install.sh | sh \
    && ln -sf /root/.local/bin/uv /usr/local/bin/uv

RUN uv python install "${PYTHON_VERSION}" \
    && uv venv /opt/venv --python "${PYTHON_VERSION}"

# Host-prefetched only (no git-in-Docker fallback).
COPY third_party/ColabDesign /opt/ColabDesign

# ColabDesign --no-deps so PyPI does not pull jax≥0.5; then one jax[cuda12] pin.
RUN set -eux; \
    test -f /opt/ColabDesign/setup.py || test -f /opt/ColabDesign/pyproject.toml; \
    uv pip install --no-deps /opt/ColabDesign; \
    uv pip install --index-url "${UV_INDEX_URL}" \
      "numpy==${NUMPY_VERSION}" \
      "biopython==${BIOPYTHON_VERSION}" \
      "pyyaml" \
      "py3Dmol" "absl-py" "chex" "dm-tree" "immutabledict" \
      "ml-collections" "pandas" "scipy" "optax" "joblib" "matplotlib" \
      "jmp" "tabulate" "flax==${FLAX_VERSION}"; \
    uv pip install --index-url "${UV_INDEX_URL}" --no-deps \
      "dm-haiku==${DM_HAIKU_VERSION}"; \
    uv pip install \
      "jax[cuda12]==${JAX_VERSION}" \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html \
    || (echo "ERROR: jax[cuda12]==${JAX_VERSION} failed (needs GCS)." >&2; exit 1); \
    python -c "import sys, jax, haiku, colabdesign, Bio; \
assert sys.version_info[:2]==(3,10), sys.version; \
assert jax.__version__=='${JAX_VERSION}', jax.__version__; \
print('python', sys.version.split()[0], 'jax', jax.__version__, \
      'haiku', getattr(haiku,'__version__','?'), 'Bio', Bio.__version__)"

WORKDIR /work
ENTRYPOINT ["python"]
