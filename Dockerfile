# Multi-stage Dockerfile for FIN (Nanopore Isoform Detection)
# Supports both GPU (CUDA) and CPU-only builds

ARG CUDA_VERSION=11.8.0
ARG UBUNTU_VERSION=22.04

# =============================================================================
# Base stage with common dependencies
# =============================================================================
FROM nvidia/cuda:${CUDA_VERSION}-devel-ubuntu${UBUNTU_VERSION} AS base

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libhdf5-dev \
    hdf5-tools \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN groupadd -r finuser && useradd -r -g finuser finuser

# Set working directory
WORKDIR /app

# =============================================================================
# GPU stage - with CUDA support
# =============================================================================
FROM base AS gpu

# Install CUDA-specific dependencies
RUN apt-get update && apt-get install -y \
    cuda-nvcc-${CUDA_VERSION} \
    libcublas-dev \
    libcusparse-dev \
    && rm -rf /var/lib/apt/lists/*

# Set CUDA environment variables
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=${CUDA_HOME}/bin:${PATH}
ENV LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}

# Install PyTorch with CUDA support for OpenDBA compilation
RUN pip3 install --no-cache-dir torch torchvision --index-url https://download.pytorch.org/whl/cu118

# Copy source code
COPY --chown=finuser:finuser . .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Install the package with CUDA extensions
RUN pip3 install --no-cache-dir -e .

# Test CUDA availability
RUN python3 -c "from fin._opendba import opendba_cuda; print('CUDA extension loaded successfully')"

# =============================================================================
# CPU stage - without CUDA
# =============================================================================
FROM base AS cpu

# Install CPU-only PyTorch (won't be used but prevents import errors)
RUN pip3 install --no-cache-dir torch torchvision --index-url https://download.pytorch.org/whl/cpu

# Copy source code
COPY --chown=finuser:finuser . .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Install without CUDA (CPU-only)
RUN OPENDBA_USE_CUDA=0 pip3 install --no-cache-dir -e .

# Test CPU fallback
RUN python3 -c "from fin.core.dtw_gpu import OPENDBA_AVAILABLE; print(f'CUDA available: {OPENDBA_AVAILABLE}')"

# =============================================================================
# Runtime stage for GPU
# =============================================================================
FROM nvidia/cuda:${CUDA_VERSION}-runtime-ubuntu${UBUNTU_VERSION} AS runtime-gpu

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    libhdf5-103 \
    && rm -rf /var/lib/apt/lists/*

# Copy installed package from GPU build stage
COPY --from=gpu /usr/local/lib/python3.*/dist-packages/fin* /usr/local/lib/python3/dist-packages/
COPY --from=gpu /usr/local/lib/python3.*/dist-packages/fin-*.egg-info /usr/local/lib/python3/dist-packages/
COPY --from=gpu /usr/local/lib/python3.*/dist-packages/fin*so /usr/local/lib/python3/dist-packages/
COPY --from=gpu /app/fin /app/fin

WORKDIR /app

# Set up non-root user
RUN groupadd -r finuser && useradd -r -g finuser finuser
RUN chown -R finuser:finuser /app
USER finuser

# Default command
CMD ["python3", "-c", "import fin; print('FIN GPU version ready')"]

# =============================================================================
# Runtime stage for CPU
# =============================================================================
FROM ubuntu:${UBUNTU_VERSION} AS runtime-cpu

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    libhdf5-103 \
    && rm -rf /var/lib/apt/lists/*

# Copy installed package from CPU build stage
COPY --from=cpu /usr/local/lib/python3.*/dist-packages/fin* /usr/local/lib/python3/dist-packages/
COPY --from=cpu /usr/local/lib/python3.*/dist-packages/fin-*.egg-info /usr/local/lib/python3/dist-packages/
COPY --from=cpu /usr/local/lib/python3.*/dist-packages/fin*so /usr/local/lib/python3/dist-packages/
COPY --from=cpu /app/fin /app/fin

WORKDIR /app

# Set up non-root user
RUN groupadd -r finuser && useradd -r -g finuser finuser
RUN chown -R finuser:finuser /app
USER finuser

# Default command
CMD ["python3", "-c", "import fin; print('FIN CPU version ready')"]

# Use CPU runtime as default
FROM runtime-cpu
