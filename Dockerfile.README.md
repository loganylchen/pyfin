# Docker Setup for FIN

FIN provides separate Docker images for CPU-only and GPU-enabled usage.

## Images

### CPU-Only Image (Recommended for most users)
- **Image**: `loganylchen/pyfin:latest`
- **Dockerfile**: `Dockerfile.cpu`
- **Size**: ~500MB (much smaller than GPU version)
- **Platforms**: linux/amd64, linux/arm64
- **Use case**: General use, laptops, clusters without GPU

```bash
docker pull loganylchen/pyfin:latest
docker run -it --rm loganylchen/pyfin:latest
```

### GPU-Enabled Image
- **Image**: `loganylchen/pyfin:latest-gpu`
- **Dockerfile**: `Dockerfile.gpu`
- **Size**: ~2.5GB (includes CUDA and PyTorch)
- **Platforms**: linux/amd64
- **Use case**: High-performance computing with CUDA-enabled GPU

```bash
docker pull loganylchen/pyfin:latest-gpu
docker run -it --rm --gpus all loganylchen/pyfin:latest-gpu
```

## Building Images Locally

### CPU Image
```bash
docker build -f Dockerfile.cpu -t pyfin:cpu-local .
```

### GPU Image
```bash
docker build -f Dockerfile.gpu -t pyfin:gpu-local .
```

## Key Differences

### CPU Dockerfile (`Dockerfile.cpu`)
- Based on Ubuntu
- No CUDA dependencies
- No PyTorch installation
- Pure CPU implementation for DTW calculations
- Smaller image size (faster pull, lower storage)

### GPU Dockerfile (`Dockerfile.gpu`)
- Based on NVIDIA CUDA base image
- Includes CUDA 11.8 development tools
- Installs PyTorch with CUDA support (used for compilation)
- GPU-accelerated DTW via OpenDBA CUDA kernels
- Faster computation for large datasets

## Runtime Requirements

### CPU Image
- Docker or compatible container runtime

### GPU Image
- NVIDIA GPU with compute capability 3.5+
- NVIDIA Container Toolkit installed
- Docker with GPU support

## Testing Images

After building, test the images:

```bash
# Test GPU image
docker run --rm --gpus all loganylchen/pyfin:latest-gpu \
  python3 -c "from fin._opendba import opendba_cuda; print('CUDA extension loaded successfully')"

# Test CPU image
docker run --rm loganylchen/pyfin:latest \
  python3 -c "import fin; from fin.core.dtw_gpu import OPENDBA_AVAILABLE; print(f'CUDA available: {OPENDBA_AVAILABLE}')"
```

## CI/CD

GitHub Actions automatically builds and pushes both images on pushes to main branch:
- `docker-gpu` job builds and pushes GPU image
- `docker-cpu` job builds and pushes CPU image (with multi-arch support)

Both images are tagged appropriately:
- `latest` → CPU image
- `latest-gpu` → GPU image
- Version-specific tags from git commits and branches
