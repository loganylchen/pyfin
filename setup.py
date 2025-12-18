"""
Setup script for py-fin package with f5c integration

This script builds the f5c eventalign module as a Python extension.
CUDA extension (_dtw) is optional and requires nvcc.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# f5c source directory (relative to setup.py)
# Use relative paths for all sources to comply with setuptools requirements


def apply_f5c_compile(ext):
    import numpy

    ext.extra_compile_args += [
        "-O3",  # Optimize compilation
        "-std=c99",  # C99 standard (matches your C code)
        "-Wall",  # Show warnings (debug)
    ]
    ext.extra_link_args += ["-lm"]  # Example link flag for Type1
    ext.include_dirs += [numpy.get_include()]
    ext.language = "c"

def apply_align_cpp_compile(ext):
    import numpy
    ext.extra_compile_args += [
        "-x", "c++",
        "-O3",
        "-std=c++17",
        "-Wall",
    ]
    ext.extra_link_args += ["-lm"]
    ext.include_dirs += [numpy.get_include()]
    ext.language = "c++"

def find_cuda_home():
    """Find CUDA installation directory"""
    # Check environment variable first
    cuda_home = os.environ.get("CUDA_HOME") or os.environ.get("CUDA_PATH")
    if cuda_home and os.path.exists(cuda_home):
        return cuda_home

    # Try to find nvcc
    nvcc_path = shutil.which("nvcc")
    if nvcc_path:
        # nvcc is typically in $CUDA_HOME/bin/nvcc
        return os.path.dirname(os.path.dirname(os.path.realpath(nvcc_path)))

    # Common installation paths
    for path in ["/usr/local/cuda", "/usr/cuda", "/opt/cuda"]:
        if os.path.exists(path):
            return path

    return None


# --------------------------
# Custom Build Extension
# --------------------------


class CUDAExtension(Extension):
    """Custom Extension class that uses nvcc to compile CUDA code"""

    pass


class MultiExt(build_ext):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def build_extensions(self):
        # Filter out CUDA extensions if CUDA is not available
        cuda_available = find_cuda_home() is not None and shutil.which("nvcc") is not None

        extensions_to_build = []
        for ext in self.extensions:
            if not hasattr(ext, "ext_type"):
                raise ValueError(f"Extension {ext.name} must have 'ext_type' (f5c/dtw/f5c_cuda/align/align_cuda)!")

            if ext.ext_type == "dtw":
                if not cuda_available:
                    print(f"WARNING: Skipping CUDA extension {ext.name} - nvcc not found")
                    print("Install CUDA Toolkit to enable GPU acceleration features")
                    continue
                self._configure_cuda_extension(ext)
            elif ext.ext_type == "f5c_cuda":
                if not cuda_available:
                    print(
                        f"WARNING: Skipping CUDA eventalign extension {ext.name} - nvcc not found"
                    )
                    print("GPU-accelerated eventalign will not be available")
                    continue
                self._configure_f5c_cuda_extension(ext)
            elif ext.ext_type == "f5c":
                apply_f5c_compile(ext)
            elif ext.ext_type == "align":
                apply_align_cpp_compile(ext)
            elif ext.ext_type == "align_cuda":
                if not cuda_available:
                    print(
                        f"WARNING: Skipping CUDA eventalign extension {ext.name} - nvcc not found"
                    )
                    print("GPU-accelerated eventalign will not be available")
                    continue
                self._configure_align_cuda_extension(ext)
            else:
                raise ValueError(f"Unknown ext_type: {ext.ext_type} (must be f5c/align/dtw/f5c_cuda/align_cuda)")

            extensions_to_build.append(ext)

        self.extensions = extensions_to_build
        super().build_extensions()

    def _configure_f5c_cuda_extension(self, ext):
        """Configure f5c CUDA extension compilation (eventalign with GPU)"""
        cuda_home = find_cuda_home()
        if not cuda_home:
            raise RuntimeError("CUDA_HOME not found")

        import sysconfig
        import numpy

        python_include = sysconfig.get_path("include")
        numpy_include = numpy.get_include()

        # Add CUDA, Python, and NumPy include paths
        ext.include_dirs.append(os.path.join(cuda_home, "include"))
        ext.include_dirs.append(python_include)
        ext.include_dirs.append(numpy_include)
        ext.include_dirs.append("fin/_f5c")

        ext.library_dirs = [os.path.join(cuda_home, "lib64")]
        ext.libraries = ["cudart"]

        # Compiler flags for nvcc with CUDA_ENABLED defined
        ext.extra_compile_args = [
            "--compiler-options",
            "-fPIC",
            "-std=c++14",
            "-O3",
            "-DCUDA_ENABLED",  # Enable CUDA code paths in eventalign.c
            "--generate-code=arch=compute_80,code=sm_80",  # Ampere
        ]

        ext.extra_link_args = [
            f'-L{os.path.join(cuda_home, "lib64")}',
            "-lcudart",
        ]

    def _configure_cuda_extension(self, ext):
        """Configure CUDA extension compilation"""
        cuda_home = find_cuda_home()
        if not cuda_home:
            raise RuntimeError("CUDA_HOME not found")

        # Get Python and NumPy include directories
        import sysconfig
        import numpy

        python_include = sysconfig.get_path("include")
        numpy_include = numpy.get_include()

        # Add CUDA, Python, and NumPy include paths
        ext.include_dirs.append(os.path.join(cuda_home, "include"))
        ext.include_dirs.append(python_include)
        ext.include_dirs.append(numpy_include)
        ext.include_dirs.append("fin/_dtw")  # Add local include directory

        ext.library_dirs = [os.path.join(cuda_home, "lib64")]
        ext.libraries = ["cudart"]

        # Set compiler flags for nvcc
        ext.extra_compile_args = [
            "-x",
            "cu",  # Treat input as CUDA
            "--compiler-options",
            "-fPIC",
            "-std=c++11",
            "-O3",
            "--generate-code=arch=compute_80,code=sm_80",  # Ampere architecture
        ]

        ext.extra_link_args = [
            f'-L{os.path.join(cuda_home, "lib64")}',
            "-lcudart",
        ]

    def build_extension(self, ext):
        # Use nvcc for CUDA extensions
        if hasattr(ext, "ext_type") and ext.ext_type in ("dtw", "f5c_cuda"):
            self._compile_cuda_extension(ext)
        else:
            # Build non-CUDA extensions normally
            super().build_extension(ext)

    def _compile_cuda_extension(self, ext):
        """Compile CUDA extension using nvcc directly"""
        nvcc_path = shutil.which("nvcc")
        if not nvcc_path:
            raise RuntimeError("nvcc not found!")

        # Get output paths
        build_temp = Path(self.build_temp)
        build_lib = Path(self.build_lib)

        # Create directories
        build_temp.mkdir(parents=True, exist_ok=True)
        ext_path = build_lib / self.get_ext_filename(ext.name)
        ext_path.parent.mkdir(parents=True, exist_ok=True)

        # Compile each source file to object file
        objects = []
        for source in ext.sources:
            source_path = Path(source)
            # Include suffix to avoid name collisions (align.c vs align.cu)
            obj_name = source_path.stem + source_path.suffix.replace(".", "_") + ".o"
            obj_path = build_temp / obj_name

            # Determine file type and use appropriate compiler
            suffix = source_path.suffix.lower()

            if suffix == ".cu":
                # CUDA source file - compile with nvcc
                nvcc_cmd = [
                    nvcc_path,
                    "-c",
                    str(source_path),
                    "-o",
                    str(obj_path),
                ]
                # Add include directories
                for inc_dir in ext.include_dirs:
                    nvcc_cmd.extend(["-I", inc_dir])
                # Add compile flags
                nvcc_cmd.extend(ext.extra_compile_args)

                print(f"Compiling {source} with nvcc...")
                print(" ".join(nvcc_cmd))
                result = subprocess.run(nvcc_cmd, capture_output=True, text=True)

            elif suffix == ".c":
                # C source file - compile with gcc, add CUDA_ENABLED define
                import sysconfig

                cc = os.environ.get("CC", "gcc")

                c_cmd = [
                    cc,
                    "-c",
                    "-fPIC",
                    "-O3",
                    "-DCUDA_ENABLED",  # Enable CUDA code paths
                    str(source_path),
                    "-o",
                    str(obj_path),
                ]
                # Add include directories
                for inc_dir in ext.include_dirs:
                    c_cmd.extend(["-I", inc_dir])

                print(f"Compiling {source} with {cc}...")
                print(" ".join(c_cmd))
                result = subprocess.run(c_cmd, capture_output=True, text=True)

            else:
                # C++ source file - compile with nvcc treating as CUDA
                nvcc_cmd = [
                    nvcc_path,
                    "-x",
                    "cu",  # Treat as CUDA
                    "-c",
                    str(source_path),
                    "-o",
                    str(obj_path),
                ]
                # Add include directories
                for inc_dir in ext.include_dirs:
                    nvcc_cmd.extend(["-I", inc_dir])
                # Add compile flags
                nvcc_cmd.extend(ext.extra_compile_args)

                print(f"Compiling {source} with nvcc...")
                print(" ".join(nvcc_cmd))
                result = subprocess.run(nvcc_cmd, capture_output=True, text=True)

            if result.returncode != 0:
                print(result.stdout)
                print(result.stderr)
                raise RuntimeError(f"Compilation failed for {source}")

            objects.append(str(obj_path))

        # Link objects into shared library
        link_cmd = [
            nvcc_path,
            "-shared",
            "-o",
            str(ext_path),
        ] + objects

        # Add library directories and libraries
        if hasattr(ext, "library_dirs"):
            for lib_dir in ext.library_dirs:
                link_cmd.extend(["-L", lib_dir])

        if hasattr(ext, "libraries"):
            for lib in ext.libraries:
                link_cmd.append(f"-l{lib}")

        # Add link flags
        link_cmd.extend(ext.extra_link_args)

        print(f"Linking {ext.name}...")
        print(" ".join(link_cmd))
        result = subprocess.run(link_cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)
            raise RuntimeError(f"nvcc linking failed for {ext.name}")


# --------------------------
# Path Configuration
# --------------------------

F5C_DIR = os.path.join("fin", "_f5c")
ALIGN_DIR = os.path.join("fin", "_align")

align_extension = Extension(
    name="fin._align._align",
    sources=[
        os.path.join(ALIGN_DIR, "align_python.c"),
        os.path.join(ALIGN_DIR, "eventalign.c"),
        os.path.join(ALIGN_DIR, "align.c"),
        os.path.join(ALIGN_DIR, "hmm.c"),
        os.path.join(ALIGN_DIR, "model.c"),
        os.path.join(ALIGN_DIR, "nanopolish_read_db.c"),
        # Add other necessary source files here
    ],
    depends=[
        os.path.join(ALIGN_DIR, "f5c.h"),
        os.path.join(ALIGN_DIR, "f5cmisc.h"),
        os.path.join(ALIGN_DIR, "matrix.h"),
        os.path.join(ALIGN_DIR, "str.h"),
        os.path.join(ALIGN_DIR, "logsum.h"),
        os.path.join(ALIGN_DIR, "error.h"),
        os.path.join(ALIGN_DIR, "model.h"),
        os.path.join(ALIGN_DIR, "nanopolish_read_db.h"),
    ],
    include_dirs=[ALIGN_DIR],
    libraries=["hts", "z"],
)
align_extension.ext_type = "f5c_cpp"

f5c_extension = Extension(
    name="fin._f5c._event",
    sources=[
        os.path.join(F5C_DIR, "f5c_python.c"),
        os.path.join(F5C_DIR, "event_detection_simple.c"),
    ],
    depends=[os.path.join(F5C_DIR, "event_detection_simple.h")],
    include_dirs=[F5C_DIR],
)
f5c_extension.ext_type = "f5c"

# Eventalign extension (CPU only)
eventalign_extension = Extension(
    name="fin._f5c._eventalign",
    sources=[
        os.path.join(F5C_DIR, "eventalign.c"),
        os.path.join(F5C_DIR, "event_detection_simple.c"),
        os.path.join(F5C_DIR, "align.c"),  # CPU alignment with 3-state HMM
    ],
    depends=[
        os.path.join(F5C_DIR, "event_detection_simple.h"),
        os.path.join(F5C_DIR, "align_common.h"),
        os.path.join(F5C_DIR, "model.h"),
    ],
    include_dirs=[F5C_DIR],
)
eventalign_extension.ext_type = "f5c"

# CUDA Eventalign extension (GPU accelerated) - built only if CUDA is available
# Note: Needs align.c for profile_hmm_align (CPU implementation)
cuda_eventalign_extension = Extension(
    name="fin._f5c._eventalign_cuda",
    sources=[
        os.path.join(F5C_DIR, "eventalign.c"),
        os.path.join(F5C_DIR, "event_detection_simple.c"),
        os.path.join(F5C_DIR, "align.c"),  # CPU profile_hmm_align function
        os.path.join(F5C_DIR, "align.cu"),  # GPU align_with_flanking_gpu
    ],
    depends=[
        os.path.join(F5C_DIR, "event_detection_simple.h"),
        os.path.join(F5C_DIR, "align_common.h"),
        os.path.join(F5C_DIR, "model.h"),
    ],
    include_dirs=[F5C_DIR],
)
cuda_eventalign_extension.ext_type = "f5c_cuda"

# DTW/CUDA extension with Python bindings
OPENDBA_DIR = os.path.join("fin", "_dtw")
cuda_dtw_extension = Extension(
    name="fin._dtw._cuda_dtw",
    sources=[
        os.path.join(OPENDBA_DIR, "dtw_api.cpp"),
        os.path.join(OPENDBA_DIR, "multithreading.cpp"),
    ],
    depends=[
        os.path.join(OPENDBA_DIR, "cuda_utils.hpp"),
        os.path.join(OPENDBA_DIR, "dtw_api.h"),
        os.path.join(OPENDBA_DIR, "dtw.hpp"),
        os.path.join(OPENDBA_DIR, "limits.hpp"),
    ],
)
cuda_dtw_extension.ext_type = "dtw"


# Main setup configuration
def main():

    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()

    setup(
        name="py-fin",
        version="0.1.0",
        author="loganylchen",
        author_email="yuelong.chen.btr@gmail.com",
        description="A Python package for nanopore Direct RNA-seq data analysis and transcriptome assembly",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/loganylchen/pyfin",
        packages=["fin", "fin.io", "fin.utils", "fin.analysis", "fin._f5c", "fin._dtw", "fin._align"],
        package_dir={"fin": "fin"},
        package_data={
            "fin": ["*.py", "*.c", "*.cu", "*.h", "*.yaml", "*.yml"],
        },
        include_package_data=True,
        ext_modules=[
            f5c_extension,
            eventalign_extension,
            cuda_eventalign_extension,
            cuda_dtw_extension,
            align_extension,
        ],
        cmdclass={"build_ext": MultiExt},
        python_requires=">=3.8",
        install_requires=[
            "numpy>=1.21.0",
            "pandas>=1.3.0",
            "scipy>=1.7.0",
            "pysam>=0.21.0",
            "ont-fast5-api>=4.0.0",
            "pyslow5>=0.3.0",
            "pod5>=0.2.0",
            "h5py>=3.0.0",
            "click>=8.0.0",
            "tqdm>=4.62.0",
            "matplotlib>=3.4.0",
            "seaborn>=0.11.0",
            "pyyaml>=6.0",
        ],
        extras_require={
            "dev": [
                "pytest>=6.0.0",
                "pytest-cov",
                "black",
                "flake8",
            ],
            "gpu": [
                "cupy>=12.0.0",
                "numba>=0.56.0",
            ],
        },
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
        # CLI entry point - uncomment when fin.cli module is created
        # entry_points={
        #     "console_scripts": [
        #         "fin=fin.cli:main",
        #     ],
        # },
        zip_safe=False,
    )


if __name__ == "__main__":
    main()
