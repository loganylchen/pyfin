"""
Setup script for py-fin package with f5c integration

This script builds the f5c eventalign module as a Python extension.
"""

import os
import sys
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# f5c source directory
F5C_DIR = Path(__file__).parent / "third_party" / "f5c"
SLOW5LIB_DIR = F5C_DIR / "slow5lib"

# Compiler arguments
extra_compile_args = [
    "-std=c++11",
    "-O3",
    "-g",
    "-Wall",
    "-fPIC",
    "-D HAVE_CUDA=0",  # Build without CUDA for now
    "-D DISABLE_HDF5=1",  # Disable HDF5/FAST5 support, use slow5 instead
]

# Linker arguments
extra_link_args = [
    "-lpthread",
    "-lz",
    "-lrt",
    "-ldl",
]

# f5c source files
f5c_sources = [
    str(F5C_DIR / "src" / "f5c.c"),
    str(F5C_DIR / "src" / "events.c"),
    str(F5C_DIR / "src" / "nanopolish_read_db.c"),
    str(F5C_DIR / "src" / "index.c"),
    str(F5C_DIR / "src" / "nanopolish_fast5_io.c"),
    str(F5C_DIR / "src" / "model.c"),
    str(F5C_DIR / "src" / "methmodel.c"),
    str(F5C_DIR / "src" / "align.c"),
    str(F5C_DIR / "src" / "hmm.c"),
    str(F5C_DIR / "src" / "meth.c"),
    str(F5C_DIR / "src" / "freq.c"),
    str(F5C_DIR / "src" / "eventalign.c"),
    str(F5C_DIR / "src" / "freq_merge.c"),
    str(F5C_DIR / "src" / "resquiggle.c"),
    str(F5C_DIR / "src" / "profiles.c"),
    # Add slow5lib source files
    str(SLOW5LIB_DIR / "src" / "slow5.c"),
    str(SLOW5LIB_DIR / "src" / "slow5_index.c"),
    str(SLOW5LIB_DIR / "src" / "slow5_press.c"),
    str(SLOW5LIB_DIR / "src" / "slow5_misc.c"),
    # Python wrapper
    "fin/_f5c/f5c_python.c",
]

# Check for required libraries
def check_dependencies():
    """Check if required system libraries are available"""
    missing = []

    # Check for zlib (system library, not Python package)
    # Note: Python's zlib module doesn't guarantee the dev headers are installed
    # We mainly need the zlib development headers for compilation

    # Check for system zlib headers
    import subprocess
    try:
        result = subprocess.run(['pkg-config', '--exists', 'zlib'], capture_output=True)
        if result.returncode != 0:
            missing.append("zlib development headers (install zlib1g-dev)")
    except FileNotFoundError:
        # pkg-config not available, just warn
        pass

    # Don't check for pysam here - it will be installed automatically by pip
    # and checking now would give false warnings during installation
    # The real pysam import check happens at runtime when someone tries to use the modules

    # However, we can check if htslib headers are available for compilation
    try:
        result = subprocess.run(['pkg-config', '--exists', 'htslib'], capture_output=True)
        if result.returncode != 0:
            missing.append("htslib development headers (install libhts-dev or libhts-devel)")
    except FileNotFoundError:
        # pkg-config not available, assume headers are present (or will be handled by pysam)
        pass

    if missing:
        print(f"Warning: Some build dependencies may be missing: {', '.join(missing)}")
        print("These are only needed for compiling the f5c extension.")
        print("If you're just installing the Python package, you can ignore this.")
        print("\nTo install build dependencies on Ubuntu/Debian:")
        print("  sudo apt-get install zlib1g-dev libhts-dev")
        print("\nOn CentOS/RHEL:")
        print("  sudo yum install zlib-devel htslib-devel")
        print("\nThe following Python packages will be installed automatically:")
        print("  - pysam (includes htslib Python bindings)")
        print("  - numpy, scipy, pandas, etc.")


class BuildF5CExt(build_ext):
    """Custom build extension for f5c"""

    def run(self):
        # Build slow5lib first
        self.build_slow5lib()
        super().run()

    def build_extensions(self):
        # Add numpy include dirs here (after numpy is installed)
        try:
            import numpy as np
            numpy_include = np.get_include()
            for ext in self.extensions:
                ext.include_dirs.append(numpy_include)
        except ImportError:
            print("Warning: numpy not available, extension may not compile correctly")
            print("Make sure numpy is installed before building")

        super().build_extensions()

    def build_slow5lib(self):
        """Build slow5lib static library"""
        import subprocess

        slow5_build_dir = SLOW5LIB_DIR / "build"
        slow5_build_dir.mkdir(exist_ok=True)

        print("Building slow5lib...")

        # Configure and build slow5lib
        env = os.environ.copy()
        env["CFLAGS"] = "-fPIC -O3"

        try:
            subprocess.run(
                ["make", "-C", str(SLOW5LIB_DIR), "lib"],
                env=env,
                check=True
            )
            print("slow5lib built successfully")
        except subprocess.CalledProcessError as e:
            print(f"Failed to build slow5lib: {e}")
            raise


# Initialize include_dirs base (numpy will be added during build)
include_dirs_base = [
    str(F5C_DIR / "include"),
    str(F5C_DIR / "src"),
    str(SLOW5LIB_DIR / "include"),
]

# f5c extension module
f5c_module = Extension(
    "fin._f5c",
    sources=f5c_sources,
    include_dirs=include_dirs_base + [
        str(F5C_DIR / "slow5lib" / "include"),
    ],
    libraries=["z", "pthread", "m"],
    library_dirs=[],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args + [
        str(SLOW5LIB_DIR / "lib" / "libslow5.a")
    ],
    language="c++",
)

# Main setup configuration
def main():
    check_dependencies()

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
        packages=[
            "fin",
            "fin.io",
            "fin.utils",
            "fin.core",
        ],
        package_dir={"fin": "fin"},
        package_data={
            "fin": ["*.py", "*.c", "*.h", "*.yaml", "*.yml"],
        },
        ext_modules=[f5c_module],
        cmdclass={"build_ext": BuildF5CExt},
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
        entry_points={
            "console_scripts": [
                "fin=fin.cli:main",
            ],
        },
        zip_safe=False,
    )


if __name__ == "__main__":
    main()
