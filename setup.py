"""
Setup script for py-fin package with f5c integration

This script builds the f5c eventalign module as a Python extension.
"""

import os
import sys
import numpy as np
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# f5c source directory (relative to setup.py)
# Use relative paths for all sources to comply with setuptools requirements
NUMPY_INCLUDE = np.get_include()
if not os.path.exists(os.path.join(NUMPY_INCLUDE, "numpy", "arrayobject.h")):
    raise RuntimeError(
        f"NumPy headers not found at {NUMPY_INCLUDE}! "
        "Ensure numpy is installed in the current Python environment."
    )

# --------------------------
# Path Configuration
# --------------------------
PACKAGE_ROOT = os.path.dirname(os.path.abspath(__file__))
F5C_DIR = os.path.join(PACKAGE_ROOT, "fin", "_f5c")


f5c_extension = Extension(
    name="fin._f5c",
    sources=[
        os.path.join(F5C_DIR,'f5c_python.c'),
        os.path.join(F5C_DIR,'event_detection_simple.c')
    ],
    inlcude_dirs=[F5C_DIR,np.get_include() ],
    extra_compile_args=["-O3",          # Optimize compilation
        "-std=c99",     # C99 standard (matches your C code)
        "-Wall"         # Show warnings (debug)
        ],
    extra_link_args=["-lm"],
    language="c"
)
   


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
        packages=[
            "fin",
            "fin.io",
            "fin.utils",
            "fin.core",
            "fin._f5c",
        ],
        package_dir={"fin": "fin"},
        package_data={
            "fin": ["*.py", "*.c", "*.h", "*.yaml", "*.yml"],
        },
        ext_modules=[f5c_extension], 
        # cmdclass={"build_ext": BuildF5CExt},
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
