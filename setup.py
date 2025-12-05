"""
Setup script for py-fin package with f5c integration

This script builds the f5c eventalign module as a Python extension.
"""

import os
import sys
import shutil
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# f5c source directory (relative to setup.py)
# Use relative paths for all sources to comply with setuptools requirements


def apply_f5c_compile(ext):
    import numpy
    ext.extra_compile_args+=["-O3",          # Optimize compilation
        "-std=c99",     # C99 standard (matches your C code)
        "-Wall"         # Show warnings (debug)
        ]
    ext.extra_link_args += ["-lm"]  # Example link flag for Type1
    ext.include_dirs += [numpy.get_include()]
    ext.language = 'c'


def apply_dtw_compile(ext):
    nvcc_path = shutil.which('nvcc')
    if not nvcc_path:
        raise RuntimeError("nvcc not found! Please set CUDA_PATH or install CUDA Toolkit.")
    
    cuda_path = os.path.dirname(os.path.dirname(nvcc_path))
    include_path = os.path.join(cuda_path,'include')
    
    ext.extra_compile_args += [
        "-std=c++11",
        "-O3",
        "-Xcompiler", "-fPIC",
        "-Xcuda", "-gencode=arch=compute_80,code=sm_80",  #  80 = Ampere
    ]
    ext.extra_link_args += [
        "-lcudart",
        f"-L{cuda_path}/lib64",  # Linux/macOS
    ]
    ext.include_dirs.append(include_path)
    ext.language='c++'
    ext.compiler='nvcc'

# --------------------------
# Custom Build Extension
# --------------------------


class MultiExt(build_ext):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def build_extensions(self):
        # Apply type-specific logic to each extension
        for ext in self.extensions:
            if not hasattr(ext, "ext_type"):
                raise ValueError(f"Extension {ext.name} must have 'ext_type' (f5c/dtw)!")
            if ext.ext_type == "f5c":
                apply_f5c_compile(ext)
            elif ext.ext_type == "dtw":
                apply_dtw_compile(ext)
            else:
                raise ValueError(f"Unknown ext_type: {ext.ext_type} (must be f5c/dtw)")
        super().build_extensions()
    


# --------------------------
# Path Configuration
# --------------------------

F5C_DIR = os.path.join("fin", "_f5c")


f5c_extension = Extension(
    name="fin._f5c._event",
    sources=[
        os.path.join(F5C_DIR,'f5c_python.c'),
        os.path.join(F5C_DIR,'event_detection_simple.c')
    ],
    depends = [
      os.path.join(F5C_DIR,'event_detection_simple.h')  
    ],
    include_dirs=[F5C_DIR ],
    
)
f5c_extension.ext_type='f5c'
   
OPENDBA_DIR = os.path.join("fin",'_dtw')
cuda_dtw_extension = Extension(
        name="fin._dtw._cuda_dtw",  #
        sources=[
            os.path.join(OPENDBA_DIR,"dtw_api.cpp"),  
            
        ],
        depends = [
            os.path.join(OPENDBA_DIR,'cuda_utils.hpp'),
            os.path.join(OPENDBA_DIR,'dtw_api.h'),
            os.path.join(OPENDBA_DIR,'dtw.hpp'),
            os.path.join(OPENDBA_DIR,'limits.hpp'),
        ],
        
)
cuda_dtw_extension.ext_type='dtw'

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
            'fin._dtw'
        ],
        package_dir={"fin": "fin"},
        package_data={
            "fin": ["*.py", "*.c", "*.h", "*.yaml", "*.yml"],
        },
        include_package_data=True,
        ext_modules=[f5c_extension,cuda_dtw_extension], 
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
        entry_points={
            "console_scripts": [
                "fin=fin.cli:main",
            ],
        },
        zip_safe=False,
    )


if __name__ == "__main__":
    main()
