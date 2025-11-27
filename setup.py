from setuptools import setup, find_packages, Extension
import sys

# Try to import CUDA-enabled build tools
try:
    from torch.utils.cpp_extension import BuildExtension, CUDAExtension
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

# Try to import numpy - may not be available during initial setup
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    np = None
    NUMPY_AVAILABLE = False
    # Only error if we're actually building extensions
    if 'build_ext' in sys.argv or 'install' in sys.argv or 'bdist_wheel' in sys.argv:
        print("ERROR: numpy is required to build extensions. Install it first with: pip install numpy")
        sys.exit(1)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

# Delay numpy import for include directories
def get_numpy_include():
    global np
    global NUMPY_AVAILABLE
    if not NUMPY_AVAILABLE:
        import numpy as _np
        np = _np
        NUMPY_AVAILABLE = True
    return np.get_include()

# f5c C extension - using simplified standalone event detection
f5c_extension = Extension(
    'fin._f5c.f5c_python',
    sources=[
        'fin/_f5c/f5c_python.c',
        'fin/_f5c/event_detection_simple.c',
    ],
    include_dirs=[
        get_numpy_include(),
        'fin/_f5c',
    ],
    libraries=['m'],
    extra_compile_args=[
        '-O3',
        '-std=c99',
        '-Wall',
        '-Wno-strict-prototypes',
    ],
    define_macros=[
        ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
    ]
)

# OpenDBA CUDA extension
if TORCH_AVAILABLE:
    # Use PyTorch's CUDA build system if available
    opendba_extension = CUDAExtension(
        'fin._opendba.opendba_cuda',
        sources=[
            'fin/_opendba/opendba_python.cpp',
            'fin/_opendba/openDBA.cu',
        ],
        include_dirs=[
            np.get_include(),
            'fin/_opendba',
        ],
        extra_compile_args={
            'cxx': ['-O3'],
            'nvcc': ['-O3'],
        },
        define_macros=[
            ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
        ]
    )
else:
    # CPU-only version
    opendba_extension = Extension(
        'fin._opendba.opendba_cuda',
        sources=['fin/_opendba/opendba_python.cpp'],
        include_dirs=[
            np.get_include(),
            'fin/_opendba',
        ],
        extra_compile_args=['-O3', '-std=c++11'],
        define_macros=[
            ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
        ],
        language='c++'
    )

setup(
    name="fin",
    version="0.1.0",
    author="fin authors",
    author_email="",
    description="A Python tool for detecting RNA modifications using nanopore Direct RNA-seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fin/fin",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    setup_requires=['numpy'],
    ext_modules=[f5c_extension, opendba_extension],
    cmdclass={'build_ext': BuildExtension} if TORCH_AVAILABLE else {},
    entry_points={
        "console_scripts": [
            "FIN.py=fin.cli:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
