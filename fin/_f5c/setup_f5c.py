"""
Setup script for building f5c Python extensions.
"""

from setuptools import Extension, setup
import numpy as np

# Define the C extension
f5c_extension = Extension(
    'fin._f5c.f5c_python',
    sources=['f5c_python.c'],
    include_dirs=[
        np.get_include(),
        '.',  # Include current directory for headers
    ],
    libraries=[
        'm',  # Math library
    ],
    extra_compile_args=[
        '-O3',      # Optimize for speed
        '-std=c99',  # C99 standard
        '-Wall',     # Enable all warnings
        '-Wno-strict-prototypes',  # NumPy compatibility
    ],
    define_macros=[
        ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
    ]
)

setup(
    name='fin-f5c',
    version='0.1.0',
    ext_modules=[f5c_extension],
    include_dirs=[np.get_include()],
    python_requires='>=3.7',
)
