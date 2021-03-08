#!/usr/bin/env python
# Copyright (c) 2021, Lino Galiana
#
# Distributed under the 3-clause BSD license, see accompanying file LICENSE
# or https://github.com/scikit-hep/package for details.

# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.extension import Extension

from Cython.Build import cythonize
import numpy as np


try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = { }
ext_modules = [ ]

if use_cython:
    ext_modules = cythonize([ 
        Extension("btbpy.algo_smooth", ["btbpy/algo_smooth.pyx"], include_dirs=[np.get_include()])
    ])
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules = Extension("btbpy.algo_smooth", [ "btbpy/algo_smooth.cpp" ])



setup(
    cmdclass = cmdclass,
    ext_modules = ext_modules
)
