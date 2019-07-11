# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.extension import Extension

from Cython.Build import cythonize
import numpy as np

with open("README.md", "r") as fh:
    long_description = fh.read()

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
    name="btbpy",
    version="0.0.5",
    author="Julien Jamme, Arlindo Dos Santos, François Sémécurbe",
    maintainer="Julien Jamme",
    maintainer_email="julien.jamme@protonmail.com",
    description="Smooth geographical data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/mrteste/btbpy",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        'Programming Language :: Cython',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Operating System :: OS Independent",
    ],
    packages = find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    license = 'GPL3',
    keywords="smooth greographic kernel smoothing",
    install_requires = ['numpy'],
    zip_safe=False,
    setup_requires = ['Cython'],
    cmdclass = cmdclass,
    ext_modules = ext_modules,
    package_data={
        'btbpy': ['data/*.csv']
    }
)
