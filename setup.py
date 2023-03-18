"""SOMD: A SIESTA oriented Molecular Dynamics package

SOMD is an ab-initio molecular dynamics (AIMD) package designed for the SIESTA
(https://departments.icmab.es/leem/siesta/) code.
"""

import os
import sys
import versioneer
from setuptools import setup
from setuptools import Extension
from setuptools import find_packages
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Affero General Public License v3
Programming Language :: C++
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Chemistry
Operating System :: POSIX
Operating System :: Unix
"""

extensions = [
    Extension(
        'somd.core._lib',
        ['./somd/core/src/lib.pyx', './somd/core/src/math_utils.cxx'],
        include_dirs=['./somd/core/src'],
        extra_compile_args=['-std=c++11', '-Wall', '-fopenmp', '-O4',
                            '-march=native'],
        extra_link_args=['-fopenmp']),
    Extension(
        'somd.potentials._nepwrapper',
        ['./somd/potentials/src/nep.pyx'],
        include_dirs=['./somd/potentials/src'],
        extra_compile_args=['-std=c++11', '-Wall', '-fopenmp', '-O4',
                            '-march=native'],
        extra_link_args=['-fopenmp'])]

metadata = dict(
    name='somd',
    author='github.com/initqp',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=versioneer.get_version(),
    license='AGPLv3',
    url='https://www.github.com/initqp/somd',
    platforms=['Linux'],
    classifiers=CLASSIFIERS.splitlines(),
    packages=find_packages(),
    cmdclass=versioneer.get_cmdclass({'build_ext': build_ext}),
    install_requires=['Cython', 'mdtraj', 'toml'],
    zip_safe=False,
    ext_modules=cythonize(extensions),
    entry_points={'console_scripts': ['somd = somd.apps.cli:main']},
)

if __name__ == '__main__':
    setup(**metadata)
    if 'clean' in sys.argv[1:]:
        if os.path.isfile('somd/core/src/lib.cpp'):
            os.remove('somd/core/src/lib.cpp')
        if os.path.isfile('somd/potentials/src/nep.cpp'):
            os.remove('somd/potentials/src/nep.cpp')
