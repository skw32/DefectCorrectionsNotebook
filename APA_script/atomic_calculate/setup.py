from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules=[ Extension("atomic_V",
              ["atomic_V.pyx"],
              extra_compile_args = ["-ffast-math"])]

setup(
  name = "atomic_V",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules,
  include_dirs = [numpy.get_include()] 
)
