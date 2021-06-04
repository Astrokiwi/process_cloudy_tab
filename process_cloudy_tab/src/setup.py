from distutils.core import setup, Extension
import numpy

extension_mod = Extension("_tab_interp", ["_tab_interp_module.cc", "dust_temp_interp.cpp"])

setup(name = "tab_interp", ext_modules=[extension_mod],include_dirs = [numpy.get_include(),'.'])
