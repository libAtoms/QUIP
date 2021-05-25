import sysconfig
from setuptools import setup, Extension

with open('VERSION') as fin:
    version='0.9+git'+fin.readline().strip().replace('-dirty', '.dirty')
print('version:', version)

platform = sysconfig.get_platform() + "-" + sysconfig.get_python_version()
ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")

setup(
    name='quippy-ase',
    version=version,
    author='James Kermode <james.kermode@gmail.com>',
    install_requires=['numpy', 'f90wrap', 'ase'],
    packages=['quippy'],
    package_data={'quippy': [f'../_quippy{ext_suffix}']},
    ext_modules=[Extension('_quippy', [])],
)
