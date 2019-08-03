import subprocess
from setuptools import setup, Extension

version = subprocess.getoutput(["../../bin/gitversion"]).strip()
print('version:', version)

setup(
    name='quippy',
    version=version,
    author='James Kermode <james.kermode@gmail.com>',
    packages=['quippy'],
    package_data={'quippy': ['../_quippy.*.so']},
    ext_modules=[Extension('_quippy', [])],
)
