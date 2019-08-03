import glob
import subprocess
# from setuptools import setup
from distutils.core import setup

version = subprocess.getoutput(["../../bin/gitversion"]).strip()
print('version:', version)

setup(
    name='quippy',
    version=version,
    author='James Kermode <james.kermode@gmail.com>',
    packages=['quippy'],
    package_data={'quippy': glob.glob('_quippy.*.so')},
)
