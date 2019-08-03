import subprocess
import glob
from setuptools import setup

version = subprocess.getoutput(["../../bin/gitversion"]).strip()
print('version:', version)

setup(
    name='quippy',
    version=version,
    author='James Kermode <james.kermode@gmail.com>',
    packages=['', 'quippy'],
    package_data={'': glob.glob('_quippy.*.so')},
)
