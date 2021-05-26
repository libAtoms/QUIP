import sysconfig
from setuptools import setup, Extension
import re

major_version = '0.9'

# check for match with semantic version of the form a.b.c, plus optional suffix
semver_re = re.compile(r"^(\d+\.)?(\d+\.)?(\*|\d+)(-[A-Za-z0-9]+)*$")

with open('VERSION') as fin:
    version_string = fin.readline().strip()
    if version_string.startswith('v'):
        version_string = version_string[1:]
    if semver_re.match(version_string):
        version = version_string
    else:
        version_string = version_string.replace('-dirty', '.dirty')
        version = major_version + '+git' + version_string
print('version:', version)

platform = sysconfig.get_platform() + "-" + sysconfig.get_python_version()
ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='quippy-ase',
    version=version,
    maintainer='James Kermode',
    maintainer_email='james.kermode@gmail.com',
    description = 'ASE-compatible Python bindings for the QUIP and GAP codes',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',

        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'License :: Public Domain',
        'License :: Other/Proprietary License',

        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    url='https://github.com/libAtoms/QUIP',
    install_requires=['numpy>=1.16', 'f90wrap', 'ase'],
    python_requires=">=3.6",
    packages=['quippy'],
    package_data={'quippy': [f'../_quippy{ext_suffix}']},
    ext_modules=[Extension('_quippy', [])],
)
