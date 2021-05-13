from setuptools import setup, Extension

with open('VERSION') as fin:
    version='0.9+git'+fin.readline().strip().replace('-dirty', '.dirty')
print('version:', version)

setup(
    name='quippy',
    version=version,
    author='James Kermode <james.kermode@gmail.com>',
    packages=['quippy'],
    package_data={'quippy': ['../_quippy.*.so']},
    ext_modules=[Extension('_quippy', [])],
)
