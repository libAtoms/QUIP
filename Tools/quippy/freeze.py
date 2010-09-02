from cx_Freeze import setup, Executable

buildOptions = dict(excludes=['tcl','tk','Tkinter','ase','matplotlib'])

setup(name='convert',
      executables= [ Executable('scripts/convert.py') ],
      options = dict(build_exe=buildOptions))
