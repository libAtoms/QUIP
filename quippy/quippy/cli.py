import os
import sys
import subprocess as sp
import quippy
import argparse

def _run_command(name):
    """Run a bundled QUIP executable."""
    path = quippy.__path__[0]
    command = os.path.join(path, name)
    if not os.path.exists(command):
        print(f"Error: '{name}' executable not found at {command}", file=sys.stderr)
        print("This may indicate an incomplete installation.", file=sys.stderr)
        sys.exit(1)
    return sp.call([command] + sys.argv[1:])

def gap_fit():
    sys.exit(_run_command('gap_fit'))

def quip():
    sys.exit(_run_command('quip'))

def md():
    sys.exit(_run_command('md'))

def vasp_driver():
    sys.exit(_run_command('vasp_driver'))

def quip_config():
    parser = argparse.ArgumentParser(description='Configuration tool for QUIP')
    parser.add_argument('--libs', action='store_true', help="Arguments to link to libquip")
    args = parser.parse_args()

    if args.libs:
        libdir = quippy.__path__[0]
        print(f'-L{libdir} -lquip')
