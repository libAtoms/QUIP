#!/usr/bin/env python3
"""Install QUIP executables into the quippy package directory."""
import os
import sys
import shutil
import stat

def main():
    # Arguments: source_dir, dest_dir
    if len(sys.argv) < 3:
        print("Usage: install_executables.py <source_dir> <dest_dir>")
        sys.exit(1)

    source_dir = sys.argv[1]
    dest_dir = sys.argv[2]

    # Handle DESTDIR for meson install (used for staged installs / wheel building)
    destdir = os.environ.get('DESTDIR', '')
    if destdir:
        # dest_dir is absolute, so we need to handle the join carefully
        if dest_dir.startswith('/'):
            dest_dir = destdir + dest_dir
        else:
            dest_dir = os.path.join(destdir, dest_dir)

    executables = ['quip', 'gap_fit', 'md']

    os.makedirs(dest_dir, exist_ok=True)

    for exe in executables:
        src = os.path.join(source_dir, exe)
        dst = os.path.join(dest_dir, exe)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            # Ensure executable permissions
            os.chmod(dst, os.stat(dst).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            print(f"Installed {exe} to {dst}")
        else:
            print(f"Warning: {exe} not found at {src}", file=sys.stderr)

if __name__ == '__main__':
    main()
