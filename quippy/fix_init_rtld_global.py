#!/usr/bin/env python3
"""
Fix the import order in f90wrap-generated __init__.py to set RTLD_GLOBAL
before importing _quippy.

F90wrap generates imports at the top of __init__.py, but we need to set
RTLD_GLOBAL before importing _quippy.so to ensure symbol visibility for
runtime linking.
"""

import sys
import re
from pathlib import Path


def fix_init_rtld_global(init_file):
    """Move RTLD_GLOBAL setup to before the import quippy._quippy line."""

    content = init_file.read_text()

    # Find the RTLD_GLOBAL setup block
    rtld_pattern = r'(# IMPORTANT: Set RTLD_GLOBAL.*?sys\.setdlopenflags\(os\.RTLD_NOW \| os\.RTLD_GLOBAL\))'
    rtld_match = re.search(rtld_pattern, content, re.DOTALL)

    if not rtld_match:
        print(f"  ✗ RTLD_GLOBAL setup not found in {init_file}")
        return False

    rtld_block = rtld_match.group(1)

    # Find the restore block
    restore_pattern = r'(# Restore original dlopen flags.*?sys\.setdlopenflags\(_quippy_dlopen_flags\))'
    restore_match = re.search(restore_pattern, content, re.DOTALL)

    if not restore_match:
        print(f"  ✗ RTLD_GLOBAL restore not found in {init_file}")
        return False

    restore_block = restore_match.group(1)

    # Check if already in correct order (RTLD setup before import _quippy)
    import_quippy_pattern = r'^import quippy\._quippy$'
    import_match = re.search(import_quippy_pattern, content, re.MULTILINE)

    if not import_match:
        print(f"  ✗ import quippy._quippy not found in {init_file}")
        return False

    # Check if RTLD_GLOBAL is already before the import
    if rtld_match.start() < import_match.start():
        print(f"  - RTLD_GLOBAL already in correct position in {init_file}")
        return False

    # Remove both blocks from their current positions
    content_without_rtld = content.replace(rtld_block, '')
    content_without_rtld = content_without_rtld.replace(restore_block, '')

    # Find the import line again in the modified content
    import_match = re.search(import_quippy_pattern, content_without_rtld, re.MULTILINE)

    if not import_match:
        print(f"  ✗ import quippy._quippy disappeared after removing RTLD blocks")
        return False

    # Split at the import line
    before_import = content_without_rtld[:import_match.start()]
    import_line = import_match.group(0)
    after_import = content_without_rtld[import_match.end():]

    # Reconstruct with RTLD_GLOBAL before the import
    new_content = before_import + rtld_block + '\n\n' + import_line + '\n\n' + restore_block + after_import

    # Write back
    init_file.write_text(new_content)
    print(f"  ✓ Fixed RTLD_GLOBAL order in {init_file}")
    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: fix_init_rtld_global.py <__init__.py>")
        sys.exit(1)

    init_file = Path(sys.argv[1])
    if not init_file.exists():
        print(f"Error: {init_file} not found")
        sys.exit(1)

    try:
        if fix_init_rtld_global(init_file):
            sys.exit(0)
        else:
            sys.exit(0)  # No changes needed is not an error
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
