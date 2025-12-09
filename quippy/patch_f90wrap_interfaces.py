#!/usr/bin/env python3
"""
Post-process f90wrap-generated Python files to fix bugs in f90wrap output.

Fixes:
1. Overloaded interface shadowing bug - F90wrap generates overloaded interfaces
   that shadow the methods they're trying to call.
2. Finalizer initialise/finalise bug - F90wrap's fallback finalizer incorrectly
   uses _initialise instead of _finalise due to a pattern matching bug.
"""

import re
import sys
from pathlib import Path


def fix_finalizer_initialise_bug(content):
    """Fix finalizer bug where _initialise is used instead of _finalise.

    F90wrap's pywrapgen.py has a bug in write_constructor where it tries to
    replace "__initialise" with "__finalise" for the fallback finalizer, but
    the actual pattern in subroutine names is "_initialise" (single underscore).

    This results in finalizers calling the initialise function instead of
    finalise, causing errors on cleanup when garbage collection runs.

    Bug pattern:
        self._finalizer = weakref.finalize(self, quippy._quippy.f90wrap_..._initialise, self._handle)

    Fixed pattern:
        self._finalizer = weakref.finalize(self, quippy._quippy.f90wrap_..._finalise, self._handle)
    """
    # Pattern matches finalizer lines that incorrectly use _initialise
    # We need to be careful to only match finalizer lines, not other uses of _initialise
    # The pattern looks for: self._finalizer = weakref.finalize(..._initialise...)
    pattern = r'(self\._finalizer\s*=\s*weakref\.finalize\([^)]*?)_initialise(\s*,|\s*\))'

    def replacer(match):
        return match.group(1) + '_finalise' + match.group(2)

    return re.sub(pattern, replacer, content)


def fix_init_alloc(content):
    """Fix __init__ methods to set self._alloc = True in else clause."""
    # Pattern: else clause that gets handle but doesn't set _alloc
    pattern = r'(else:\s+result = \w+\._quippy\.f90wrap_\w+\(\)\s+self\._handle = result\[0\] if isinstance\(result, tuple\) else result)\n'
    replacement = r'\1\n            self._alloc = True\n'
    return re.sub(pattern, replacement, content)


def fix_overloaded_interface(content):
    """Fix overloaded interface pattern to avoid method shadowing."""

    # Pattern to match overloaded interface definitions
    # Looks for: for proc in [ClassName.method1, ClassName.method2]:
    pattern = r'(\s+)(def (\w+)\(\*args, \*\*kwargs\):.*?\n\s+""".*?""")\n(\s+)(for proc in \[([^\]]+)\]:)'
    
    def replacer(match):
        indent = match.group(1)
        method_def = match.group(2)
        method_name = match.group(3)
        for_indent = match.group(4)
        for_line = match.group(5)
        proc_list = match.group(6)
        
        # Parse the proc list to find methods that will be shadowed
        procs = [p.strip() for p in proc_list.split(',')]
        
        # Check if any proc has the same name as the method being defined
        shadowed = []
        new_procs = []
        for i, proc in enumerate(procs):
            # Extract method name from expressions like "ClassName.method_name"
            if '.' in proc:
                proc_method = proc.split('.')[-1]
                if proc_method == method_name:
                    # This will be shadowed - create a saved reference
                    saved_name = f"_{method_name}_{i}"
                    shadowed.append((proc_method, saved_name))
                    # Replace in proc list
                    new_procs.append(proc.rsplit('.', 1)[0] + '.' + saved_name)
                else:
                    new_procs.append(proc)
            else:
                new_procs.append(proc)
        
        if not shadowed:
            # No shadowing, return original
            return match.group(0)
        
        # Generate the fix
        save_refs = '\n'.join([f"{indent}# Save references to original methods before overloading"] + 
                              [f"{indent}{saved} = {orig}" for orig, saved in shadowed] +
                              [f"{indent}"])
        
        new_proc_list = ', '.join(new_procs)
        
        return f"{save_refs}\n{indent}{method_def}\n{for_indent}{for_line.replace(proc_list, new_proc_list)}"
    
    return re.sub(pattern, replacer, content, flags=re.DOTALL)


def patch_file(filepath):
    """Patch a single Python file."""
    print(f"Patching {filepath}...")

    content = filepath.read_text()
    original_content = content

    # Fix finalizer bug where _initialise is incorrectly used instead of _finalise
    content = fix_finalizer_initialise_bug(content)

    # Note: fix_init_alloc is now handled by f90wrap itself (pywrapgen.py)
    # We only need to fix overloaded interfaces here
    content = fix_overloaded_interface(content)

    if content != original_content:
        filepath.write_text(content)
        print(f"  ✓ Patched {filepath}")
        return True
    else:
        print(f"  - No changes needed for {filepath}")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: patch_f90wrap_interfaces.py <file1.py> [file2.py ...]")
        sys.exit(1)
    
    patched_count = 0
    for arg in sys.argv[1:]:
        filepath = Path(arg)
        if filepath.exists():
            if patch_file(filepath):
                patched_count += 1
        else:
            print(f"Warning: {filepath} not found")
    
    print(f"\nPatched {patched_count} file(s)")


if __name__ == '__main__':
    main()
