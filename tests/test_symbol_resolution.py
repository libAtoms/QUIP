#!/usr/bin/env python3
"""
Test to verify f90wrap_abort symbol resolution works correctly.

This test ensures that:
1. The weak symbol in libAtoms doesn't conflict with the real implementation
2. Error handling works correctly (errors are raised, not silently ignored)
3. Symbol resolution is deterministic (not flaky)
"""

import unittest
import subprocess
import sys
import os
import platform


class TestSymbolResolution(unittest.TestCase):
    """Test that f90wrap_abort symbol resolution works correctly"""

    def _find_libatoms(self, lib_extension):
        """Find libAtoms library file."""
        import quippy
        quippy_dir = os.path.dirname(quippy.__file__)

        potential_paths = [
            os.path.join(quippy_dir, '..', '..', '..', 'builddir', 'src', 'libAtoms', f'liblibAtoms{lib_extension}'),
            os.path.join(quippy_dir, '..', '..', 'builddir', 'src', 'libAtoms', f'liblibAtoms{lib_extension}'),
        ]

        for path in potential_paths:
            abs_path = os.path.abspath(path)
            if os.path.exists(abs_path):
                return abs_path
        return None

    def _check_weak_symbol(self, lib_extension, nm_flags, weak_check_func):
        """Common implementation for weak symbol verification."""
        libatoms_path = self._find_libatoms(lib_extension)
        if libatoms_path is None:
            self.skipTest(f"Could not find liblibAtoms{lib_extension}")

        try:
            result = subprocess.run(
                ['nm'] + nm_flags + [libatoms_path],
                capture_output=True,
                text=True,
                check=True
            )

            abort_lines = [line for line in result.stdout.split('\n')
                          if 'f90wrap_abort' in line]

            self.assertTrue(len(abort_lines) > 0,
                          "f90wrap_abort symbol not found in libAtoms")

            weak_check_func(abort_lines)

        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.skipTest(f"Could not run nm: {e}")

    def test_repeated_calculations_are_deterministic(self):
        """
        Run the same calculation multiple times and verify results are identical.
        This would catch flakiness from symbol resolution issues.
        """
        from quippy.potential import Potential
        from ase.io import read
        import numpy as np

        # Load test data
        test_dir = os.path.dirname(os.path.abspath(__file__))
        pot = Potential("IP GAP", param_filename=os.path.join(test_dir, "GAP.xml"))
        atoms = read(os.path.join(test_dir, 'gap_sample.xyz'))
        atoms.calc = pot

        # Run calculation 5 times
        energies = []
        forces_list = []

        for i in range(5):
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces().copy()
            energies.append(energy)
            forces_list.append(forces)

        # All energies should be identical (not just close, but exact)
        for i in range(1, 5):
            self.assertEqual(energies[0], energies[i],
                           f"Energy in run {i+1} differs from run 1: "
                           f"{energies[i]} != {energies[0]}")

        # All forces should be identical
        for i in range(1, 5):
            np.testing.assert_array_equal(
                forces_list[0], forces_list[i],
                err_msg=f"Forces in run {i+1} differ from run 1"
            )

    def test_error_handling_raises_exceptions(self):
        """
        Test that errors in Fortran code properly propagate to Python.
        If f90wrap_abort stub is used instead of real implementation,
        errors might be silently ignored.
        """
        from quippy.potential import Potential

        # This should raise an exception, not crash or silently fail
        with self.assertRaises(Exception) as cm:
            # Try to load non-existent file
            pot = Potential('IP GAP', param_filename='/tmp/nonexistent_file_12345.xml')

        # Verify we got a meaningful exception (not a crash or silent failure)
        self.assertTrue(len(str(cm.exception)) > 0,
                       "Exception message should not be empty")

    @unittest.skipUnless(platform.system() == 'Darwin', "macOS-specific test")
    def test_weak_symbol_on_macos(self):
        """On macOS, verify f90wrap_abort is a weak symbol."""
        def check_weak(lines):
            has_weak = any('weak' in line.lower() for line in lines)
            self.assertTrue(has_weak,
                          f"f90wrap_abort should be weak symbol on macOS.\n"
                          f"Symbol info: {lines}")

        self._check_weak_symbol('.dylib', ['-m'], check_weak)

    @unittest.skipUnless(platform.system() == 'Linux', "Linux-specific test")
    def test_weak_symbol_on_linux(self):
        """On Linux, verify f90wrap_abort is a weak symbol."""
        def check_weak(lines):
            has_weak = any(' W ' in line for line in lines)
            has_strong = any(' T ' in line for line in lines)
            self.assertTrue(has_weak,
                          f"f90wrap_abort should be weak symbol (W) on Linux.\n"
                          f"Symbol info: {lines}")
            self.assertFalse(has_strong,
                           f"f90wrap_abort should NOT be strong symbol (T) on Linux.\n"
                           f"Symbol info: {lines}")

        self._check_weak_symbol('.so', ['-D'], check_weak)


if __name__ == '__main__':
    unittest.main()
