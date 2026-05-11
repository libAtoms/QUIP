#!/usr/bin/env python3
"""
Regression test for https://github.com/libAtoms/QUIP/issues/724

Invalid XML passed to quippy.potential.Potential(param_filename=...) must raise
a catchable Python exception. Previously, FoX's error handler called Fortran
`stop` directly, silently killing the Python process before any try/except
could run.

The test runs the user's reproducer in a subprocess because, in the RED state,
the bug terminates the interpreter; an in-process assertRaises would crash the
whole test runner instead of failing one test.
"""

import os
import subprocess
import sys
import tempfile
import textwrap
import unittest


class TestFoxXmlErrors(unittest.TestCase):

    def _run_reproducer(self, args_str):
        """Run a tiny quippy script in a subprocess and return its result.

        The script writes a non-XML file, then tries to load it as a Potential
        param_filename. It prints CAUGHT:<exc-type> on a Python exception,
        NO_EXCEPTION if Potential() returned without raising.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            bad_xml = os.path.join(tmpdir, "bad.xml")
            with open(bad_xml, "w") as fh:
                fh.write("pair_style lj/cut 3.5\n")

            code = textwrap.dedent(f"""
                import sys
                from quippy.potential import Potential
                try:
                    pot = Potential({args_str!r}, param_filename={bad_xml!r})
                    print("NO_EXCEPTION")
                except BaseException as e:
                    print("CAUGHT:" + type(e).__name__ + ":" + str(e)[:200])
                    sys.exit(0)
            """)

            return subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True,
                text=True,
                timeout=120,
            )

    def test_invalid_xml_raises_python_exception_sw(self):
        """IP SW with a non-XML param file must raise, not silently exit."""
        result = self._run_reproducer("IP SW")
        self.assertIn(
            "CAUGHT:",
            result.stdout,
            msg=(
                "Expected Python exception from invalid XML.\n"
                f"returncode={result.returncode}\n"
                f"stdout={result.stdout!r}\n"
                f"stderr={result.stderr!r}"
            ),
        )
        self.assertNotIn("NO_EXCEPTION", result.stdout)


if __name__ == "__main__":
    unittest.main()
