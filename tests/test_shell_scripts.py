import unittest
import os

from pathlib import Path
from subprocess import Popen, PIPE

import quippytest

class TestShellScripts(quippytest.QuippyTestCase):
    def test_shell_scripts(self):
        for f in Path(__file__).parent.glob('test_*.sh'):
            print("shell script test", f)

            p = Popen([f], stdout=PIPE, stderr=PIPE)
            stdout, stderr = p.communicate()

            if 'QUIP_WHEEL_TEST' in os.env and (p.returncode == 2 and f.endswith('test_xyz_6vector.sh')):
                print(f'Skipping {f} when testing built wheels')
                continue

            self.assertEqual(p.returncode, 0, f'Shell script test {f} failed with error code '
                                              f'{p.returncode} stderr {stderr.decode()} stdout {stdout.decode()}')


if __name__ == '__main__':
    unittest.main()
