import os
from pathlib import Path
from subprocess import Popen, PIPE

def test_shell_scripts(tmpdir):
    os.chdir(tmpdir)
    for f in Path(__file__).parent.glob('test_*.sh'):
        print(f)

        p = Popen([f], stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        if p.returncode != 0:
            raise RuntimeError(f'Shell script test {f} failed with error code {p.returncode} stderr {stderr.decode()} stdout {stdout.decode()}')
