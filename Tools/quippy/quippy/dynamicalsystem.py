# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from quippy import FortranDynamicalSystem

class DynamicalSystem(FortranDynamicalSystem):

    __doc__ = FortranDynamicalSystem.__doc__

    def run(self, pot, dt, n_steps, hook_interval=None, write_interval=None, connect_interval=None, trajectory=None, args_str=None, hook=None,
            save_interval=None):
        if hook is None:
            if hook_interval is not None:
                raise ValueError('hook_interval not permitted when hook is not present. save_interval is used instead')
            traj = []
            FortranDynamicalSystem.run(self, pot, dt, n_steps, lambda:traj.append(self.atoms.copy()), hook_interval=save_interval, write_interval=write_interval,
                                       connect_interval=connect_interval, trajectory=trajectory, args_str=args_str)
            return traj
        else:
            FortranDynamicalSystem.run(self, pot, dt, n_steps, hook, hook_interval, write_interval,
                                       connect_interval, trajectory, args_str)
