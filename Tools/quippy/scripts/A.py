#!/usr/bin/env python

from qlab import *
import optparse

p = optparse.OptionParser(usage='%prog [options] <trajectory file>...')

p.add_option('-f', '--frame', action='append', help='Initial frame to show. Can be given separately for each trajectory file. Default is last frame in file.', type=int)
p.add_option('-p', '--property', action='append', help='Property to show. Can be given separately for each trajectory file')
p.add_option('-a', '--arrows', action='append', help='Property to use to draw arrows. Can be given separately for each trajectory file')
p.add_option('-l', '--load-view', action='store', help='Load view from AtomEye command script')
p.add_option('-W', '--width', action='store', help="""Width of output movie, in pixels.""", type='int')
p.add_option('-H', '--height', action='store', help="""Height of output movie, in pixels.""", type='int')
p.add_option('-A', '--aspect', action='store', help="""Aspect ratio. Used if only one of --width or --height is given. Default 0.75.""", default=0.75, type='float')
p.add_option('-R', '--rcut', action='append', help="""Following three arguments should be SYM1 SYM2 INCREMENT, to increment cutoff distance for SYM1-SYM2 bonds.""", nargs=3)
p.add_option('-s', '--single', action='append', help="""Read next argument as a single frame, not a trajectory""", default=[])
p.add_option('-L', '--loadall', action='store_true', help="""Load all frames from all trajectories into memory on startup.""")

opt, args = p.parse_args()

viewers = []
for filename in opt.single:
    viewers.append(atoms(filename), single=True)
for filename in args:
    viewers.append(atoms(filename), loadall=opt.loadall)

if opt.frame is None:
    opt.frame = [-1 for viewer in viewers]

show_args_list = [{} for viewer in viewers]
for arg in ['frame', 'property', 'arrows']:
    values = getattr(opt, arg)
    if values is None:
        continue
    if len(values) == 1:
        values = [values[0] for traj in viewers]

    if len(values) != len(viewers):
        p.error('Number of -%s/--%s options does not match number of trajectory files' % (arg[0], arg))

    for show_args, value in zip(show_args_list, values):
        show_args[arg] = value

for traj, show_args in zip(viewers, show_args_list):
    if opt.load_view is not None:
        print 'Loading view script %s' % opt.load_view
        traj.load_script(opt.load_view)

    if opt.rcut is not None:
        print 'Applying rcut patches %r' % opt.rcut
        for (sym1, sym2, rcut) in opt.rcut:
            traj.rcut_patch(sym1, sym2, float(rcut))

    if opt.width is not None or opt.height is not None:
        if opt.width  is None: opt.width = int(opt.height/opt.aspect)
        if opt.height is None: opt.height = int(opt.width*opt.aspect)
        print 'Setting viewer size to %d x %d pixels' % (opt.width, opt.height)
        traj.resize(opt.width, opt.height)

    if len(show_args) != 0:
        print 'Applying show_args=%r to trajectory %s' % (show_args, traj.name)
        traj.show(**show_args)

del traj, show_args, filename, viewers, show_args_list, p, opt, args, arg, optparse

from IPython import embed
embed()

