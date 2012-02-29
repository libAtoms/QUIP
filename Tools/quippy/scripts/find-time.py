from quippy import *
import sys

traj = AtomsReader(sys.argv[1], range=(1,2)) # only read the first atom
time = float(sys.argv[2])

t_start = traj[0].time
t_end = traj[len(traj)-1].time

if time < t_start or time > t_end:
    sys.stderr.write('Time %f is out of range of this file, which contains %f <= t <= %f' % (time, t_start, t_end))

delta_t = (t_end - t_start)/len(traj)
frame = int((time - t_start)/delta_t)

print frame, traj[frame].time
