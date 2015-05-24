import sys
import os
import optparse

import numpy as np
import matplotlib.pyplot as plt

from quippy.io import AtomsReader
from quippy.farray import farray, fzeros
from quippy.ringstat import distance_map, count_sp_rings
from quippy.structures import supercell
from quippy.util import parse_slice

def tetrahedra_to_bonds(q):
    si_si = []
    si_o = []

    max_si_si_dist = 0.0
    
    for i in q.indices:
        if q.z[i] != 14:
            continue

        for nj in q.neighbours[i]:
            j = nj.j

            if q.z[j] != 8:
                continue

            # i and j are an Si-O neighbour pair
            si_o.append((i, j, tuple(nj.shift), nj.distance))

            for nk in q.neighbours[j]:
                k = nk.j
                
                if not i < k:
                    continue

                shift = nj.shift + nk.shift

                if q.z[k] != 14 or (k == i and all(shift == 0)):
                    continue

                # i and k are an Si-Si pair with a common O neighbour
                r_ij = float(np.sqrt(((nj.diff + nk.diff)**2).sum()))
                if r_ij > max_si_si_dist:
                    max_si_si_dist = r_ij
                    
                si_si.append((i, k, tuple(shift), r_ij))

    print '%d Si-Si, %d Si-O' % (len(si_si), len(si_o))

    q.connect.wipe()
    for (i, j, shift, r_ij) in si_si:
        q.connect.add_bond(q.pos, q.lattice, i, j, shift)

    return si_si, si_o, max_si_si_dist

p = optparse.OptionParser(usage='%prog: [OPTIONS] INFILE...')

p.add_option('-c', '--si-si-cutoff', action='store', type=float, help='Si-Si cutoff. Default is longest Si-Si distance found.')
p.add_option('-C', '--si-o-cutoff', action='store', type=float, default=2.0, help='Si-O cutoff. Default is 2.0 A.')
p.add_option('-m', '--max-ring-size', action='store', type=int, default=12, help='Maximum ring size to search for. Default 12.')
p.add_option('-r', '--range', action='store', help="""Range of frames to include. Should be either a single frame
number or a slice [start]:[stop][:step]. If -r is omitted,
default is all frames. Frames start from zero. Negative indices count backwards from the end of the file,
with -1 being the last frame.""")
p.add_option('-s', '--supercell', action='store', type=int, nargs=3, help="Make a supercell from input structure before computing ring statistics.")
p.add_option('-t', '--tetra', action='store_true', help="Convert Si-O-Si bonds in tetrahedra to Si-Si bonds.")
p.add_option('-p', '--plot', action='store_true')
p.add_option('-n', '--no-clear', action='store_true', help="Do not clear figure before plotting. Useful for comparisons")

opt, infiles = p.parse_args()

# check for proper format of --range
if opt.range is not None:
    try:
        opt.range = parse_slice(opt.range)
    except:
        p.error('Cannot parse range "%s" - should be in format [start]:[stop][:step]' % opt.range)
else:
    # Default is all frames
    opt.range = slice(0, None, None)

if isinstance(opt.range, int):
    if opt.range >= 0:
        opt.range = slice(opt.range, opt.range+1,+1)
    else:
        opt.range = slice(opt.range, opt.range-1,-1)    

for infile in infiles:
    basename = os.path.basename(os.path.splitext(infile)[0])

    configs = AtomsReader(infile, start=opt.range.start, stop=opt.range.stop, step=opt.range.step)

    for frame, q0 in enumerate(configs):
        print 'Read %d atoms from file %s frame %d' % (q0.n, infile, frame)

        diameter = farray(0, dtype=np.int32)
        ring_sizes = range(1,opt.max_ring_size+1)

        if opt.supercell is None:
            q = q0
        else:
            print 'Forming %d x %d x %d supercell' % tuple(opt.supercell)
            q = supercell(q0, opt.supercell[0], opt.supercell[1], opt.supercell[2])
        q.map_into_cell()

        q.set_cutoff(opt.si_o_cutoff)
        q.calc_connect()

        if opt.tetra:
            print 'Converting topology from Si-O-Si to Si-Si'
            si_si, si_o, si_si_cutoff = tetrahedra_to_bonds(q)

        if opt.si_si_cutoff is None:
            opt.si_si_cutoff = si_si_cutoff

        dm = distance_map(q, q.n, q.n, diameter=diameter)
        print 'Distance map diameter %d' % diameter
        print 'Max ring size %d' % opt.max_ring_size

        assert(diameter >= opt.max_ring_size)

        print 'Using Si-Si cutoff of %.3f' % opt.si_si_cutoff

        ring_counts = fzeros(opt.max_ring_size, dtype=np.int32)
        count_sp_rings(q, opt.si_si_cutoff, dm, opt.max_ring_size, ring_counts)

        rings_per_si = np.array(ring_counts.astype(float)/(q.z == 14).sum())

        # do it again saving the rings
        n_rings = ring_counts.sum()
        ring_counts[:] = 0
        rings_array = fzeros((opt.max_ring_size+1, n_rings), dtype=np.int32)
        dr = fzeros((3,opt.max_ring_size, n_rings))
        count_sp_rings(q, opt.si_si_cutoff, dm, opt.max_ring_size, ring_counts, rings_out=rings_array, dr_out=dr)

        ring_lens = []
        rings = []
        ring_centers = []

        for ring_len, ring in zip(rings_array[1,:], rings_array[2:,:]):
            ring_lens.append(ring_len)
            rings.append(list(ring[:ring_len]))
            ring_center = (q.pos[:,ring[ring_len]] +  dr[:,:ring_len,1].T).mean(axis=1)
            ring_centers.append(ring_center)

        ring_lens = np.array(ring_lens)
        ring_centers = np.array(ring_centers)

        print
        title = '%-10s %10s %10s' % ('Ring Size', 'Count', 'Rings per Si')
        print title
        print '-'*len(title)
        for size, count, frac in zip(ring_sizes, ring_counts, rings_per_si):
            print '%-10d %10d %10.3f' % (size, count, frac)
        print '-'*len(title)
        print '%-10s %10d %10.3f' % ('Total', ring_counts.sum(), rings_per_si.sum())
        print


        if len(configs) > 1:
            outfile = basename+'.%05d.out' % frame
            figfile = basename+'.%05d.pdf' % frame
        else:
            outfile = basename+'.out'
            figfile = basename+'.pdf'

        print 'Saving results to %s' % outfile
        np.savetxt(outfile, np.c_[rings_array.T, ring_centers])

        if opt.plot:
            if not opt.no_clear:
                plt.clf()
            plt.plot(ring_sizes, rings_per_si, 'o-', label=basename)
            plt.xlabel('Ring size')
            plt.ylabel('Rings per Si')
            plt.draw()
            plt.savefig(figfile)
