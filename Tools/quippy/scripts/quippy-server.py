#!/usr/bin/env python

import sys
import socket
import subprocess
import multiprocessing
import threading
import SocketServer
from Queue import Queue

import numpy as np

from quippy.atoms import Atoms
from quippy.farray import fzeros, farray
from quippy.lotf import iter_atom_centered_clusters
from quippy.potential import Potential
from quippy.structures import diamond, supercell

MSG_LEN_SIZE = 8
MSG_END_MARKER = 'done.\n'
MSG_END_MARKER_SIZE = len(MSG_END_MARKER)
MSG_INT_SIZE = 4
MSG_FLOAT_SIZE = 25
MSG_FLOAT_FORMAT = '%25.16f'
MSG_INT_FORMAT = '%4d'

N_JOBS = int(sys.argv[1])

def pack_atoms_to_reftraj_str(at):
    data = ''
    data += MSG_INT_FORMAT % at.n + '\n'
    data += ''.join(MSG_INT_FORMAT % z for z in at.z) + '\n'
    for i in (1, 2, 3):
        data += (3*MSG_FLOAT_FORMAT) % tuple(at.lattice[i, :]) + '\n'
    for i in at.indices:
        data += (3*MSG_FLOAT_FORMAT) % tuple(np.dot(at.g, at.pos[:, i])) + '\n'

    est_data_len = (at.n + 1)*MSG_INT_SIZE + MSG_FLOAT_SIZE*(3*3 + 3*at.n) + at.n+5

    # preceed message by its length
    data_length = ('%8d' % len(data)).encode('ascii')
    data = data_length + data.encode('ascii')
    return data

def unpack_reftraj_str_to_atoms(data):
    lines = data.split('\n')
    z = [int(field) for field in lines[0].split()]
    at = Atoms(n=len(z), lattice=np.eye(3))
    at.set_atoms(z)
    for i in (1, 2, 3):
        at.lattice[i, :] = [float(x) for x in lines[i].split()]
    at.set_lattice(at.lattice)
    for i, line in zip(at.indices, lines[4:]):
        t = [float(x) for x in line.split()]
        at.pos[:, i] = np.dot(at.lattice, t)
    return at

def pack_results_to_reftraj_output_str(at):
    data = ''
    data += MSG_INT_FORMAT % at.n + '\n'
    data += MSG_FLOAT_FORMAT % at.energy + '\n'
    for i in at.indices:
        data += (3*MSG_FLOAT_FORMAT) % tuple(at.force[:, i]) + '\n'
    data += (6*MSG_FLOAT_FORMAT) % (at.virial[1,1], at.virial[2,2], at.virial[3,3], at.virial[1,2], at.virial[2,3], at.virial[1,3])

    # preceed message by its length
    data_length = ('%8s' % len(data)).encode('ascii')
    data = data_length + data
    return data

def unpack_reftraj_output_str_to_results(data):
    lines = data.strip().split('\n')
    natoms = int(lines[0])
    energy = float(lines[1])
    force = farray(np.loadtxt(lines[2:-1])).T
    v6 = [float(v) for v in lines[-1].split()]
    virial = fzeros((3,3))
    virial[1,1], virial[2,2], virial[3,3], virial[1,2], virial[2,3], virial[1,3] = v6
    virial[2,1] = virial[1,2]
    virial[3,2] = virial[2,3]
    virial[3,1] = virial[1,3]
    return (natoms, energy, force, virial)


class QuippyRequestHandler(SocketServer.StreamRequestHandler):

    def handle(self):
        ip, port = self.client_address

        # receive request code and client ID
        request_str = self.rfile.read(MSG_LEN_SIZE)
        request = request_str[0]
        client_id = int(request_str[1:])
        
        if client_id > N_JOBS-1:
            raise RuntimeError('Unknown client ID %d outside of range 0 < ID < %d' % (client_id, N_JOBS-1))

        print '"%s" request from %s:%d client %d' % (request, ip, port, client_id)
        print 'input queue lengths ', ''.join(['%d:%d ' % (i,q.qsize()) for (i,q) in enumerate(input_qs)])
        print 'output queue lengths ', ''.join(['%d:%d ' % (i,q.qsize()) for (i,q) in enumerate(output_qs)])

        if request == 'A':
            # client is ready for some more work to do
            data = input_qs[client_id].get()
            self.wfile.write(data)

        elif request == 'R':
            # results are available from client
            data_size = int(self.rfile.read(MSG_LEN_SIZE))
            data = self.rfile.read(data_size)
            output_qs[client_id].put(data)
            input_qs[client_id].task_done()
                                    
        else:
            raise RuntimeError('Unknown request code "%s"' % request)

        # say goodbye to this client
        self.wfile.write(MSG_END_MARKER)
        print 'done processing "%s" request from client %d' % (request, client_id)
        
class ThreadedQuippyServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
    request_queuesize = 4*N_JOBS
    allow_reuse_address = True

server = ThreadedQuippyServer(('', 8888), QuippyRequestHandler)
ip, port = server.server_address

server_thread = threading.Thread(target=server.serve_forever)
server_thread.daemon = True
print 'Starting threaded quippy server on %s:%d with N_JOBS=%d' % (ip, port, N_JOBS)
server_thread.start()

# we need an input and an output Queue for each client: this is so
# that we can exploit wavefunction reuse by sending consecutive
# clusters belonging to the same atom to the same QM partition
input_qs = [Queue() for i in range(N_JOBS)]
output_qs = [Queue() for i in range(N_JOBS)]
cluster_list = [ [] for i in range(N_JOBS)]

# setup clusters to be calculated
d = diamond(5.44, 14)
at = supercell(d, 2, 2, 2)
at.rattle(0.05)
at.calc_connect()
clusters = list(iter_atom_centered_clusters(at, buffer_hops=2, randomise_buffer=False))

pot = Potential('IP SW')

for i, c in enumerate(clusters):
    data = pack_atoms_to_reftraj_str(c)
    qi = i % N_JOBS
    cluster_list[qi].append(c)
    input_qs[qi].put(data)
    pot.calc(c, args_str='energy force virial')

# write initial input file for each client
for i, (input_q, cs) in enumerate(zip(input_qs, cluster_list)):
    data = input_q.get()
    c = cs[0]
    c.write('atoms.%d.xyz' % i)

# spawn clients as background processes
clients = []
stdouts = []
for i in range(N_JOBS):
    stdout = open('client.stdout.%d' % i, 'w')
    client = subprocess.Popen(['/Users/jameskermode/Code/QUIP/build.darwin_x86_64_gfortran_openmp/socktest', ip, str(port), str(i)],
                              stderr=stdout)
    clients.append(client)

print 'All calculations queued, waiting for results.'

# wait for input queues to empty
for input_q in input_qs:
    input_q.join()

print 'Input queues drained. Shutting down clients.'

# stop the clients by sending them a calculation with zero atoms
dummy_at = Atoms(n=0, lattice=np.eye(3))
dummy_data = pack_atoms_to_reftraj_str(dummy_at)
for input_q in input_qs:
    input_q.put(dummy_data)

# wait for them all to shutdown
for client in clients:
    client.wait()

print 'Clients terminated.'   

for stdout in stdouts:
    stdout.flush()
    stdout.close()

print 'Client logs flushed.'

print 'Collecting results.'

print 'Output queue lengths:', [q.qsize() for q in output_qs]

# collect results
result_at = at.copy()
result_at.add_property('cluster_force', 0.0, n_cols=3)
for i, c in enumerate(clusters):
    data = output_qs[i % N_JOBS].get()
    natoms, energy, force, virial = unpack_reftraj_output_str_to_results(data)

    # compare returned values to the locally-computed energy, force and virial
    assert natoms == c.n
    assert abs(energy - c.energy) < 1e-7
    assert abs(force - c.force).max() < 1e-6
    assert abs(virial - c.virial).max() < 1e-6

    # we only need force on first atom
    result_at.cluster_force[:, i+1] = force[:, 1]
    output_qs[i % N_JOBS].task_done()

print 'Output queues drained.'
    
# reference calculation of full system
pot.calc(result_at, args_str='force')

# compare cluster and reference forces
max_force_err = abs(result_at.force - result_at.cluster_force).max()
print 'max_force_error = %.3e eV/A' % max_force_err
assert max_force_err < 1.0e-6

print 'Shutting down quippy server.'
server.shutdown()



