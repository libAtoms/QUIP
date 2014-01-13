#!/usr/bin/env python

import sys
import socket
import threading
import SocketServer
import Queue

from quippy.atoms import Atoms
from quippy.lotf import iter_atom_centered_clusters
from quippy.potential import Potential
from quippy.structures import diamond, supercell

MSG_LEN_SIZE = 8
MSG_END_MARKER = 'done.\n'
MSG_END_MARKER_SIZE = len(MSG_END_MARKER)

QUIPPY_HOST = ''
QUIPPY_PORT = 8888

N_JOBS = int(sys.argv[1])

class QuippyRequestHandler(SocketServer.StreamRequestHandler):

    def handle(self):
        ip, port = self.client_address

        # receive request code and client ID
        request_str = self.rfile.read(MSG_LEN_SIZE)
        request = request_str[0]
        client_id = int(request_str[1:])
        
        if client_id > N_JOBS-1:
            raise RuntimeError('Unknown client ID %d outside of range 0 < ID < %d' % (client_id, N_JOBS-1))

        print '"%s" request from %s:%d client %d qsize %d' % (request, ip, port, client_id, input_qs[client_id].qsize())

        if request == 'A':
            # client is ready for some more work to do
            at = input_qs[client_id].get()
            data = at.write('string').encode('ascii')
            data_length = ('%8d' % len(data)).encode('ascii')
            self.wfile.write(data_length)
            self.wfile.write(data)

        elif request == 'R':
            # results are available from client
            data_size = int(self.rfile.read(MSG_LEN_SIZE))
            data = self.rfile.read(data_size)
            at = Atoms(data, format='string')
            input_qs[client_id].task_done()
            output_qs[client_id].put(at.force[:, 1].copy())
                        
        else:
            raise RuntimeError('Unknown request code "%s"' % request)

        # say goodbye to this client
        self.wfile.write(MSG_END_MARKER)
        
class ThreadedQuippyServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
    pass

server = ThreadedQuippyServer((QUIPPY_HOST, QUIPPY_PORT), QuippyRequestHandler)

server_thread = threading.Thread(target=server.serve_forever)
server_thread.daemon = True
print 'Starting threaded quippy server on %s:%d with N_JOBS=%d' % (QUIPPY_HOST, QUIPPY_PORT, N_JOBS)
server_thread.start()

# we need an input and an output Queue for each client: this is so
# that we can exploit wavefunction reuse by sending consecutive
# clusters belonging to the same atom to the same QM partition
input_qs = [Queue.Queue() for i in range(N_JOBS)]
output_qs = [Queue.Queue() for i in range(N_JOBS)]

# setup clusters to be calculated
d = diamond(5.44, 14)
at = supercell(d, 2, 2, 2)
at.rattle(0.05)
at.calc_connect()
n_clusters = 0
for i, c in enumerate(iter_atom_centered_clusters(at, buffer_hops=3)):
    n_clusters += 1
    input_qs[i % N_JOBS].put(c)

# wait for all clients to finish processing clusters
for input_q in input_qs:
    input_q.join()

# collect results
result_at = at.copy()
result_at.add_property('cluster_force', 0.0, n_cols=3)
for i in range(n_clusters):
    result = output_qs[i % N_JOBS].get()
    result_at.cluster_force[:, i+1] = result
    output_qs[i % N_JOBS].task_done()
    
# reference calculation of full system
pot = Potential('IP SW')
pot.calc(result_at, args_str='force')

# compare cluster and reference forces
max_force_err = abs(result_at.force - result_at.cluster_force).max()
print 'max_force_error = %.3e eV/A' % max_force_err
assert max_force_err < 1.0e-6

print 'Shutting down quippy server'
server.shutdown()



