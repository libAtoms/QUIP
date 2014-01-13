#!/usr/bin/env python

import socket
import threading
import SocketServer
import Queue

from quippy.atoms import Atoms
from quippy.clusters import create_hybrid_weights, create_cluster_simple
from quippy.potential import Potential
from quippy.structures import diamond, supercell
from quippy.util import args_str

MSG_LEN_SIZE = 8
MSG_END_MARKER = 'done.\n'
MSG_END_MARKER_SIZE = len(MSG_END_MARKER)

QUIPPY_HOST = ''
QUIPPY_PORT = 8888

N_JOBS = 1

class QuippyRequestHandler(SocketServer.StreamRequestHandler):

    def handle(self):
        ip, port = self.client_address

        # receive request code and client ID
        request_str = self.rfile.read(MSG_LEN_SIZE)
        request = request_str[0]
        client_id = int(request_str[1:])
        
        if client_id > N_JOBS-1:
            raise RuntimeError('Unknown client ID %d outside of range 0 < ID < %d' % (client_id, N_JOBS-1))

        print 'Serving "%s" request from %s:%d client %d' % (request, ip, port, client_id)

        if request == 'A':
            # client is ready for some more work to do
            print 'Queue size: %d' % input_qs[client_id].qsize()
            at = input_qs[client_id].get()
            data = at.write('string').encode('ascii')
            data_length = ('%8d' % len(data)).encode('ascii')
            self.wfile.write(data_length)
            self.wfile.write(data)
            print 'Sent Atoms to client %s' % client_id

        elif request == 'R':
            # results are available from client
            data_size = int(self.rfile.read(MSG_LEN_SIZE))
            data = self.rfile.read(data_size)
            at = Atoms(data, format='string')
            print 'Got Atoms from client %s, energy=%.3f' % (client_id, at.energy)
            input_qs[client_id].task_done()
            output_qs[client_id].put(at)
                        
        else:
            raise RuntimeError('Unknown request code "%s"' % request)

        # say goodbye to this client
        self.wfile.write(MSG_END_MARKER)
        
class ThreadedQuippyServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
    pass

server = ThreadedQuippyServer((QUIPPY_HOST, QUIPPY_PORT), QuippyRequestHandler)

server_thread = threading.Thread(target=server.serve_forever)
server_thread.daemon = True
print 'Starting threaded quippy server on %s:%d...' % (QUIPPY_HOST, QUIPPY_PORT)
server_thread.start()

input_qs = [Queue.Queue() for i in range(N_JOBS)]
output_qs = [Queue.Queue() for i in range(N_JOBS)]

def iter_little_clusters(at, **cluster_args):
    saved_hybrid_mark = None
    if hasattr(at, 'hybrid_mark'):
        saved_hybrid_mark = at.properties['hybrid_mark'].copy()
        indices = saved_hybrid_mark.nonzero()[0]
    else:
        indices = at.indices
    
    at.add_property('hybrid_mark', 0, overwrite=True)

    for i in indices:
        at.hybrid_mark[:] = 0
        at.hybrid_mark[i] = True
        create_hybrid_weights(at, args_str=args_str(cluster_args))
        c = create_cluster_simple(at,
                                  args_str=args_str(cluster_args),
                                  mark_name='hybrid_mark')
        yield c

    if saved_hybrid_mark is not None:
        at.add_property('hybrid_mark', saved_hybrid_mark, overwrite=True)
    else:
        del at.properties['hybrid_mark']


d = diamond(5.44, 14)
at = supercell(d, 2, 2, 2)
at.rattle(0.05)
at.calc_connect()
clusters = list(iter_little_clusters(at, buffer_hops=3))
        
for i, c in enumerate(clusters):
    input_qs[i % N_JOBS].put(c)

for input_q in input_qs:
    input_q.join()


print 'Shutting down server'
server.shutdown()



