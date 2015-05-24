#!/usr/bin/env python

import multiprocessing
import sys
import socket
import threading
import SocketServer

import numpy as np

from quippy.atoms import Atoms
from quippy.potential import Potential

MSG_LEN_SIZE = 8
MSG_END_MARKER = 'done.'
MSG_END_MARKER_SIZE = len(MSG_END_MARKER)
MSG_FORMAT = '%16.10f'
MSG_INT_SIZE = 4
MSG_FLOAT_SIZE = 16


def pack_atoms_to_reftraj_str(at):
    data = ''
    data += ' '.join(str(z) for z in at.z) + '\n'
    for i in (1, 2, 3):
        data += (3*MSG_FORMAT) % tuple(at.lattice[i, :]) + '\n'
    for i in at.indices:
        data += (3*MSG_FORMAT) % tuple(np.dot(at.g, at.pos[:, i])) + '\n'

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
    data += MSG_FORMAT % at.energy + '\n'
    for i in at.indices:
        data += (3*MSG_FORMAT) % tuple(at.force[:, i]) + '\n'
    data += (6*MSG_FORMAT) % (at.virial[1,1], at.virial[2,2], at.virial[3,3], at.virial[1,2], at.virial[2,3], at.virial[1,3])
    data += '\n'

    print 'len(data) =', len(data), MSG_FLOAT_SIZE*(1 + 3*at.n + 6) + at.n+2
    
    # preceed message by its length
    data_length = ('%8s' % len(data)).encode('ascii')
    data = data_length + data
    return data

def unpack_reftraj_output_str_to_results(data):
    lines = data.split('\n')
    energy = float(lines[0])
    force = farray(np.loadtxt(lines[1:-1])).T
    v6 = [float(v) for v in lines[-1].split()]
    virial = fzeros((3,3))
    virial[1,1], virial[2,2], virial[3,3], virial[1,2], virial[2,3], virial[1,3] = v6
    virial[2,1] = virial[1,2]
    virial[3,2] = virial[2,3]
    virial[1,3] = virial[3,1]
    return (energy, force, virial)


def quippy_client(ip, port, client_id, maxsteps):
    pot = Potential('IP SW')
    args_str = 'force energy virial'

    print 'Starting client with ID %d maxsteps %d' % (client_id, maxsteps)
    
    for step in range(maxsteps):
        
        # open connection to server
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((ip, port))
        try:
            # say hello and ask server for some Atoms (status 'A')
            # identify ourselves with an ID number, formatted as an 8-byte ascii string
            id_str = ('A%7d' % client_id).encode('ascii')
            totalsent = 0
            while totalsent < MSG_LEN_SIZE:
                sent = sock.send(id_str[totalsent:])
                if sent == 0:
                    raise RuntimeError('socket connection broken while sending client ID')
                totalsent += sent

            # now receive the length of the data, formatted as an 8-byte ascii string
            data_length = ''
            while len(data_length) < MSG_LEN_SIZE:
                chunk = sock.recv(MSG_LEN_SIZE - len(data_length))
                if not chunk:
                    raise RuntimeError('socket connection broken while reading length')
                data_length += chunk
            data_length = int(data_length)

            # now receive the data itself
            data = ''
            while len(data) < data_length:
                chunk = sock.recv(data_length - len(data))
                if not chunk:
                    raise RuntimeError('socket connection broken while reading data')
                data += chunk

            # and finally receive the end marker
            marker = ''
            while len(marker) < MSG_END_MARKER_SIZE:
                chunk = sock.recv(MSG_END_MARKER_SIZE - len(marker))
                if not chunk:
                    raise RuntimeError('socket connection broken while reading end marker')
                marker += chunk

            assert marker == MSG_END_MARKER

        finally:
            sock.close()

        # construct Atoms object from data string in REFTRAJ format
        at = unpack_reftraj_str_to_atoms(data)

        # do the calculation - this will take a long time in real use-case
        pot.calc(at, args_str=args_str)

        print 'client %d calculation on %d atoms complete' % (client_id, len(at))

        # format the results into a string, cf. vasp_driver REJTRAJ_OUTPUT format
        data = pack_results_to_reftraj_output_str(at)

        # open a new connection to the server to send back the results
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((ip, port))
        try:
            # say hello and tell server we've got some results (status 'R')
            # identify ourselves with an ID number, formatted as an 8-byte ascii string
            id_str = ('R%7d' % client_id).encode('ascii')
            totalsent = 0
            while totalsent < MSG_LEN_SIZE:
                sent = sock.send(id_str[totalsent:])
                if sent == 0:
                    raise RuntimeError('socket connection broken while sending client ID')
                totalsent += sent

            # send back results of the calculation
            totalsent = 0
            while totalsent < len(data):
                sent = sock.send(data[totalsent:])
                if sent == 0:
                    raise RuntimeError('socket connection broken while sending results')
                totalsent += sent

            # wait to receive the end marker
            marker = ''
            while len(marker) < MSG_END_MARKER_SIZE:
                chunk = sock.recv(MSG_END_MARKER_SIZE - len(marker))
                if not chunk:
                    raise RuntimeError('socket connection broken while reading end marker')
                marker += chunk

            assert marker == MSG_END_MARKER

        finally:
            sock.close()
            


N_ATOMS = 64 # hard-coded for this test

ip = sys.argv[1]
port = int(sys.argv[2])
n_jobs = int(sys.argv[3])
client_id = int(sys.argv[4])
n_steps = N_ATOMS/n_jobs

quippy_client(ip, port, client_id, n_steps)
            
