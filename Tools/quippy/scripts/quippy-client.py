#!/usr/bin/env python

import sys
import socket
import threading
import SocketServer

from quippy.atoms import Atoms
from quippy.potential import Potential

MSG_LEN_SIZE = 8
MSG_END_MARKER = 'done.\n'
MSG_END_MARKER_SIZE = len(MSG_END_MARKER)

ATOMS_HOST = ''
ATOMS_PORT = 8888

CLIENT_ID = int(sys.argv[1])

pot = Potential('IP SW')
args_str = 'force energy virial'

for step in range(64):
    
    # open connection to server
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect((ATOMS_HOST, ATOMS_PORT))
    try:
        # say hello and ask server for some Atoms (status 'A')
        # identify ourselves with an ID number, formatted as an 8-byte ascii string
        id_str = ('A%7d' % CLIENT_ID).encode('ascii')
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
        print 'Data length: %d' % data_length

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

    # construct an Atoms object from data string
    at = Atoms(data, format='string')
    print 'Received Atoms with %d atoms' % len(at)

    # do the calculation
    pot.calc(at, args_str=args_str)

    # format the results into a string
    data = at.write('string')
    data_length = ('%8s' % len(data)).encode('ascii')

    # preceed message by its length
    data = data_length + data

    # open a new connection to the server to send back the results
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect((ATOMS_HOST, ATOMS_PORT))
    try:
        # say hello and tell server we've got some results (status 'R')
        # identify ourselves with an ID number, formatted as an 8-byte ascii string
        id_str = ('R%7d' % CLIENT_ID).encode('ascii')
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

    print 'Sent Atoms with %d atoms in %d bytes' % (len(at), len(data))
